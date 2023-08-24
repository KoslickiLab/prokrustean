
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include "algorithms.procedures.hpp"
#include <algorithm>
#include <numeric>
#include <thread>
#include <future>
#include <list>

using namespace std;

Prokrustean build_prokrustean(FmIndex &fm_idx, uint64_t Lmin=1, bool recover_sequences=false){
    // step1: collect representative suffix array
    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats;
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    // print_repeats(repeats);
    cout << "repeat(" << repeats.size() << ") collection completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
    
    // step2 get repr structure
    start = std::chrono::steady_clock::now();
    auto repr_annotation = ReprSuffixAnnotation();
    repr_annotation.initialize_rank(fm_idx.size(), repeats);
    repr_annotation.initialize_repr_sa(repeats);
    cout << "repr conversion completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;

    // step3 build structure
    // parallelize by sequences
    start = std::chrono::steady_clock::now();
    Prokrustean pk;
    pk.set_sizes(fm_idx.seq_cnt(), repeats.size(), recover_sequences);
    for(uint64_t i=0; i < fm_idx.seq_cnt(); i++){
        SequenceAnnotation annot = get_sequence_annotations(i, fm_idx, repr_annotation, repeats, recover_sequences);
        // print_positions(annot.positions);
        vector<MinCover> mcs = get_min_covers(annot);
        for(auto mc: mcs){
            if(mc.is_rep){
                pk.rep_mcs[mc.id] = mc;
            } else {
                pk.seq_mcs[mc.id] = mc;
            }
        }
        if(recover_sequences)
        pk.sequences.value()[i] = annot.sequence.value();
    }
    cout << "pk construction completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
    return pk;
}

Prokrustean build_prokrustean_parallel(FmIndex &fm_idx, unsigned int num_threads, uint64_t Lmin=1, bool recover_sequences=false){
    // step1: collect representative suffix array
    int subtree_depth = 3+(int)(num_threads/10);
    subtree_depth = subtree_depth > Lmin? Lmin :subtree_depth;
    auto subtree_roots = collect_nodes(get_root(fm_idx), fm_idx, subtree_depth);
    std::mutex m;
    vector<MaximalRepeatAnnotation> repeats;
    auto func__collect_repeats = [](vector<SuffixArrayNode> &_subtree_roots, uint64_t _Lmin, FmIndex &_fm_idx, std::mutex &_m) {
            vector<MaximalRepeatAnnotation> repeats;
            while(true){
                _m.lock();
                SuffixArrayNode root = _subtree_roots[_subtree_roots.size()-1];
                _subtree_roots.pop_back();
                _m.unlock();

                navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, _Lmin, _fm_idx, repeats);
                if(_subtree_roots.size()==0){
                    break;
                }
            }
            return repeats;
        };

    auto start = std::chrono::steady_clock::now();
    vector<future<vector<MaximalRepeatAnnotation>>> futures;
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__collect_repeats, ref(subtree_roots), Lmin, ref(fm_idx), ref(m)));
    }
    for(int i=0; i<num_threads; i++){
        auto r = futures[i].get();
        repeats.insert(repeats.end(), r.begin(), r.end());
    }
    futures.clear();
    cout << "repeat(" << repeats.size() << ") collection completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
    
    // step2 get repr structure
    start = std::chrono::steady_clock::now();
    auto repr_annotation = ReprSuffixAnnotation();
    repr_annotation.initialize_rank(fm_idx.size(), repeats);
    repr_annotation.initialize_repr_sa(repeats);
    cout << "repr conversion completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;

    // step3 build structure
    auto func__collect_min_covers = //lambda function for distributing tasks to threads
    []( vector<SeqId> &seq_ids,
        Prokrustean &pk,  
        FmIndex &fm_idx, 
        ReprSuffixAnnotation &repr_sa, 
        vector<MaximalRepeatAnnotation> &rep_annots, 
        bool recover_sequences,
        mutex &m) {
            while(true){
                m.lock();
                auto i = seq_ids[seq_ids.size()-1];
                seq_ids.pop_back();
                m.unlock();

                /* parallelized body */
                SequenceAnnotation annot = get_sequence_annotations(i, fm_idx, repr_sa, rep_annots, recover_sequences);
                // print_positions(annot.positions);
                vector<MinCover> mcs = get_min_covers(annot);
                for(auto mc: mcs){
                    /* this part cannot have any collision due to primary repr sa structure*/
                    if(mc.is_rep) pk.rep_mcs[mc.id] = mc;
                    else pk.seq_mcs[mc.id] = mc;
                }
                if(recover_sequences)
                pk.sequences.value()[i] = annot.sequence.value();
                /* parallelized body end */

                if(seq_ids.size()==0){
                    break;
                }
            }
        };

    start = std::chrono::steady_clock::now();
    Prokrustean pk;
    pk.set_sizes(fm_idx.seq_cnt(), repeats.size(), recover_sequences);
    vector<SeqId> seq_ids(fm_idx.seq_cnt());
    iota(seq_ids.begin(), seq_ids.end(), 0);
    // run threads
    vector<future<void>> fut3;
    for(int i=0; i<num_threads; i++){
        fut3.push_back(std::async(std::launch::async, func__collect_min_covers, 
        ref(seq_ids), ref(pk), ref(fm_idx), ref(repr_annotation), ref(repeats), recover_sequences, ref(m)));
    }
    for(int i=0; i<num_threads; i++){
        fut3[i].wait();
    }

    cout << "pk construction completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
    return pk;
}
#endif