
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include "algorithms.procedures.hpp"
#include <algorithm>

using namespace std;

Prokrustean build_prokrustean(FmIndex &fm_idx, uint64_t Lmin=1, bool recover_sequences=false){
    // step1: collect representative suffix array
    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats = navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx);
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

#endif