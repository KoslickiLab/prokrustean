
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include <algorithm>
#include <numeric>
#include <thread>
#include <future>
#include <list>
#include "algorithms.stage1_projection.hpp"
#include "algorithms.stage2_stratification.hpp"

using namespace std;

struct Configuration{
    // 
    bool collect_ext_count=false;
};

void construct_prokrustean(FmIndex &fm_idx, Prokrustean &prokrustean, uint64_t Lmin=1, ProkrusteanEnhancement* opt=nullptr){
    auto start = std::chrono::steady_clock::now();
    cout << "step1: ";
    SuffixArrayNode root = get_root(fm_idx);
    StratumProjectionWorkspace workspace_step1(prokrustean, fm_idx, opt);
    navigate_strata<StratumProjectionWorkspace, report_representative_locations>(root, Lmin, fm_idx, workspace_step1);
    cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratum_count << endl;

    start = std::chrono::steady_clock::now();
    cout << "step2: ";
    workspace_step1.prepare_prokrustean_spaces();
    StratificationWorkSpace workspace_step2;
    uint64_t cnt=prokrustean.sequence_count;
    for(uint64_t i=0; i<cnt; i++){
        workspace_step2.update_contexts_for_seq(i, fm_idx, workspace_step1, prokrustean);
        build_prokrustean(workspace_step2, prokrustean);
    }
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

auto func__stage1_projection = [](vector<SuffixArrayNode> &roots, FmIndex &fm_idx, int Lmin, StratumProjectionWorkspace &output, atomic<int> &root_idx_gen, int thread_idx) {
    while(true){
        auto idx=root_idx_gen.fetch_add(1);
        if(idx>=roots.size()) 
        break;
        navigate_strata<StratumProjectionWorkspace, report_representative_locations>(roots[idx], Lmin, fm_idx, output, thread_idx);
    }
};

auto func__stage1_flip = [](StratumProjectionWorkspace &output, atomic<int> &block_idx_gen) {
    while(true){
        auto idx=block_idx_gen.fetch_add(1);
        if(idx>=output.block_count) 
        break;
        output.set_block(idx);
    }
};

auto func__stage2_stratifiaction = [](FmIndex &fm_index, Prokrustean &prokrustean, StratumProjectionWorkspace &output, atomic<int> &seq_idx_gen) {
    StratificationWorkSpace workspace;
    while(true){
        auto idx = seq_idx_gen.fetch_add(1);
        if(idx>=prokrustean.sequence_count)
        break;
        workspace.update_contexts_for_seq(idx, fm_index, output, prokrustean);
        build_prokrustean(workspace, prokrustean);
    }
};

void construct_prokrustean_parallel(FmIndex &fm_idx, Prokrustean &prokrustean, int num_threads, int Lmin=10, ProkrusteanEnhancement* opt=nullptr){
    assert(Lmin>=1);
    assert(num_threads>0);
    int root_depth = 5; // is there a clever way to decide the scale of parallelism?
    StratumProjectionWorkspace workspace_step1(prokrustean, fm_idx, opt, 1);
    vector<future<void>> futures;
    atomic<int> root_idx_gen;
    atomic<int> block_idx_gen;
    atomic<int> seq_id_iter;

    cout << "step1: ";
    auto start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode> roots = collect_roots_while_navigate_strata<StratumProjectionWorkspace, report_representative_locations>(Lmin, fm_idx, workspace_step1, root_depth);
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage1_projection, ref(roots), ref(fm_idx), Lmin, ref(workspace_step1), ref(root_idx_gen), 0));
    }
    for (auto &f : futures) {
        f.wait();
    }
    // workspace_step1.dispose();
    cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratum_count << endl;
    cout << "step1.5: ";
    futures.clear();
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage1_flip, ref(workspace_step1), ref(block_idx_gen)));
    }
    for (auto &f : futures) {
        f.wait();
    }
    futures.clear();
    cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratum_count << endl;
    
    // sleep for memory check
    // std::this_thread::sleep_for(std::chrono::seconds(1000));

    start = std::chrono::steady_clock::now();
    cout << "step2: ";
    workspace_step1.prepare_prokrustean_spaces();
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage2_stratifiaction, ref(fm_idx), ref(prokrustean), ref(workspace_step1), ref(seq_id_iter)));
    }
    for (auto &f : futures) {
        f.wait();
    }
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

auto func__recover_texts = [](FmIndex &fm_index, vector<string> &output, atomic<int> &idx_gen) {
    while(true){
        auto idx = idx_gen.fetch_add(1);
        if(idx>=fm_index.seq_cnt()){
            break;
        }
        output[idx]=fm_index.recover_text(idx);
    }
};

void recover_sequences_parallel(FmIndex &fm_idx, vector<string> &sequences, int num_threads){
    sequences.resize(fm_idx.seq_cnt());
    vector<future<void>> futures;
    atomic<int> seq_id_iter;
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__recover_texts, ref(fm_idx), ref(sequences), ref(seq_id_iter)));
    }
    for (auto &f : futures) {
        f.wait();
    }
}
#endif