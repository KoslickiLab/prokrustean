
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

void construct_prokrustean(FmIndex &fm_idx, Prokrustean &prokrustean, uint64_t Lmin=1, ProkrusteanExtension* opt=nullptr){
    auto start = std::chrono::steady_clock::now();
    cout << "step1 collect strata, prepare suffix array annotation: ";
    prokrustean.lmin=Lmin;
    SuffixArrayNode root = get_root(fm_idx);
    StratumProjectionWorkspace workspace_step1(prokrustean, fm_idx, opt);
    navigate_strata<StratumProjectionWorkspace, report_representative_locations>(root, Lmin, fm_idx, workspace_step1);
    for(int i=0; i< workspace_step1.block_count; i++){
        workspace_step1.set_block(i);
    }
    cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratums__size.size() << endl;

    start = std::chrono::steady_clock::now();
    cout << "step2 build prokrustean: ";
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

void construct_prokrustean_parallel(FmIndex &fm_idx, Prokrustean &prokrustean, int num_threads, int Lmin=10, ProkrusteanExtension* opt=nullptr){
    assert(Lmin>=1);
    assert(num_threads>0);
    int root_depth = 5; // is there a clever way to decide the scale of parallelism?
    StratumProjectionWorkspace workspace_step1(prokrustean, fm_idx, opt, 1);
    vector<future<void>> futures;
    atomic<int> root_idx_gen;
    atomic<int> block_idx_gen;
    atomic<int> seq_id_iter;

    cout << "step1 collect strata, prepare suffix array annotation: ";
    auto start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode> roots = collect_roots_while_navigate_strata<StratumProjectionWorkspace, report_representative_locations>(Lmin, fm_idx, workspace_step1, root_depth);
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage1_projection, ref(roots), ref(fm_idx), Lmin, ref(workspace_step1), ref(root_idx_gen), 0));
    }
    for (auto &f : futures) {
        f.wait();
    }
    futures.clear();
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage1_flip, ref(workspace_step1), ref(block_idx_gen)));
    }
    for (auto &f : futures) {
        f.wait();
    }
    futures.clear();
    cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratums__size.size() << endl;

    start = std::chrono::steady_clock::now();
    cout << "step2 build prokrustean: ";
    workspace_step1.prepare_prokrustean_spaces();
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__stage2_stratifiaction, ref(fm_idx), ref(prokrustean), ref(workspace_step1), ref(seq_id_iter)));
    }
    for (auto &f : futures) {
        f.wait();
    }
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

#endif