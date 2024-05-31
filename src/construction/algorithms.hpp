
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include <algorithm>
#include <numeric>
#include <thread>
#include <future>
#include <list>
#include "algorithms.stage1_projection.hpp"
#include "algorithms.stage2_annotation.hpp"
#include "algorithms.stage3_build.hpp"

using namespace std;

void construct_prokrustean_single_thread(FmIndex &fm_idx, Prokrustean &prokrustean, uint64_t Kmin=1){
    prokrustean.kmin=Kmin;
    StratumProjectionWorkspace workspace_projection(prokrustean, fm_idx.seq_cnt(), fm_idx.size());
    SuffixAnnotationWorkspace workspace_annotation(fm_idx.size());
    StratificationWorkSpace workspace_stratification;
    
    auto start = std::chrono::steady_clock::now();
    cout << "stage1 collect strata ... ";
    
    /* stage1 - strata are collected. stratified regions are projected as raw data */
    SuffixArrayNode root = get_root(fm_idx);
    navigate_strata<StratumProjectionWorkspace, report_representative_locations>(root, Kmin, fm_idx, workspace_projection);
    workspace_projection.prepare_prokrustean_spaces();

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << " stratum " << prokrustean.stratums__size.size() << endl;
    start = std::chrono::steady_clock::now();
    cout << "stage2 setup suffix array annotations ... ";
    
    /* stage2 - raw data collected in stage1 are arrange to occupy less space and be queried in the next stage */
    for(int i=0; i< workspace_projection.block_count; i++){
        workspace_annotation.set_block(i, workspace_projection);
    }

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();
    cout << "step3 build prokrustean ... ";

    /* stage3 - prokrustean is built from the projected regions */
    for(uint64_t i=0; i<prokrustean.sequence_count; i++){
        workspace_stratification.update_contexts_for_seq(i, fm_idx, workspace_annotation, prokrustean);
        build_prokrustean(workspace_stratification, prokrustean);
    }
    prokrustean.total_sequence_region_count=workspace_stratification.total_seq_region_cnt;
    prokrustean.total_strata_region_count=workspace_stratification.total_strata_region_cnt;

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void _strata_projection_parallel(FmIndex &fm_idx, int Kmin, int thread_cnt, int parallel_roots_collection_depth, StratumProjectionWorkspace &workspace){
    vector<future<void>> futures;
    atomic<int> root_idx_gen;
    auto func_ = [](vector<SuffixArrayNode> &roots, FmIndex &fm_idx, int Kmin, StratumProjectionWorkspace &output, atomic<int> &root_idx_gen) {
        while(true){
            auto idx=root_idx_gen.fetch_add(1);
            if(idx>=roots.size()) 
            break;
            navigate_strata<StratumProjectionWorkspace, report_representative_locations>(roots[idx], Kmin, fm_idx, output);
        }
    };
    vector<SuffixArrayNode> roots = collect_roots_while_navigate_strata<StratumProjectionWorkspace, report_representative_locations>(Kmin, fm_idx, workspace, parallel_roots_collection_depth);

    for(int i=0; i<thread_cnt; i++){ futures.push_back(std::async(std::launch::async, func_, ref(roots), ref(fm_idx), Kmin, ref(workspace), ref(root_idx_gen)));}
    for (auto &f : futures) {f.wait();}

    workspace.prepare_prokrustean_spaces();
}

void _suffix_annotation_parallel(StratumProjectionWorkspace &workspace_projection, SuffixAnnotationWorkspace &workspace_annotation, int thread_cnt){
    vector<future<void>> futures;
    atomic<int> block_idx_gen;
    auto func_ = [](StratumProjectionWorkspace &workspace_projection, SuffixAnnotationWorkspace &workspace_annotation, atomic<int> &block_idx_gen) {
        while(true){
            auto idx=block_idx_gen.fetch_add(1);
            if(idx>=workspace_projection.block_count) 
            break;
            workspace_annotation.set_block(idx, workspace_projection);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(std::async(std::launch::async, func_, ref(workspace_projection), ref(workspace_annotation), ref(block_idx_gen)));}
    for (auto &f : futures) {f.wait();}
}

void _build_prokrustean_parallel(FmIndex &fm_index, Prokrustean &prokrustean, SuffixAnnotationWorkspace &workspace_annotation, int thread_cnt){
    vector<future<void>> futures;
    atomic<int> seq_idx_gen;
    SpinLock lock;
    auto func_ = [](FmIndex &fm_index, Prokrustean &prokrustean, SuffixAnnotationWorkspace &workspace_annotation, atomic<int> &seq_idx_gen, SpinLock &lock) {
        StratificationWorkSpace workspace;
        while(true){
            auto idx = seq_idx_gen.fetch_add(1);
            if(idx>=prokrustean.sequence_count)
            break;
            workspace.update_contexts_for_seq(idx, fm_index, workspace_annotation, prokrustean);
            build_prokrustean(workspace, prokrustean);
        }
        lock.lock();
        prokrustean.total_sequence_region_count+=workspace.total_seq_region_cnt;
        prokrustean.total_strata_region_count+=workspace.total_strata_region_cnt;
        lock.unlock();
    };  
    for(int i=0; i<thread_cnt; i++){futures.push_back(std::async(std::launch::async, func_, ref(fm_index), ref(prokrustean), ref(workspace_annotation), ref(seq_idx_gen), ref(lock)));}
    for (auto &f : futures) {f.wait();}
}


void construct_prokrustean_parallel(FmIndex &fm_idx, Prokrustean &prokrustean, int thread_cnt, int Kmin=10){
    assert(Kmin>=1);
    assert(thread_cnt>0);
    int root_depth = 5; // is there a clever way to decide the scale of parallelism?
    prokrustean.kmin=Kmin;
    StratumProjectionWorkspace* workspace_projection = new StratumProjectionWorkspace(prokrustean, fm_idx.seq_cnt(), fm_idx.size());
    // StratumProjectionWorkspace workspace_projection(prokrustean, fm_idx.seq_cnt(), fm_idx.size(), opt);
    SuffixAnnotationWorkspace workspace_annotation(fm_idx.size());
    StratificationWorkSpace workspace_stratification;

    cout << "stage1 collect strata ... ";
    auto start = std::chrono::steady_clock::now();
    
    /* stage1 - strata are collected. stratified regions are projected as raw data */
    _strata_projection_parallel(fm_idx, Kmin, thread_cnt, root_depth, *workspace_projection);

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();
    cout << "stage2 setup suffix array annotations ... ";
    
    /* stage2 - raw data collected in stage1 are arrange to occupy less space and be queried in the next stage */
    _suffix_annotation_parallel(*workspace_projection, workspace_annotation, thread_cnt);
    delete workspace_projection;

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();
    cout << "stage3 build prokrustean ... ";

    /* stage3 - prokrustean is built from the projected regions */
    _build_prokrustean_parallel(fm_idx, prokrustean, workspace_annotation, thread_cnt);

    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

#endif