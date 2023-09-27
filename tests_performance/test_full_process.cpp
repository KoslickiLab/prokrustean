#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "util.cpp"	
#include "../src/construction/algorithms.step1_project_stratums.hpp"
#include "../src/construction/algorithms.step2_build_prokrustean.hpp"
#include "../src/prokrustean.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/ssa.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/fm_index/tree_new.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

void test_step1_push(){
    int Lmin=30;
    auto num_threads=12;
    auto sampling_factor=4;

    // auto str=WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    WaveletString str(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
    // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    FmIndex fm_idx(str);
    SampledSuffixArray ssa(fm_idx, sampling_factor);
    Prokrustean prokrustean;
    
    auto start = std::chrono::steady_clock::now();
    vector<future<void>> futures;
    // collect blocks
    auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    // consume blocks
    auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    cout << "finished sampling: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
    SuffixArrayNode_NEW root = get_root_new(fm_idx);
    StratumProjectionOutput output(prokrustean, fm_idx.seq_cnt());
    
    fm_idx.set_sampled_suffix_array(ssa);
    atomic<int> idx_gen;
    start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
    auto func__navigate = [](vector<SuffixArrayNode_NEW> &roots, FmIndex &fm_idx, int Lmin, StratumProjectionOutput &output, atomic<int> &idx_gen) {
        while(true){
            auto idx = idx_gen.fetch_add(1);
            if(idx>=roots.size()){
                break;
            } else {
                navigate_maximals<StratumProjectionOutput, report_representative_locations>(roots[idx], Lmin, fm_idx, output);
            }
        }
    };
    
    // futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));}
    for (auto &f : futures) {f.wait();}
    
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void test_full_process_push(){
    int Lmin=30;
    auto num_threads=12;
    auto sampling_factor=4;

    WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT2_BWT, '$');
    // auto str=WaveletString(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
    // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    FmIndex fm_idx(str);
    SampledSuffixArray ssa(fm_idx, sampling_factor);
    Prokrustean prokrustean;
    
    auto start = std::chrono::steady_clock::now();
    vector<future<void>> futures;
    // collect blocks
    auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    // consume blocks
    auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    cout << "finished sampling: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
    SuffixArrayNode_NEW root = get_root_new(fm_idx);
    StratumProjectionOutput workspace(prokrustean, fm_idx.seq_cnt());
    
    ssa.set_sequences(prokrustean);
    fm_idx.set_sampled_suffix_array(ssa);
    atomic<int> idx_gen;
    start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
    auto func__navigate = [](vector<SuffixArrayNode_NEW> &roots, FmIndex &fm_idx, int Lmin, StratumProjectionOutput &output, atomic<int> &idx_gen) {
        while(true){
            auto idx = idx_gen.fetch_add(1);
            if(idx>=roots.size()){
                break;
            } else {
                navigate_maximals<StratumProjectionOutput, report_representative_locations>(roots[idx], Lmin, fm_idx, output);
            }
        }
    };
    
    // futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(workspace), ref(idx_gen)));}
    for (auto &f : futures) {f.wait();}
    
    cout << "step1 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    workspace.prokrustean->stratums__region.resize(workspace.prokrustean->stratums__size.size());
    workspace.prokrustean->stratums__region_cnt.resize(workspace.prokrustean->stratums__size.size(), 0);

    StratificationWorkSpace step2_workspace;
    for(int i=0;i<fm_idx.seq_cnt(); i++){
        auto &regions = workspace.sequence_regions[i];
        step2_workspace.update_contexts_by_seq(i, regions, prokrustean.stratums__size);
        build_prokrustean(step2_workspace, prokrustean);
    }
    
}

void main_performance_full() {
    // test_step1_push();
    test_full_process_push();
}
