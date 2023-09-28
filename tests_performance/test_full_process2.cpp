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

void test_step1_push2(){
    int Lmin=20;
    auto num_threads=12;
    int sleep = 10;
    // auto sampling_factor=8;

    // auto str=WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    WaveletString str(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
    // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    // WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    FmIndex fm_idx(str);
    // SampledSuffixArray ssa(fm_idx, sampling_factor);
    Prokrustean prokrustean;
    cout << "wavelte "  << endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
    
    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode_NEW root = get_root_new(fm_idx);
    
    StratumProjectionOutput output(prokrustean, fm_idx.seq_cnt(), fm_idx.size());
    
    atomic<int> idx_gen;
    start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
    vector<future<void>> futures;
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
    cout << "cardinality: " << output.get_cardinality() << endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
}

void test_full_process_push(){
    int Lmin=20;
    auto num_threads=12;
    // auto sampling_factor=16;
    int sleep=10;
    auto start = std::chrono::steady_clock::now();

    // WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    auto str=WaveletString(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
    // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    // auto str=WaveletString(PATH3_PERFORMANCE_SREAD_GUT_GRLBWT_BWT, '$');
    cout << "wavelet string: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    if(sleep>0) std::cout << "1. Wavelete string. sleeping... " << sleep << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
    if(sleep>0) std::cout << "2. start " << sleep << std::endl;
    FmIndex fm_idx(str);
    // SampledSuffixArray ssa(fm_idx, sampling_factor);
    Prokrustean prokrustean;
    
    start = std::chrono::steady_clock::now();
    vector<future<void>> futures;

    SuffixArrayNode_NEW root = get_root_new(fm_idx);
    StratumProjectionOutput workspace(prokrustean, fm_idx.seq_cnt(), fm_idx.size());
    
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
    cout << "cardinality: " << workspace.get_cardinality() << endl;
    if(sleep>0) std::cout << "4. navigated sleeping... " << sleep << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
    if(sleep>0) std::cout << "start " << std::endl;
    
    start = std::chrono::steady_clock::now();
    
    workspace.setup_prokrustean();

    auto func__build = [](FmIndex &fm_index, Prokrustean &prokrustean, StratumProjectionOutput &output, atomic<int> &idx_gen) {
        StratificationWorkSpace workspace;
        while(true){
            auto idx = idx_gen.fetch_add(1);
            if(idx>=prokrustean.sequence_count()){
                break;
            }
            workspace.update_contexts_for_seq(idx, fm_index, output, prokrustean.stratums__size);
            build_prokrustean(workspace, prokrustean);
        }
    };
    atomic<int> seq_id_iter;
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__build, ref(fm_idx), ref(prokrustean), ref(workspace), ref(seq_id_iter)));}
    for (auto &f : futures) {f.wait();}
    cout << "step2 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    cout << "cardinality:" << prokrustean.get_cardinality() << endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
}

void main_performance_full2() {
    // test_step1_push2();
    test_full_process_push();
}
