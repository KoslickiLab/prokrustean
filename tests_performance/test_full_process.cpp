#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "util.cpp"	
#include "../src/construction/models.hpp"
#include "../src/construction/algorithms.procedures_new.hpp"
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

void test_suffix_sampling_push(){
    int Lmin=30;
    auto num_threads=12;
    auto sampling_factor=8;

    // auto str=WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    // auto str=WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    auto fm_idx=FmIndex(str);
    auto ssa=SampledSuffixArray(fm_idx, sampling_factor);
    
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
    StratificationOutput output(fm_idx.seq_cnt());
    
    fm_idx.set_sampled_suffix_array(ssa);
    atomic<int> idx_gen;
    start = std::chrono::steady_clock::now();
    vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
    auto func__navigate = [](vector<SuffixArrayNode_NEW> &roots, FmIndex &fm_idx, int Lmin, StratificationOutput &output, atomic<int> &idx_gen) {
        while(true){
            auto idx = idx_gen.fetch_add(1);
            if(idx>=roots.size()){
                break;
            } else {
                navigate_tree_new<StratificationOutput, report_repr_sa>(roots[idx], Lmin, fm_idx, output);
            }
        }
    };
    
    // futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));}
    for (auto &f : futures) {f.wait();}
    
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void main_performance_full() {
    test_suffix_sampling_push();
}
