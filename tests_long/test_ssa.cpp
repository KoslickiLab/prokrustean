#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <future>
#include <random>
#include "const.cpp"	
#include "../src/prokrustean.hpp"
// #include "../src/construction/algorithms.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/ssa.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.stage1_tree_navigation.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

void test_sampling_push(){
    auto num_threads=6;
    auto sampling_factor=8;
    // auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx=FmIndex(str);
    auto ssa=SampledSuffixArray(fm_idx, sampling_factor);
    
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_basic: " << "thread:" << num_threads << ", sampling:" << sampling_factor << endl;
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
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    ssa.validate();
}


void test_sampling_still_works_if_factor_exceeds_seq_length(){
    auto num_threads=6;
    auto sampling_factor=500;// only first positions of sequences may be sampled
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx=FmIndex(str);
    auto ssa=SampledSuffixArray(fm_idx, sampling_factor);
    
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_basic: " << "thread:" << num_threads << ", sampling:" << sampling_factor << endl;
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
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    ssa.validate();
}

void test_sampling_works_the_same_for_sampling_factors(){
    int num_threads=6;
    auto str = WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_index = FmIndex(str);
    atomic<uint64_t> idx_generator;
    sdsl::bit_vector is_sampled_bv(str.size());
    auto ssa1=SampledSuffixArray(fm_index, 3);
    auto ssa2=SampledSuffixArray(fm_index, 7);
    auto ssa3=SampledSuffixArray(fm_index, 10);
    vector<SampledSuffixArray*> ssa_list;
    ssa_list.push_back(&ssa1);
    ssa_list.push_back(&ssa2);
    ssa_list.push_back(&ssa3);
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_works_the_same_for_sampling_factors: " << "thread:" << num_threads << endl;

    vector<future<void>> futures;
    for(auto &ssa: ssa_list){
        start = std::chrono::steady_clock::now();
        // collect blocks
        auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
        for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(*ssa)));}
        for (auto &f : futures) {f.wait();}
        futures.clear();
        // consume blocks
        auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
        for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(*ssa)));}
        for (auto &f : futures) {f.wait();}
        futures.clear();
        cout << "one ssa finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    }

    auto func__compare_loc = [](vector<SampledSuffixArray*> &ssa_list, atomic<uint64_t> &idx_generator, uint64_t idx_limit) {
        uint64_t idx=0;
        optional<tuple<SeqId, Pos>> loc=nullopt;
        optional<tuple<SeqId, Pos>> other_loc=nullopt;
        while(true){
            idx=idx_generator.fetch_add(1);
            if(idx>=idx_limit){
                break;
            }
            loc=nullopt;
            other_loc=nullopt;
            for(auto &ssa: ssa_list){
                loc=ssa->get_location(idx);
                if(other_loc.has_value()){
                    assert(loc==other_loc);
                } else{
                    other_loc=loc;
                }
            }
        }
    };

    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__compare_loc, ref(ssa_list), ref(idx_generator), fm_index.size()));}
    for (auto &f : futures) {f.wait();}
    cout << "location comparison finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    futures.clear();
}

void main_performance_ssa() {
    // test_sampling_push();
    // test_sampling_still_works_if_factor_exceeds_seq_length();
    test_sampling_works_the_same_for_sampling_factors();
}
