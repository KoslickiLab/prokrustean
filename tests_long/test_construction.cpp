#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "const.cpp"
#include "naive_impl.cpp"
#include "../src/construction/algorithms.stage1_projection.hpp"
#include "../src/construction/algorithms.stage3_build.hpp"
#include "../src/prokrustean.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/kmers.count.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/ssa.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

void test_full_process_push(){
    int Kmin=1;
    auto num_threads=12;
    int sleep=0;
    auto start = std::chrono::steady_clock::now();

    // WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
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

    
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, Kmin);
    // unordered_map<uint8_t, int> counts;
    // for(auto cnt: prokrustean.stratums__region_cnt){
    //     counts[cnt]++;
    // }
    // for(auto [k,v]: counts){
    //     cout << "cnt: " << (int)k << " amount: " << v << endl;
    // }
    if(sleep>0) std::cout << "sleeping... " << sleep << std::endl;
    fm_idx.dispose();
    prokrustean=Prokrustean();
    std::this_thread::sleep_for(std::chrono::seconds(1000));
    
}


void test_full_process_push_check_kmers(){
    int Kmin=1;
    auto num_threads=12;
    int sleep=0;
    auto start = std::chrono::steady_clock::now();

    // WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    cout << "wavelet string: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    if(sleep>0) std::cout << "1. Wavelete string. sleeping... " << sleep << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(sleep));
    if(sleep>0) std::cout << "2. start " << sleep << std::endl;
    FmIndex fm_idx(str);
    // SampledSuffixArray ssa(fm_idx, sampling_factor);
    Prokrustean prokrustean;
    
    start = std::chrono::steady_clock::now();

    vector<future<void>> futures;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, Kmin);

    start = std::chrono::steady_clock::now();
    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, num_threads);
    cout << "text recovery finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    vector<string> mers;
    start = std::chrono::steady_clock::now();
    vector<int> ks= {1, 4, 10, 15, 20, 30, 42, 47, 66};
    for(auto k: ks){
        int cnt = count_distinct_kmers(k, prokrustean);
        int cnt_naive = count_distinct_kmers_naive(seq_texts, k);
        assert(cnt==cnt_naive);
        cout << "distinct_kmers counted k: " << k << endl;
    }
    cout << "distinct kmers computed: " << (std::chrono::steady_clock::now()-start).count()/1000000 << " microsecond" << endl;
}




void main_construction() {
    // test_full_process_push();
    test_full_process_push_check_kmers();
}
