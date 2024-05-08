#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <set>
#include "const.cpp"	
#include "naive_impl.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.support.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/dbg.count.hpp"

using namespace std;
using namespace sdsl;

void test_left_right_extension_counting(){
    int Lmin = 1;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    count_left_right_character_extensions(ext);
    
    for(int i=0; i<prokrustean.stratum_count; i++){
        assert(ext.stratum_left_ext_count[i]==ext.prokrustean.get_left_cnt(i));
        assert(ext.stratum_right_ext_count[i]==ext.prokrustean.get_right_cnt(i));
    }
}

void test_left_right_extension_counting_parallel(){
    int Lmin = 1;
    int num_threads=4;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    count_left_right_character_extensions(ext);
    ProkrusteanExtension ext2(prokrustean);
    count_left_right_character_extensions_parallel(ext2, num_threads);
    
    for(int i=0; i<prokrustean.stratum_count; i++){
        assert(ext.stratum_left_ext_count[i]==ext2.stratum_left_ext_count[i]);
        assert(ext.stratum_right_ext_count[i]==ext2.stratum_right_ext_count[i]);
    }
}

void test_left_right_storage(){
    Prokrustean prokrustean;
    prokrustean.stratums__character_extensions_cnt.resize(1);
    CharCount left=1;
    CharCount right=5;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));

    left=0;
    right=0;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));
    
    left=3;
    right=1;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));

    left=4;
    right=4;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));
}

void test_stratum_occ_sampling_parallel(){
    int Lmin = 1;
    int num_threads=4;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seq_txts;
    fm_idx.recover_all_texts(seq_txts);

    Prokrustean prokrustean;
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, Lmin);
	
	ProkrusteanExtension ext(prokrustean);
    ProkrusteanExtension ext2(prokrustean);
	setup_stratum_example_occ(ext);
    setup_stratum_example_occ_parallel(ext2, num_threads);
    for(int i=0; i<prokrustean.stratum_count; i++){
        assert(seq_txts[ext.stratum_sample_occ_seq_id[i]].substr(ext.stratum_sample_occ_pos[i], prokrustean.stratums__size[i])
        == seq_txts[ext2.stratum_sample_occ_seq_id[i]].substr(ext2.stratum_sample_occ_pos[i], prokrustean.stratums__size[i]));
    }
}

void test_frequency_computation(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin=Lmin);

    // Compute naive frequencies
    vector<string> seqs; 
    vector<string> stratum_strings;
    fm_idx.recover_all_texts(seqs);
    setup_stratum_example_occ(ext);
    // too slow - use only 1/20 samples
    int factor=20;
    for(StratumId id=0; id<prokrustean.stratum_count; id++){
        if(id%factor==0){
            string stratum_str=seqs[ext.stratum_sample_occ_seq_id[id]].substr(ext.stratum_sample_occ_pos[id], prokrustean.stratums__size[id]);
            stratum_strings.push_back(stratum_str);
        }
    }
    auto naive_frequencies=get_pattern_frequencies_naive(seqs, stratum_strings);

    // Compute frequencies from prokrustean
    vector<uint8_t> incoming_degrees;
    vector<FrequencyCount> frequencies;
    compute_incoming_degrees(prokrustean, incoming_degrees);
    compute_strata_frequencies(ext, incoming_degrees, frequencies);
    
    // Check
    for(StratumId id=0; id<prokrustean.stratum_count; id++){
        if(id%factor==0){
            assert(naive_frequencies[id/factor]==frequencies[id]);
        }
    }
}

void test_compute_incoming_degrees_parallel(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin=Lmin);

    // Compute frequencies from prokrustean
    vector<uint8_t> incoming_degrees;
    compute_incoming_degrees(prokrustean, incoming_degrees);
    
    vector<uint8_t> incoming_degrees_parallel;
    compute_incoming_degrees_parallel(ext, 4, incoming_degrees_parallel);
    
    // Check
    for(StratumId id=0; id<prokrustean.stratum_count; id++){
        assert(incoming_degrees[id]==incoming_degrees_parallel[id]);
    }
}

void test_frequency_computation_parallel(){
    int Lmin = 1;
    int thread_cnt = 4;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin=Lmin);

    // Compute frequencies from prokrustean
    vector<uint8_t> incoming_degrees;
    vector<FrequencyCount> frequencies;
    compute_incoming_degrees(prokrustean, incoming_degrees);
    compute_strata_frequencies(ext, incoming_degrees, frequencies);
    
    vector<FrequencyCount> frequencies_parallel;
    compute_incoming_degrees(prokrustean, incoming_degrees);
    compute_strata_frequencies_parallel(ext, thread_cnt, incoming_degrees, frequencies_parallel);

    // Check
    for(StratumId id=0; id<prokrustean.stratum_count; id++){
        assert(frequencies_parallel[id]==frequencies[id]);
    }
}

void main_prokrustean_support(){
    test_left_right_extension_counting();
    test_left_right_extension_counting_parallel();
    test_left_right_storage();
    test_stratum_occ_sampling_parallel();
    test_frequency_computation();
    test_compute_incoming_degrees_parallel();
    test_frequency_computation_parallel();
}