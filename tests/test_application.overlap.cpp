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
#include "../src/application/overlap.hpp"

using namespace std;
using namespace sdsl;


void test_simple_overlap_graph_coverage_allowing_multi_edge(){
    int Lmin = 3;
    int min_length = 5;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seqs; 
    fm_idx.recover_all_texts(seqs);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    vector<uint32_t> degrees_in;
    vector<uint32_t> degrees_out;
    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    count_overlap_degrees(ext, min_length, ext.stratum_incoming_degrees, degrees_in, degrees_out, true);

    vector<uint32_t> degrees_in_naive;
    vector<uint32_t> degrees_out_naive;
    count_overlap_degrees_naive(seqs, min_length, degrees_in_naive, degrees_out_naive, true);
    
    // verification
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        // assert(degrees_in[i]+degrees_out[i] == outputs_naive[i]);
        assert(degrees_in[i] == degrees_in_naive[i]);
        assert(degrees_out[i] == degrees_out_naive[i]);
    }
}

void test_simple_overlap_graph_coverage_single_edge(){
    // the implementation is not precise, so Lmin=3 will make this test fail. However, in most datasets this will be enough.
    int Lmin = 3;
    int min_length = 4;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seqs; 
    fm_idx.recover_all_texts(seqs);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    vector<uint32_t> degrees_in;
    vector<uint32_t> degrees_out;
    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    count_overlap_degrees(ext, min_length, ext.stratum_incoming_degrees, degrees_in, degrees_out, false);

    vector<uint32_t> degrees_in_naive;
    vector<uint32_t> degrees_out_naive;
    count_overlap_degrees_naive(seqs, min_length, degrees_in_naive, degrees_out_naive, false);
    
    // verification
    // int cnt=0; 
    // int cnt_naive=0;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        // assert(degrees_in[i]+degrees_out[i] == outputs_naive[i]);
        assert(degrees_in[i] == degrees_in_naive[i]);
        assert(degrees_out[i] == degrees_out_naive[i]);
        // cnt+=(degrees_in[i]+degrees_out[i]);
        // cnt_naive+=(degrees_in_naive[i]+degrees_out_naive[i]);
    }
    // cout << "cnt " << cnt << " cnt_ naive " << cnt_naive << endl;
}

void main_application_overlap(){
    test_simple_overlap_graph_coverage_allowing_multi_edge();
    test_simple_overlap_graph_coverage_single_edge();
}