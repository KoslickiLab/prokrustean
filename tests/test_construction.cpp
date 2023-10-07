#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "const.cpp"	
#include "naive_impl.cpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;


void test_basic_construction(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
}

// if prokrustean if correct, the kmers will be perfectly collected
void test_basic_construction_w_kmers(){
    int Lmin = 1;
    auto str = WaveletString(PATH4_SREAD_PARTITIONED);
    // auto str = WaveletString(PATH3_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    prokrustean.setup_stratum_example_occ();

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    vector<string> output;
    vector<int> testing_ks={1,5,10,20};
    for(auto k: testing_ks){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        sort(output.begin(), output.end());
        assert(output==collect_distinct_kmers_naive(seq_texts, k));
    }
}

void main_construction() {
    test_basic_construction();
    test_basic_construction_w_kmers();
}
