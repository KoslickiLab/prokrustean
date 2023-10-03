#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/kmers_count.hpp"

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
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    prokrustean.setup_stratum_example_occ();

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    // vector<string> output;
    // for(int k=1; k<20; k++){
    //     get_distinct_kmers(k, prokrustean, seq_texts, output);
    //     sort(output.begin(), output.end());
        
    //     assert(output==collect_distinct_kmers_naive(seq_texts, k));
    // }

    // vector<string> output;
    // for(int k=1; k<21; k++){
    //     get_distinct_kmers(k, prokrustean, seq_texts, output);
    //     cout << "k: " << k << " " << output.size() << endl; 
    // }

    vector<uint64_t> counts;
    count_kmers_of_range(10, 20, prokrustean, counts);
    for(int i=0; i<counts.size(); i++){
        cout << counts[i] << endl;
    }
}

void main_construction() {
    test_basic_construction();
    test_basic_construction_w_kmers();
}
