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

void test_unitig_counting_single(){
    int Lmin = 1;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension enhancement(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    for(int k=2; k< 50; k++){
        if(k%8!=0){
            continue;
        }
        auto unitig_cnt = count_maximal_unitigs_single_k(k, enhancement);
        NaiveCompactedDeBruijnGraph cdbg;
        cdbg.construct_compacted(seq_texts, k);
        auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();

        assert(unitig_cnt==naive_unitig_cnt);
    }
}


void test_unitig_counting_range(){
    int Lmin = 15;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    vector<uint64_t> output;
    count_maximal_unitigs_range_of_k(30, 200, ext, output);
    for(int k=30; k<200; k++){
        auto unitig_cnt = count_maximal_unitigs_single_k(k, ext);
        assert(unitig_cnt==output[k]);
    }
}


void main_application_unitig_count() {
    test_unitig_counting_single();
    test_unitig_counting_range();
}
