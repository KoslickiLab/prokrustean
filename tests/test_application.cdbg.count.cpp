#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <set>
#include "const.cpp"	
#include "naive_impl.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.enhance.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/cdbg.count.hpp"

using namespace std;
using namespace sdsl;

void test_unitig_counting_single(){
    int Lmin = 1;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanEnhancement enhancement(prokrustean);
    enhancement.collect_left_right_extensions=true;
    construct_prokrustean(fm_idx, prokrustean, Lmin, &enhancement);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    vector<int> ks={2, 7, 30, 50};
    for(auto k: ks){
        auto unitig_cnt = count_maximal_unitigs_single_k(k, enhancement);

        NaiveCompactedDeBruijnGraph cdbg;
        cdbg.construct_compacted(seq_texts, k);
        auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();

        assert(unitig_cnt==naive_unitig_cnt);
    }
}


void main_application_unitig_count() {
    test_unitig_counting_single();
}
