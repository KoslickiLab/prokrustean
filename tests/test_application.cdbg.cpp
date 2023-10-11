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
#include "../src/application/cdbg.hpp"

using namespace std;
using namespace sdsl;

void test_cdbg_single(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanEnhancement enhancement(prokrustean);
    enhancement.collect_left_right_extensions=true;
    construct_prokrustean(fm_idx, prokrustean, Lmin, &enhancement);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    // for(auto &txt: seq_texts){
    //     cout << txt << endl;
    // }
    int k = 3;
    extract_compacted_dbg(k, enhancement, seq_texts);

    // vector<string> output;
    // get_reflectums(k-1, prokrustean, seq_texts, output);
    // for(auto &s: output){
    //     cout << "reflectum: " << s << endl;
    // }

    NaiveCompactedDeBruijnGraph cdbg;
    cdbg.construct_compacted(seq_texts, k);
    auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
    cdbg.print();
    // vector<int> ks={2, 7, 30, 50};
    // for(auto k: ks){
    //     auto unitig_cnt = count_maximal_unitigs_single_k(k, enhancement);

    //     NaiveCompactedDeBruijnGraph cdbg;
    //     cdbg.construct_compacted(seq_texts, k);
    //     auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
    //     // cout << " unitigs in naive: " << naive_unitig_cnt << endl;

    //     cout << "unitigs: " << unitig_cnt<< endl;
    //     assert(unitig_cnt==naive_unitig_cnt);
    // }
}


void main_application_cdbg() {
    test_cdbg_single();
}
