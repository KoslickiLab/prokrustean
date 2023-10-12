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
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
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
    int k = 32;
    // for(int k=2; k<300; k++){
    //     auto assembled_cnt=extract_compacted_dbg(k, enhancement, seq_texts);
    //     auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
    //     if(assembled_cnt!=calculated_unitig_cnt){
    //         cout << "k neq: " << k << endl;
    //     }
    //     // assert(assembled_cnt==calculated_unitig_cnt);
    // }
    auto assembled_cnt=extract_compacted_dbg(k, enhancement, seq_texts);
    auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
    cout << "assembled_cnt: " << assembled_cnt << endl;
    cout << "computed count: " << calculated_unitig_cnt << endl;
    

    // NaiveCompactedDeBruijnGraph cdbg;
    // cdbg.construct_compacted(seq_texts, k);
    // auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
    // // cdbg.print();
    // cout << "naive count: " << naive_unitig_cnt << endl;
    
    
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
