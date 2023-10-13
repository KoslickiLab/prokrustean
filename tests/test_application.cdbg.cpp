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

void test_cdbg_with_verifier(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanEnhancement enhancement(prokrustean);
    enhancement.collect_left_right_extensions=true;
    construct_prokrustean(fm_idx, prokrustean, Lmin, &enhancement);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<int> ks={4, 8, 12, 16, 40, 80};
    // vector<int> ks={7};
    // for(auto k: ks){
    for(int k=2; k<100; k++){
        CompactedDBGWorkspace workspace;
        extract_paritial_unitigs(k, enhancement, seq_texts, workspace);
        CdbgVerificationQuantity veritifer;
        veritifer.set(workspace, enhancement);
        veritifer.assert_result();
        auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
        if(veritifer.maximal_starting_unitig_count!=calculated_unitig_cnt){
            // cout << "k neq: " << k << " veritifer.maximal_starting_unitig_count: " << veritifer.maximal_starting_unitig_count << " calculated_unitig_cnt: " << calculated_unitig_cnt << endl;
            assert(false);
        } else {
            // cout << " veritifer.maximal_starting_unitig_count: " << veritifer.maximal_starting_unitig_count << " calculated_unitig_cnt: " << calculated_unitig_cnt << endl;
        }
    }

    // NaiveCompactedDeBruijnGraph cdbg;
    // cdbg.construct_compacted(seq_texts, k);
    // auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
    // // cdbg.print();
    // cout << "naive count: " << naive_unitig_cnt << endl;
}

void test_cdbg_construction(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanEnhancement enhancement(prokrustean);
    enhancement.collect_left_right_extensions=true;
    construct_prokrustean(fm_idx, prokrustean, Lmin, &enhancement);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    int k=10;
    CompactedDBGWorkspace workspace;
    extract_paritial_unitigs(k, enhancement, seq_texts, workspace);

    CdbgVerificationQuantity veritifer;
    veritifer.set(workspace, enhancement);
    veritifer.assert_result();
    auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
    assert(veritifer.maximal_starting_unitig_count==calculated_unitig_cnt);

    NaiveCompactedDeBruijnGraph cdbg;
    cdbg.construct_compacted(seq_texts, k);
    auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
    // cdbg.print();
    cout << "naive count: " << naive_unitig_cnt << "veritifer.maximal_starting_unitig_count: " << calculated_unitig_cnt << endl;
    assert(veritifer.maximal_starting_unitig_count==naive_unitig_cnt);
}


void main_application_cdbg() {
    test_cdbg_with_verifier();
    // test_cdbg_construction();
}
