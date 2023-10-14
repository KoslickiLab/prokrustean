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
    for(auto k: ks){
    // for(int k=2; k<100; k++){
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

void _explore(UnitigId uni_id, vector<Unitig> &unitigs, vector<string> &sequences, int k){
    int next_cnt=unitigs[uni_id].nexts.size();
    if(next_cnt>0){
        cout << "-> ";
    }
    for(auto next: unitigs[uni_id].nexts){
        cout << unitigs[next].get_string(sequences) << " ";
    }
    if(unitigs[uni_id].nexts.size()==1){
        _explore(unitigs[uni_id].nexts[0], unitigs, sequences, k);
    }
}

void build_strings(UnitigId uni_id, vector<Unitig> &unitigs, vector<string> &sequences, int k){
    assert(unitigs[uni_id].is_start_of_maximal);
    // implement string
    unitigs[uni_id].content=unitigs[uni_id].get_string(sequences);
    if(unitigs[uni_id].nexts.size()!=1){
        return;
    }
    UnitigId curr_id=uni_id;
    UnitigId next_id=unitigs[curr_id].nexts[0];
    Unitig &next_unitig=unitigs[next_id];
    while(true){
        next_unitig=unitigs[next_id];
        // normal case
        if(next_unitig.is_void_k_minus_1_unitig==false){
            if(next_unitig.is_convergence){
                unitigs[uni_id].nexts.clear();
                unitigs[uni_id].nexts.push_back(next_id);
                break;
            } else {
                // merge
                unitigs[uni_id].content+=next_unitig.get_string(sequences).substr(k-1);
                if(next_unitig.is_start_of_maximal || next_unitig.nexts.size()==0){
                    unitigs[uni_id].nexts=next_unitig.nexts;
                    break;
                } else {
                    curr_id=next_id;
                    next_id=next_unitig.nexts[0];
                }
            }
        } else {
            // void case
            if(unitigs[curr_id].is_void_k_minus_1_unitig && next_unitig.is_void_k_minus_1_unitig){
                // merge only one letter void-void
                unitigs[uni_id].content+=next_unitig.get_string(sequences).substr(k-2);
            }

            if(next_unitig.is_start_of_maximal || next_unitig.nexts.size()==0){
                unitigs[uni_id].nexts=next_unitig.nexts;
                break;
            } else {
                curr_id=next_id;
                next_id=next_unitig.nexts[0];
            }
        }
    }
}

void construct_cdbg(vector<Unitig> &unitigs, vector<string> &sequences, int k){
    int idx=0;
    cout << "----- construct ------" << endl;
    for(auto& unitig: unitigs){
        if(unitig.is_start_of_maximal){
            auto uni_id = idx;
            // cout << "[max]," << idx << " " << unitigs[idx].get_string(sequences, k);
            // _explore(uni_id, unitigs, sequences, k);
            // cout << endl;
            build_strings(uni_id, unitigs, sequences, k);
        }
        // unitig.print(sequences);
        idx++;
    }
    for(UnitigId id=0; id<unitigs.size(); id++){
        auto &unitig = unitigs[id];
        if(unitig.is_start_of_maximal){
            cout << unitig.content;
            if(unitig.nexts.size()>0){
                cout << "-> ";
            }
            for(auto next: unitig.nexts){
                cout << unitigs[next].content << " ";
            }
            cout << endl;
        }
    }
}

void test_cdbg_construction(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanEnhancement enhancement(prokrustean);
    enhancement.collect_left_right_extensions=true;
    construct_prokrustean(fm_idx, prokrustean, Lmin, &enhancement);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    int k=4;
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
    cdbg.print();
    assert(veritifer.maximal_starting_unitig_count==naive_unitig_cnt);

    update_stratum_based_loc_to_seq_based_loc(enhancement, workspace);
    construct_cdbg(workspace.unitigs, seq_texts, k);
}


void main_application_cdbg() {
    // test_cdbg_with_verifier();
    test_cdbg_construction();
}
