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
#include "../src/application/cdbg.hpp"
#include "../src/application/dbg.count.hpp"

using namespace std;
using namespace sdsl;

void test_cdbg_with_verifier(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension enhancement(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<int> ks={4, 8, 12, 16, 40, 80};
    // vector<int> ks={7};
    for(auto k: ks){
    // for(int k=2; k<100; k++){
        CompactedDBGWorkspace workspace;
        extract_paritial_unitigs(k, enhancement, seq_texts, workspace);
        CdbgInvariants veritifer;
        veritifer.set(workspace, enhancement);
        veritifer.assert_result();
        auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
        assert(veritifer.maximal_starting_unitig_count==calculated_unitig_cnt);
    }
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
    while(true){
        Unitig &next_unitig=unitigs[next_id];
        // normal case
        if(next_unitig.is_void_k_minus_1_unitig==false){
            if(next_unitig.is_convergence){
                unitigs[uni_id].nexts.clear();
                unitigs[uni_id].nexts.push_back(next_id);
                break;
            } else {
                if(unitigs[curr_id].nexts.size()<=1){
                    // merge
                    unitigs[uni_id].content+=next_unitig.get_string(sequences).substr(k-1);
                }

                if(next_unitig.nexts.size()!=1){
                    unitigs[uni_id].nexts=next_unitig.nexts;
                    break;
                } else {
                    curr_id=next_id;
                    next_id=next_unitig.nexts[0];
                }
            }
        } else {
            // void case
            if(next_unitig.nexts.size()!=1){
                unitigs[uni_id].nexts=next_unitig.nexts;
                break;
            } else if(next_unitig.nexts.size()==1){
                curr_id=next_id;
                next_id=next_unitig.nexts[0];
            }
        }
    }
}

void construct_cdbg(vector<Unitig> &unitigs, vector<string> &sequences, int k){
    for(UnitigId id=0; id<unitigs.size(); id++){
        auto &unitig = unitigs[id];
        if(unitig.is_start_of_maximal){
            // if(id!=117 && id!=16129)
            // continue;
            // cout << "[max]," << id << " " << unitigs[id].get_string(sequences);
            // if(unitigs[id].is_void_k_minus_1_unitig){
            //     cout << "(void)" << endl;
            // }
            // _explore(id, unitigs, sequences, k);
            // cout << endl;
            build_strings(id, unitigs, sequences, k);
        }
    }
}

void print_cdbg(vector<Unitig> &unitigs){
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

void assert_by_naive_cdbg(vector<Unitig> &unitigs, std::unordered_map<std::string, NaiveNode> &naive_graph){
    int cnt=0;
    for(auto &unitig: unitigs){
        if(unitig.is_start_of_maximal)
        cnt++;
    }
    assert(cnt==naive_graph.size());
    for(UnitigId id=0; id<unitigs.size(); id++){
        auto &unitig=unitigs[id];
        if(unitig.is_start_of_maximal){
            assert(naive_graph.count(unitig.content)>0);
            for(auto next: unitig.nexts){
                auto &next_node=unitigs[next].content;
                assert(naive_graph[unitig.content].outgoing.count(next_node)>0);
            }
        }
    }
}


void test_cdbg_construction(){
    int Lmin = 1;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension enhancement(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    vector<int> ks={2, 7, 30, 50};
    for(auto k: ks){
        CompactedDBGWorkspace workspace;
        extract_paritial_unitigs(k, enhancement, seq_texts, workspace);

        CdbgInvariants veritifer;
        veritifer.set(workspace, enhancement);
        veritifer.assert_result();
        auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
        assert(veritifer.maximal_starting_unitig_count==calculated_unitig_cnt);

        NaiveCompactedDeBruijnGraph cdbg;
        cdbg.construct_compacted(seq_texts, k);
        auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
        // cdbg.print();
        assert(veritifer.maximal_starting_unitig_count==naive_unitig_cnt);

        update_stratum_based_loc_to_seq_based_loc(enhancement, workspace);
        construct_cdbg(workspace.unitigs, seq_texts, k);
        
        assert_by_naive_cdbg(workspace.unitigs, cdbg.graph);
    }
}


void main_application_cdbg() {
    test_cdbg_with_verifier();
    test_cdbg_construction();
}
