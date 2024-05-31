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
    int Kmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension enhancement(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Kmin);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<int> ks={4, 8, 12, 16, 40, 80};
    // vector<int> ks={7};
    for(auto k: ks){
    // for(int k=2; k<100; k++){
        CompactedDBGWorkspace workspace;
        cdbg_stage1_extract_paritial_unitigs(k, enhancement, workspace);
        CdbgInvariants veritifer;
        veritifer.set(workspace, enhancement);
        veritifer.assert_result();
        auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, enhancement);
        assert(veritifer.maximal_starting_unitig_count==calculated_unitig_cnt);
    }
}

// void _explore(UnitigId uni_id, vector<Unitig> &unitigs, vector<string> &sequences, int k){
//     int next_cnt=unitigs[uni_id].nexts.size();
//     if(next_cnt>0){
//         cout << "-> ";
//     }
//     for(auto next: unitigs[uni_id].nexts){
//         cout << unitigs[next].get_string(sequences) << " ";
//     }
//     if(unitigs[uni_id].nexts.size()==1){
//         _explore(unitigs[uni_id].nexts[0], unitigs, sequences, k);
//     }
// }


// void print_cdbg(vector<Unitig> &unitigs){
//     for(UnitigId id=0; id<unitigs.size(); id++){
//         auto &unitig = unitigs[id];
//         if(unitig.is_start_of_maximal){
//             cout << unitig.content;
//             if(unitig.nexts.size()>0){
//                 cout << "-> ";
//             }
//             for(auto next: unitig.nexts){
//                 cout << unitigs[next].content << " ";
//             }
//             cout << endl;
//         }
//     }
// }

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
            // for(auto next: unitig.nexts){
            //     auto &next_node=unitigs[next].content;
            for(int i=0; i< unitig.next_cnt; i++){
                auto &next_node=unitigs[unitig.nexts_new[i]].content;
            //     if(naive_graph[unitig.content].outgoing.count(next_node)==0){
            //         cout << "next_node " << next_node << endl; 
            //         for(auto &outgoing: naive_graph[unitig.content].outgoing){
            //             cout << "outgoing " << outgoing << endl; 
            //         }
            //     }
                assert(naive_graph[unitig.content].outgoing.count(next_node)>0);
            }
        }
    }
}


void test_cdbg_construction(){
    int Kmin = 10;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_parallel(fm_idx, prokrustean, 8, Kmin);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, 8);
    
    vector<int> ks={20};
    for(auto k: ks){
        cout << "k: " << k << "start" << endl;
        CompactedDBGWorkspace workspace;
        cdbg_stage1_extract_paritial_unitigs(k, ext, workspace);
        cout << "k: " << k << "extracted" << endl;

        CdbgInvariants veritifer;
        veritifer.set(workspace, ext);
        veritifer.assert_result();
        auto calculated_unitig_cnt=count_maximal_unitigs_single_k(k, ext);
        assert(veritifer.maximal_starting_unitig_count==calculated_unitig_cnt);

        NaiveCompactedDeBruijnGraph cdbg;
        cdbg.construct_compacted(seq_texts, k);
        auto naive_unitig_cnt = cdbg.maximal_unitig_cnt();
        // cdbg.print();
        assert(veritifer.maximal_starting_unitig_count==naive_unitig_cnt);
        cout << "naive constructed" << endl;

        MemorySequenceAccess sequence_access(seq_texts);
        MemoryStringDataStore string_store;
        cdbg_stage2_update_stratum_based_loc_to_seq_based_loc(ext, workspace);
        cdbg_stage3_construct_cdbg(workspace.unitigs, sequence_access, string_store, k);
        cout << "constructed" << endl;
        set<string> unique_str;
        int unique_cnt=0;
        for(auto &unitig: workspace.unitigs){
            if(unitig.is_start_of_maximal){
                unique_str.insert(unitig.content);
                unique_cnt++;   
            }
        }
        assert(unique_str.size()==unique_cnt);
        assert_by_naive_cdbg(workspace.unitigs, cdbg.graph);
    }
}

void test_cdbg_construction_parallel(){
    int Kmin = 5;
    int num_threads=4;
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    prokrustean.contains_stratum_extension_count=true;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, Kmin);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, 8);
    
    vector<int> ks={10, 20};
    for(auto k: ks){
        cout << "k: " << k << "start" << endl;
        CompactedDBGWorkspace workspace1;
        cdbg_stage1_extract_paritial_unitigs_parallel(k, ext, workspace1, num_threads);
        CompactedDBGWorkspace workspace2;
        cdbg_stage1_extract_paritial_unitigs(k, ext, workspace2, num_threads);
        cout << "k: " << k << "extracted" << endl;

        assert(workspace1.unitigs.size()==workspace2.unitigs.size());

        MemorySequenceAccess sequence_access(seq_texts);
        MemoryStringDataStore string_store;

        cdbg_stage2_update_stratum_based_loc_to_seq_based_loc(ext, workspace1);
        cdbg_stage3_construct_cdbg(workspace1.unitigs, sequence_access, string_store, k);

        cdbg_stage2_update_stratum_based_loc_to_seq_based_loc_parallel(ext, workspace2);
        cdbg_stage3_construct_cdbg_parallel(workspace2.unitigs, sequence_access, string_store, k);

        cout << "constructed" << endl;
        set<string> unique_str1;
        set<string> unique_str2;
        for(int i=0; i<workspace1.unitigs.size(); i++){
            if(workspace1.unitigs[i].is_start_of_maximal){
                unique_str1.insert(workspace1.unitigs[i].content);
            }
            if(workspace2.unitigs[i].is_start_of_maximal){
                unique_str2.insert(workspace2.unitigs[i].content);
            }
        }
        assert(unique_str1==unique_str2);
    }
}


void main_application_cdbg() {
    test_cdbg_with_verifier();
    test_cdbg_construction();
    // test_cdbg_construction_parallel();
}
