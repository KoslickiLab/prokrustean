#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"

using namespace std;
using namespace sdsl;

vector<string> collect_distinct_kmers_naive(vector<string> sequences, unsigned int k){
    set<string> uniques;
    for(auto seq: sequences){
        for(int i=0; i<seq.size()-k; i++){
            string mer = seq.substr(i, k);
            uniques.insert(mer);
        }
    }
    vector<string> mers = vector(uniques.begin(), uniques.end());
    return mers;
}

void test_basic_construction(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    
    vector<string> seqs;
    fm_idx.recover_all_texts(seqs);
    for(auto &s: seqs){
        cout << s << endl;
    }

    vector<tuple<SeqId, Pos>> stratum_pos;
    prokrustean.get_stratum_example_occ(stratum_pos);
    
    vector<variant<StratifiedRegion, ReflectedRegion>> output;
    //prokrustean
    for(int i=0; i<prokrustean.sequence_count(); i++){
        auto seq=prokrustean.get_sequence(i);
        prokrustean.get_spectrum(seq, 2, output);
        auto string=fm_idx.recover_text(i);
        for(auto rgn: output){
            if(std::holds_alternative<ReflectedRegion>(rgn)){
                // auto r_rgn = std::get<ReflectedRegion>(rgn);
                // cout << "from: " << r_rgn.from << " to: " << r_rgn.to << " " << string.substr(r_rgn.from, r_rgn.size()) << endl; 
            } else {
                auto r_rgn = std::get<StratifiedRegion>(rgn);
                cout << "from: " << r_rgn.from << " to: " << r_rgn.to << " " << string.substr(r_rgn.from, r_rgn.size()) << endl; 
            }
        }
    }
    for(int i=0; i<prokrustean.stratum_count(); i++){
        auto stratum=prokrustean.get_stratum(i);
        prokrustean.get_spectrum(stratum, 2, output);
        for(auto rgn: output){
            if(std::holds_alternative<ReflectedRegion>(rgn)){
                auto r_rgn=std::get<ReflectedRegion>(rgn);
                auto seq_id = get<0>(stratum_pos[i]);
                auto pos = get<1>(stratum_pos[i]);
                auto string=fm_idx.recover_text(seq_id);
                
                cout << "reflected (stratum)" << string.substr(pos, r_rgn.size()) << endl; 
            } else {
                auto s_rgn = std::get<StratifiedRegion>(rgn);
                auto seq_id = get<0>(stratum_pos[i]);
                auto pos = get<1>(stratum_pos[i]);
                auto string=fm_idx.recover_text(seq_id);
                cout << "from: " << s_rgn.from << " to: " << s_rgn.to << " " << string.substr(s_rgn.from, s_rgn.size()) << endl; 
            }
        }
    }
    
    //naive
    for(auto &s:collect_distinct_kmers_naive(seqs, 2)){
        cout << s << endl;
    }
}

// void test_distinct_kmers(){
//     int Lmin = 1;
//     auto str = WaveletString(PATH1_BWT);
//     auto fm_idx = FmIndex(str);
//     Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
//     auto sequences = recover_text(fm_idx);
//     print_prokrustean(pk);

//     for(int k=2; k<7; k++){
//         cout << "k: " << k << endl;
//         vector<string> mers = collect_distinct_kmers(pk, k);
//         vector<string> mers_naive = collect_distinct_kmers_naive(sequences, k);
//         sort(mers.begin(), mers.end());
//         sort(mers_naive.begin(), mers_naive.end());
//         cout << "-- naive --" << endl;
//         for(auto m: mers_naive){
//             cout << m << endl;
//         }
//         cout << "-- pk --" << endl;
//         for(auto m: mers){
//             cout << m << endl;
//         }
//         vector<string>::iterator mers_itr = mers.begin();
//         vector<string>::iterator mers_naive_itr = mers_naive.begin();
//         while(mers_itr<mers.end()){
//             if(*mers_itr != *mers_naive_itr){
//                 cout << "not matched: " << *mers_itr << ", " << *mers_naive_itr << endl;
//             }
//             assert(*mers_itr == *mers_naive_itr);
//             mers_itr++;
//             mers_naive_itr++;
//         }
//         assert(mers_naive.size()==mers.size());
//     }
    
//     // for(auto s: mers_naive){
//     //     cout << s << endl;
//     // }
//     // for(auto mer: collect_distinct_kmer(pk, 2)){
//     //     cout << mer << endl;
//     // }
// }

void main_construction() {
    test_basic_construction();
    // test_distinct_kmers();
    // test_min_cover_algo();
}
