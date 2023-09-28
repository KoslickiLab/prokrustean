#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

vector<string> collect_distinct_kmers_naive(vector<string> sequences, unsigned int k){
    set<string> mers;
    for(auto seq: sequences){
        if(seq.size()<k) continue;

        for(int i=0; i<seq.size()-(k-1); i++){
            string mer = seq.substr(i, k);
            mers.insert(mer);
        }
    }
    vector<string> output(mers.begin(), mers.end());
    return output;
}

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
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    prokrustean.setup_stratum_example_occ();

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    // for(int i=0; i<prokrustean.stratum_count(); i++){
    //     prokrustean.print_stratum(i, seq_texts);
    // }
    vector<string> output;
    for(int k=1; k<20; k++){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        sort(output.begin(), output.end());
        
        assert(output==collect_distinct_kmers_naive(seq_texts, k));
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
    test_basic_construction_w_kmers();
    // test_min_cover_algo();
}
