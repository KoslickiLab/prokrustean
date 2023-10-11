#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <set>
#include "const.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;

vector<string> get_distinct_kmers_naive(vector<string> sequences, unsigned int k){
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

// if prokrustean if correct, the kmers will be perfectly collected
void test_distinct_kmers(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    // auto str = WaveletString(PATH1_BWT);
    // auto str = WaveletString(PATH2_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    ProkrusteanEnhancement ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    // for(int i=0; i<prokrustean.stratum_count(); i++){
    //     prokrustean.print_stratum(i, seq_texts);
    // }

    vector<string> output;
    for(int k=2; k<10; k++){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        sort(output.begin(), output.end());
        auto output_naive = get_distinct_kmers_naive(seq_texts, k);
        if(output!=output_naive){
            bool has_duplicates = false;
            for (size_t i = 1; i < output.size(); ++i) {
                if (output[i] == output[i-1]) {  // Compare current element with previous element
                    // Print the duplicate element, only if it hasn't been printed before
                    if(i == 1 || output[i] != output[i-2]) {
                        std::cout << "Duplicate element: " << output[i] << std::endl;
                        has_duplicates = true;
                    }
                }
            }
        }
        assert(output==output_naive);
    }
}

// if prokrustean if correct, the kmers will be perfectly collected
void test_counting_distinct_kmers(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    // auto str = WaveletString(PATH1_BWT);
    // auto str = WaveletString(PATH2_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    ProkrusteanEnhancement ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<uint64_t> counts;
    count_distinct_kmers_of_range(1, 150, prokrustean, counts);

    vector<string> output;
    for(int k=1; k<20; k++){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        assert(counts[k]==output.size());
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

void main_application_kmer() {
    test_distinct_kmers();
    // test_counting_distinct_kmers();
}
