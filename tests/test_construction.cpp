#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/application/kmers.hpp"

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
    Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    // cout << "finished with mcs: " << pk.rep_mcs.size() << endl;
    print_prokrustean(pk);
    print_bare_prokrustean(pk);
    
}

void test_distinct_kmers(){
    int Lmin = 1;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    auto sequences = recover_text(fm_idx);

    for(int k=Lmin; k<7; k++){
        cout << "k: " << k << endl;
        vector<string> mers = collect_distinct_kmers(pk, k);
        vector<string> mers_naive = collect_distinct_kmers_naive(sequences, k);
        sort(mers.begin(), mers.end());
        sort(mers_naive.begin(), mers_naive.end());
        assert(mers_naive.size()==mers.size());
        vector<string>::iterator mers_itr = mers.begin();
        vector<string>::iterator mers_naive_itr = mers_naive.begin();
        while(mers_itr<mers.end()){
            assert(*mers_itr == *mers_naive_itr);
            mers_itr++;
            mers_naive_itr++;
        }
    }
    
    // for(auto s: mers_naive){
    //     cout << s << endl;
    // }
    // for(auto mer: collect_distinct_kmer(pk, 2)){
    //     cout << mer << endl;
    // }
}


void main_construction_mc() {
    test_basic_construction();
    test_distinct_kmers();
}
