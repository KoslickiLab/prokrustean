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

vector<string> collect_kmers_naive(vector<string> sequences, unsigned int k){
    vector<string> mers;
    for(auto seq: sequences){
        for(int i=0; i<seq.size()-k; i++){
            string mer = seq.substr(i, k);
            mers.push_back(mer);
        }
    }
    sort(mers.begin(), mers.end());
    return mers;
}

vector<string> collect_kmers_prokrustean(Prokrustean pk, vector<string> sequences, unsigned int k){
    vector<string> mers;
    auto info = collect_distinct_kmer_pos(pk, k);
    for(auto pair: info){
        SeqId id = get<0>(pair);
        Pos pos = get<1>(pair);
        string mer = sequences[id].substr(pos, k);
        mers.push_back(mer);
    }
    sort(mers.begin(), mers.end());
    return mers;
}


void test_distinct_kmers(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    Prokrustean pk = build_prokrustean(fm_idx, Lmin);

    auto sequences = recover_text(fm_idx);
    for(int k=2; k<6; k++){
        vector<string> mers_naive = collect_kmers_naive(sequences, k);
        vector<string> mers_pk = collect_kmers_prokrustean(pk, sequences, k);
        assert(mers_naive.size()==mers_pk.size());
        for(int i=0; i<mers_pk.size(); i++){
            assert(mers_naive[i]==mers_pk[i]);
        }
    }
}


void main_construction_mc() {
    test_distinct_kmers();
}
