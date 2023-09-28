#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "unitig.hpp"
#include "util.hpp"

using namespace std;

void get_distinct_kmers(int k, Prokrustean &prokrustean, vector<string> &seq_texts, vector<string> &output){
    output.clear();
    vector<string> unitigs;
    get_uniform_unitigs(k, prokrustean, seq_texts, unitigs);
    
    output.reserve(unitigs.size());
    for(auto &s: unitigs){
        if(s.size()>k){
            for(int i=0; i<s.size()-(k-1); i++){
                string mer = s.substr(i, k);
                output.push_back(mer);
            }
        } else {
            output.push_back(s);
        }
    }
}

#endif