#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "unitig.hpp"
#include "util.hpp"

using namespace std;

void get_distinct_kmers(int k, Prokrustean &prokrustean, vector<string> &seq_texts, vector<string> &output){
    output.clear();
    // vector<string> unitigs;
    get_reflectums(k, prokrustean, seq_texts, output);
    int uniform_cnt=output.size();
    
    // output.reserve(unitigs.size());
    for(int i=0; i<output.size(); i++){
        if(output[i].size()==k){
            // preserve the existing kmer

        } else if(output[i].size()>k) {
            // split the uniform unitig and add from 2nds
            string s=output[i];
            // swtich mer1
            string mer1 = s.substr(0, k);
            output[i]=mer1;
            for(int p=1; p<s.size()-(k-1); p++){
                string mer = s.substr(p, k);
                output.push_back(mer);
            }
        } else {
            // reflectums of degree k cannot be shorter than k
            assert(false);
        }
    }
    // cout << "k:" << k  << " uniform: " << uniform_cnt << " kmers:" << output.size() << endl;
}

#endif