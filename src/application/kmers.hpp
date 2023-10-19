#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"

using namespace std;


void get_reflectums(int k, Prokrustean &prokrustean, vector<string> &seq_texts, vector<string> &output){
    // prokrustean.setup_stratum_example_occ();
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<Edge> spectrum;
    
    //prokrustean
    for(int i=0; i<prokrustean.sequence_count; i++){
        auto seq=prokrustean.get_sequence(i);
        if(seq.size<k){
            continue;
        }
        prokrustean.get_spectrum(seq, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
            } else {
                output.push_back(seq_texts[i].substr(rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[i].substr(rgn.from, rgn.size()) << endl; 
            }
        }
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        auto stratum=prokrustean.get_stratum(i);
        if(stratum.size<k){
            continue;
        }
        prokrustean.get_spectrum(stratum, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
                
            } else {
                auto seq_id = ext.stratum_sample_occ_seq_id[i];
                auto pos = ext.stratum_sample_occ_pos[i];
                output.push_back(seq_texts[seq_id].substr(pos+rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[seq_id].substr(pos+rgn.from, rgn.size()) << endl; 
            }
        }
    }
}

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