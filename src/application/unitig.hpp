#ifndef APPLICATION_UNITIG_HPP_
#define APPLICATION_UNITIG_HPP_
#include <algorithm>
#include "../prokrustean.hpp"


void get_uniform_unitigs(int k, Prokrustean &prokrustean, vector<string> &seq_texts, vector<string> &output){
    prokrustean.setup_stratum_example_occ();

    vector<Region> spectrum;
    
    //prokrustean
    for(int i=0; i<prokrustean.sequence_count(); i++){
        auto seq=prokrustean.get_sequence(i);
        prokrustean.get_spectrum(seq, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
            } else {
                output.push_back(seq_texts[i].substr(rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[i].substr(rgn.from, rgn.size()) << endl; 
            }
        }
    }

    for(int i=0; i<prokrustean.stratum_count(); i++){
        auto stratum=prokrustean.get_stratum(i);
        prokrustean.get_spectrum(stratum, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
                
            } else {
                auto seq_id = get<0>(prokrustean.stratum_occ_samples[i]);
                auto pos = get<1>(prokrustean.stratum_occ_samples[i]);
                output.push_back(seq_texts[seq_id].substr(pos+rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[seq_id].substr(pos+rgn.from, rgn.size()) << endl; 
            }
        }
    }
}

#endif