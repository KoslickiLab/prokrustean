#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "util.hpp"

using namespace std;

vector<string> collect_distinct_kmers(Prokrustean pk, unsigned int k){
    /* don't need to use set to check uniqueness */
    assert(pk.sequences.has_value());
    vector<Occurrence> rep_occs;
    // vector<Occurrence> rep_occs = collect_rep_example_occurrences(pk);
    vector<string> kmers;
    // from seq min covers
    for(int i=0;i<pk.seqs.size();i++){
        auto string = pk.sequences.value()[i];
        auto gaps = get_gaps_seq(pk, i, k);
        for(auto gap: gaps){
            for(auto interval: gap.chop(k)){
                assert(interval.from+interval.size() <= string.size());
                auto mer = string.substr(interval.from, interval.size());
                for(auto m: kmers){
                    assert(mer != m);
                }
                kmers.push_back(mer);
            }
        }
    }
    // from rep min covers
    for(int i=0;i<pk.stratums.size();i++){
        if(pk.stratums[i].size<k){
            continue;
        }
        Occurrence rep_occ = rep_occs[i];
        assert(rep_occ.from+rep_occ.size() <= pk.sequences.value()[rep_occ.seq_id].size());
        auto string = pk.sequences.value()[rep_occ.seq_id].substr(rep_occ.from, rep_occ.size());
        // cout<< string << endl;
        auto gaps = get_gaps_rep(pk, i, k);
        for(auto gap: gaps){
            for(auto interval: gap.chop(k)){
                assert(interval.from+interval.size()<=string.size());
                auto mer = string.substr(interval.from, interval.size());
                for(auto m: kmers){
                    if(mer == m){
                        cout << mer << " of " << string << endl;
                    }
                    assert(mer != m);
                }
                kmers.push_back(mer);
            }
        }
    }
    
    return kmers;
}

// vector<tuple<SeqId, Pos>> collect_distinct_kmer_pos(Prokrustean pk, unsigned int k){

// }

#endif