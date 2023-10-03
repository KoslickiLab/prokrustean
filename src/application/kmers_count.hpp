#ifndef APPLICATION_KMER_COUNT_HPP_
#define APPLICATION_KMER_COUNT_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "unitig.hpp"
#include "util.hpp"

using namespace std;

int64_t _count_k_mers(uint64_t size, vector<StratifiedRegion> &regions, int k){
    if(size<k) return 0;

    int rgn_cnt=regions.size();
    vector<int> indices;
    indices.reserve(rgn_cnt);
    for(int i=0; i<rgn_cnt; i++){
        if(regions[i].size()>=k) indices.push_back(i);
    }

    rgn_cnt=indices.size();
    if(rgn_cnt==0){
        return size-k+1;
    }

    int64_t cnt=0;    
    for(int i=0; i<rgn_cnt; i++){
        if(i==0 && regions[indices[i]].from>0){
            cnt+=regions[indices[i]].from;
        }

        if(i<rgn_cnt-1){
            if(regions[indices[i]].to - regions[indices[i+1]].from >= k-1){
            } else {
                // regions[i+1].from-regions[i].to+2*(k-1)-k+1
                cnt+=regions[indices[i+1]].from-regions[indices[i]].to+(k-1);
            }
        }

        if(i==rgn_cnt-1 && regions[indices[i]].to<size){
            cnt+=size-regions[indices[i]].to;
        }
    }
    return cnt;
}

void _stratified_pop_event(vector<StratifiedRegion> &regions, int k, vector<int64_t> &partial_C){
    for(int i=0; i<regions.size(); i++){
        if(regions[i].size()<k) continue;

        partial_C[regions[i].size()+1]++;
    }
}

void _during_no_stratified_event(uint64_t size, vector<StratifiedRegion> &regions, int k, vector<int64_t> &partial_partial_C){
    int max_stra_rgn = k;
    for(int i=0; i<regions.size(); i++){
        if(regions[i].size()<k) continue;

        if(max_stra_rgn<regions[i].size()){
            max_stra_rgn=regions[i].size();
        }
    }
    partial_partial_C[max_stra_rgn+1]--;

    if(size+1<partial_partial_C.size()){
        partial_partial_C[size+1]++;
    }
}

void _during_intersection_exists_event(uint64_t size, vector<StratifiedRegion> &regions, int k, vector<int64_t> &partial_partial_C){
    optional<int> prev_i;
    for(int i=0; i<regions.size(); i++){
        if(regions[i].size()<k) continue;

        if(prev_i.has_value() && regions[i].from<regions[prev_i.value()].to){
            int l=regions[prev_i.value()].to-regions[i].from;
            partial_partial_C[l+2]++;
            int m = min(regions[prev_i.value()].size(), regions[i].size());
            partial_partial_C[m+1]++;
        }

        prev_i=i;
    }
}

void count_kmers_of_range(uint64_t from, uint64_t to, Prokrustean &prokrustean, vector<uint64_t> &output){
// Definitions:
// - Def partial C(k): The contribution to the change of distinct k-mers count.
// - Def partial partial C(k): The contribution to the change of the changing rate of distinct k-mers count.

// Iterate through all v in V such that size(v) >= k:
//     - Let stratified regions in v have sizes n_1, n_2,...,n_t. Then, partial C(n_1+1)++, partial C(n_2+1)++, ..., partial C(n_t+1)++.
//     - Let the maximum size m among stratified regions in v, and m = 0 if no stratified region exists. Then, partial partial C(m+1)-- and partial partial C(size(v)+1)++.
//     - Let each overlapping stratified region pair with the intersection of length l and smaller side n, partial partial C(l+2)++ and partial partial C(n+1)--.

// Def C(l) = distinct l-mers, dC(l) = 0

// Iterate k = l+1 to r:
//     - dC(k) = dC(k-1) + partial partial C(k)
//     - C(k) = C(k-1) + partial C(k) + dC(k)

    output.clear();
    output.resize(to+1, 0);

    vector<int64_t> partial_C(to+1, 0);
    vector<int64_t> partial_partial_C(to+1, 0);
    vector<int64_t> dOutput(to+1, 0);

    for(int i=0; i<prokrustean.sequence_count(); i++){
        auto sequence=prokrustean.get_sequence(i);
        output[from]+=_count_k_mers(sequence.size, sequence.regions, from);
        // _stratified_pop_event(sequence.regions, from, partial_C);
        // _during_no_stratified_event(sequence.size, sequence.regions, from, partial_partial_C);
        // _during_intersection_exists_event(sequence.size, sequence.regions, from, partial_partial_C);
    }

    for(int i=0; i<prokrustean.stratum_count(); i++){
        auto stratum=prokrustean.get_stratum(i);
        output[from]+=_count_k_mers(stratum.size, stratum.regions, from);
        // _stratified_pop_event(stratum.regions, from, partial_C);
        // _during_no_stratified_event(stratum.size, stratum.regions, from, partial_partial_C);
        // _during_intersection_exists_event(stratum.size, stratum.regions, from, partial_partial_C);
    }
    // for(int k=from+1; k<to+1; k++){
    //     dOutput[k]=dOutput[k-1] + partial_partial_C[k];
    //     output[k]=output[k-1] + partial_C[k] + dOutput[k];
    // }
}

#endif