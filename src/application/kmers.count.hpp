#ifndef APPLICATION_KMER_COUNT_HPP_
#define APPLICATION_KMER_COUNT_HPP_
#include <algorithm>
#include "../prokrustean.hpp"

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


void _viable_kmer_decreases(uint64_t size, int from, vector<StratifiedRegion> &regions, vector<int64_t> &partial_partial_C){
    if(size<from){
        return;
    }
    if(from+1 <partial_partial_C.size()){
        // cout << "decrease from " << from+1 << " for size " << size << " having " << regions.size() << endl;
        partial_partial_C[from+1]--;
    }

    if(size+2<partial_partial_C.size()){
        // cout << "no more decrease from " << size+2 << " for size " << size << endl;
        partial_partial_C[size+2]++;
    }
}

void _each_stra_rgn_range_decided_by_intersection(uint64_t size, int k, vector<StratifiedRegion> &regions, vector<int64_t> &partial_partial_C){
    if(size<k){
        return;
    }
    int rgn_cnt=regions.size();
    vector<int> indices;
    indices.reserve(rgn_cnt);
    for(int i=0; i<rgn_cnt; i++){
        if(regions[i].size()>=k){
            indices.push_back(i);
        }
    }

    rgn_cnt=indices.size();
    if(rgn_cnt==0){
        return;
    }

    for(int i=0; i<rgn_cnt; i++){
        int increase_k_from;
        if(i<rgn_cnt-1){
            increase_k_from=max<Pos>(2,regions[indices[i]].to-regions[indices[i+1]].from+2);
        } else {
            increase_k_from=2;
        }
        if(increase_k_from<k+1){
           increase_k_from=k+1; 
        }
        int increase_k_to=regions[indices[i]].size()+2;

        if(increase_k_from<partial_partial_C.size()){
            // cout << "increase k from " << increase_k_from << " for size " << size << " range i " << indices[i] << endl;
            partial_partial_C[increase_k_from]++;
            if(increase_k_to<partial_partial_C.size()){
                partial_partial_C[increase_k_to]--;
            }
        }
    }
}


uint64_t count_distinct_kmers(uint64_t k, Prokrustean &prokrustean){
    uint64_t cnt=0;
    for(int i=0; i<prokrustean.sequence_count; i++){
        auto sequence=prokrustean.get_sequence(i);
        cnt+=_count_k_mers(sequence.size, sequence.s_edges, k);
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        auto stratum=prokrustean.get_stratum(i);
        cnt+=_count_k_mers(stratum.size, stratum.s_edges, k);
    }
    return cnt;
}

void count_distinct_kmers_of_range(uint64_t from, uint64_t to, Prokrustean &prokrustean, vector<uint64_t> &output){
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
    assert(from>0 && from<=to);
    output.clear();
    output.resize(to+1, 0);

    vector<int64_t> partial_C(to+1, 0);
    vector<int64_t> partial_partial_C(to+1, 0);
    vector<int64_t> dOutput(to+1, 0);

    for(int i=0; i<prokrustean.sequence_count; i++){
        auto sequence=prokrustean.get_sequence(i);
        output[from]+=_count_k_mers(sequence.size, sequence.s_edges, from);
        // _stratified_pop_event(sequence.size, sequence.regions, from, partial_C);
        // _during_no_stratified_event(sequence.size, sequence.regions, from, partial_C, partial_partial_C);
        // _during_intersection_exists_event(sequence.size, sequence.regions, from, partial_partial_C);
        _viable_kmer_decreases(sequence.size, from, sequence.s_edges, partial_partial_C);
        _each_stra_rgn_range_decided_by_intersection(sequence.size, from, sequence.s_edges, partial_partial_C);
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        auto stratum=prokrustean.get_stratum(i);
        output[from]+=_count_k_mers(stratum.size, stratum.s_edges, from);
        // _stratified_pop_event(stratum.size, stratum.regions, from, partial_C);
        // _during_no_stratified_event(stratum.size, stratum.regions, from, partial_C, partial_partial_C);
        // _during_intersection_exists_event(stratum.size, stratum.regions, from, partial_partial_C);
        _viable_kmer_decreases(stratum.size, from, stratum.s_edges, partial_partial_C);
        _each_stra_rgn_range_decided_by_intersection(stratum.size, from, stratum.s_edges, partial_partial_C);
    }
    for(int k=from+1; k<to+1; k++){
        dOutput[k]=dOutput[k-1] + partial_partial_C[k];
        output[k]=output[k-1] + dOutput[k];
    }
}

#endif