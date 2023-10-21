#ifndef APPLICATION_KMER_COUNT_HPP_
#define APPLICATION_KMER_COUNT_HPP_
#include <algorithm>
#include "../prokrustean.hpp"

using namespace std;

int64_t _count_k_mers(Vertex &vertex, int k){
    /* _count_k_mers
    count kmers directly from the vertex.
    It would be straightforward to work with the cover, 
    but for efficiency we compute in-line.
    assumption - no s_edge of length under k is included, 
    which should be cared in the construction of the vertex.
    */
    if(vertex.size<k) return 0;

    int rgn_cnt=vertex.s_edges.size();
    if(rgn_cnt==0){
        return vertex.size-k+1;
    }
    int64_t cnt=0;    
    for(int i=0; i<rgn_cnt; i++){
        // refracted at front
        if(i==0 && vertex.s_edges[i].from>0){
            cnt+=vertex.s_edges[i].from;
        }
        // refracted at the middle (between stratifieds)
        if(i<rgn_cnt-1){
            if(vertex.s_edges[i].to - vertex.s_edges[i+1].from >= k-1){
            } else {
                cnt+=vertex.s_edges[i+1].from-vertex.s_edges[i].to+(k-1);
            }
        }
        // refracted at back
        if(i==rgn_cnt-1 && vertex.s_edges[i].to<vertex.size){
            cnt+=vertex.size-vertex.s_edges[i].to;
        }
    }
    return cnt;
}

void _refracted_contributions_when_no_stratified_exists(Vertex &vertex, int from, vector<uint64_t> &partial_partial_C){
    if(vertex.size<from){
        return;
    }
    if(vertex.s_edges.size()==0){
        if(from+1<partial_partial_C.size()){
            // cout << "decrease from " << from+1 << " for size " << size << " having " << regions.size() << endl;
            partial_partial_C[from+1]--;
        }
    } else {
        auto max_s_length=0;
        for(auto s: vertex.s_edges){
            if(max_s_length < s.size()){
                max_s_length=s.size();
            }
        }
        auto k_mer_starts_decreasing_at=2+max_s_length;
        if(k_mer_starts_decreasing_at<from+1){
            k_mer_starts_decreasing_at=from+1;
        }
        
        if(k_mer_starts_decreasing_at<partial_partial_C.size()){
            // cout << "decrease from " << from+1 << " for size " << size << " having " << regions.size() << endl;
            partial_partial_C[k_mer_starts_decreasing_at]--;
        }
    }
    // kmer stops decreasing. Until vertex.size+1, k mer contributed by the vertex decreases
    if(vertex.size+2<partial_partial_C.size()){
        // cout << "no more decrease from " << size+2 << " for size " << size << endl;
        partial_partial_C[vertex.size+2]++;
    }
}

void _refracted_contributions_between_two_strata(Vertex &vertex, int from, vector<uint64_t> &partial_partial_C){
    if(vertex.size<from){
        return;
    }
    int rgn_cnt=vertex.s_edges.size();
    if(rgn_cnt<2){
        return;
    }
    auto maximum_length_on_right=0;
    for(int i=rgn_cnt-2; i>=0; i--){
        if(maximum_length_on_right<vertex.s_edges[i+1].size()){
            maximum_length_on_right=vertex.s_edges[i+1].size();
        }
        // find when a refracted region starts growing in the edge i 
        auto refracted_appears_at=0;
        if(vertex.s_edges[i].to<=vertex.s_edges[i+1].from){
            // no intersection
            refracted_appears_at=from+1;
        } else {
            // intersection
            refracted_appears_at=vertex.s_edges[i].to - vertex.s_edges[i+1].from + 2;
            if(refracted_appears_at<from + 1){
                refracted_appears_at=from + 1;
            }
        }
        if(refracted_appears_at<partial_partial_C.size()){
            partial_partial_C[refracted_appears_at]++;
        }
        // find when a refracted region does not grow or does not start in the edge.
        auto refracted_stops_increasing_at=2 + std::min<SequenceSize>(maximum_length_on_right, vertex.s_edges[i].size());
        if(refracted_stops_increasing_at<partial_partial_C.size()){
            partial_partial_C[refracted_stops_increasing_at]--;
        }
    }
}

uint64_t count_distinct_kmers(uint64_t k, Prokrustean &prokrustean){
    uint64_t cnt=0;
    for(int i=0; i<prokrustean.sequence_count; i++){
        auto sequence=prokrustean.get_sequence(i, k);
        cnt+=_count_k_mers(sequence, k);
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        auto stratum=prokrustean.get_stratum(i, k);
        cnt+=_count_k_mers(stratum, k);
    }
    return cnt;
}

void count_distinct_kmers_of_range(uint64_t from, uint64_t to, Prokrustean &prokrustean, vector<uint64_t> &output){
    // Definitions:
    // - Def partial partial C(k): The contribution to the change of the changing rate of distinct k-mers count.

    assert(from>0 && from<=to);
    output.clear();
    output.resize(to+1, 0);

    vector<uint64_t> partial_partial_C(to+1, 0);
    vector<uint64_t> dOutput(to+1, 0);
    Vertex vertex;
    for(int i=0; i<prokrustean.sequence_count; i++){
        prokrustean.get_sequence(i, vertex, from);
        output[from]+=_count_k_mers(vertex, from);
        _refracted_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
        _refracted_contributions_between_two_strata(vertex, from, partial_partial_C);
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        prokrustean.get_stratum(i, vertex, from);
        output[from]+=_count_k_mers(vertex, from);
        _refracted_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
        _refracted_contributions_between_two_strata(vertex, from, partial_partial_C);
    }
    for(int k=from+1; k<to+1; k++){
        dOutput[k]=dOutput[k-1] + partial_partial_C[k];
        output[k]=output[k-1] + dOutput[k];
    }
}

void count_distinct_kmers_of_range_parallel(uint64_t from, uint64_t to, int thread_cnt, Prokrustean &prokrustean, vector<uint64_t> &output){
    // Definitions:
    // - Def partial partial C(k): The contribution to the change of the changing rate of distinct k-mers count.

    assert(from>0 && from<=to);
    vector<vector<uint64_t>> outputs(thread_cnt);
    vector<vector<uint64_t>> partial_partial_Cs(thread_cnt);
    for(int i=0; i<thread_cnt; i++){
        outputs[i].resize(to+1,0);
        partial_partial_Cs[i].resize(to+1,0);
    }

    vector<future<void>> futures;
    atomic<int> seq_idx_gen;
    atomic<int> stratum_idx_gen;
    auto func_ = [](Prokrustean &prokrustean, uint64_t from, uint8_t thread_idx, uint8_t thread_cnt, atomic<int> &seq_idx_gen, atomic<int> &stratum_idx_gen, vector<uint64_t> &output, vector<uint64_t> &partial_partial_C) {
        Vertex vertex;
        for(uint64_t i=thread_idx; i<prokrustean.sequence_count; i+=thread_cnt){
            prokrustean.get_sequence(i, vertex, from);
            output[from]+=_count_k_mers(vertex, from);
            _refracted_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
            _refracted_contributions_between_two_strata(vertex, from, partial_partial_C);
        }
        for(uint64_t i=thread_idx; i<prokrustean.stratum_count; i+=thread_cnt){
            prokrustean.get_stratum(i, vertex, from);
            output[from]+=_count_k_mers(vertex, from);
            _refracted_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
            _refracted_contributions_between_two_strata(vertex, from, partial_partial_C);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(
        std::async(std::launch::async, func_, ref(prokrustean), from, i, thread_cnt, ref(seq_idx_gen), ref(stratum_idx_gen), ref(outputs[i]), ref(partial_partial_Cs[i]))
    );}
    for (auto &f : futures) {f.wait();}
    output.clear();
    output.resize(to+1, 0);
    vector<uint64_t> partial_partial_C(to+1, 0);
    vector<uint64_t> dOutput(to+1, 0);
    for(int i=0; i<thread_cnt; i++){
        for(int k=0; k<to+1; k++){
            output[k]+=outputs[i][k];    
            partial_partial_C[k]+=partial_partial_Cs[i][k];
        }
    }
    for(int k=from+1; k<to+1; k++){
        dOutput[k]=dOutput[k-1] + partial_partial_C[k];
        output[k]=output[k-1] + dOutput[k];
    }
}
#endif