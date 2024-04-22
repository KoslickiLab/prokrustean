#ifndef APPLICATION_BRAYCURTIS_HPP_
#define APPLICATION_BRAYCURTIS_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"
#include "../util/string.access.hpp"
#include "../util/data.store.hpp"

using namespace std;

/*
Algorithm
For metagenome read set i=1,...,N,
For each metagenome pair (i,j), compute two information.
Bray-Curtis
1) The nominator: min(N_i(mer),N_j(mer)) for mer occurs in both datasets i and j.
2) The denominator: N_i(mer) + N_j(mer) for all mers 

Each stratum S has frequencies f_i(S) meaning the stratum occurs f_i times in the read set i. The sum of such f_i is the total frequency of the stratum.
At a k, when the number of co-occurring k-mers in S is C(S,k),
the nominator is C(S,k)*min(f_i, f_j) for all S, and the denominator is C(S,k)*(f_i + f_j)

For a range of k's, we define dC and C as in counting k-mers, for both nominator and denominator
In the operation of dC and C, each min(f_i, f_j) and (f_i + f_j) is used for computation.

Output matrix: dC and C for each pair. O(|Prokrustean| + |pairs|*|k-range|)
Time complexity: O(|Prokrustean|*|pairs| + |pairs|*|k-range|)
If comparing two datasets and trivially |Prokrustean|>>|k-range|, time and space complexity: O(|Prokrustean|)=O(|m_A|+|m_B|) 
*/
struct Counting {
    // for nominator, (i,j) represents the comparison.
    // for denominator, i represents the index and j is void.
    DatasetId id1; 
    DatasetId id2;
    vector<int64_t> partial_partial_C;
    vector<uint64_t> C;
};

struct BrayCurtisOutput {
    uint64_t k;
    DatasetId id1;
    DatasetId id2;
    uint64_t nominator; // min frequency sum
    uint64_t denominator; // frequency sum
    float value;
};

struct BrayCurtisIntermediate {
    // config
    uint64_t from;
    uint64_t to;

    // dataset info
    vector<DatasetId> dataset_ids_by_sequence;
    DatasetId dataset_count;

    // frequency info
    vector<vector<FrequencyCount>> stratum_frequencies_by_dataset; // frequency_by_stratum_by_dataset[i][x]:= frequencies of stratum x of dataset i
    vector<SpinLock> stratum_locks;

    // counting intermediate values
    // vector<Comparison> nominator_indices;
    vector<Counting> nominator_countings; // each index at x is the dataset pair nominator_indices[x]
    vector<Counting> denominator_countings; // each index is DatasetId

    void initialize(Prokrustean &prokrustean){
        // copy stratum frequency
        stratum_frequencies_by_dataset.resize(dataset_count, vector<FrequencyCount>(prokrustean.stratum_count, 0));
        
        for(DatasetId id=0; id<dataset_count; id++){
            Counting counting;
            counting.id1 = id;
            counting.partial_partial_C.resize(to+1, 0);
            counting.C.resize(to+1, 0);

            denominator_countings.push_back(counting);
        }
        if(dataset_count==1){
            return;
        }
        for (DatasetId id1 = 0; id1 < dataset_count-1; ++id1) {
            for (DatasetId id2 = id1 + 1; id2 < dataset_count; ++id2) {
                Counting counting;
                counting.id1=id1;
                counting.id2=id2;
                counting.partial_partial_C.resize(to+1, 0);
                counting.C.resize(to+1, 0);

                nominator_countings.push_back(counting);
            }
        }
    }
    
    void set_C(StratumId stratum_id, uint64_t cooccurring_cnt, uint64_t k){
        for(auto &counting: nominator_countings){
            counting.C[k]+=cooccurring_cnt * min(stratum_frequencies_by_dataset[counting.id1][stratum_id], stratum_frequencies_by_dataset[counting.id2][stratum_id]);
            // counting.C[k]+=cooccurring_cnt;
        }
        for(auto &counting: denominator_countings){
            counting.C[k]+=cooccurring_cnt * stratum_frequencies_by_dataset[counting.id1][stratum_id];
            // counting.C[k]+=cooccurring_cnt;
        }
    }

    void set_partial_partial_C(StratumId stratum_id, int64_t idx, int value){
        // validation
        if(idx > to){
            return;
        }
        for(auto &counting: nominator_countings){
            counting.partial_partial_C[idx]+=value * (int64_t)min(stratum_frequencies_by_dataset[counting.id1][stratum_id], stratum_frequencies_by_dataset[counting.id2][stratum_id]);
        }
        for(auto &counting: denominator_countings){
            counting.partial_partial_C[idx]+=value * (int64_t)stratum_frequencies_by_dataset[counting.id1][stratum_id];
        }
    }

    void complete_countings(){
        vector<int64_t> partial_C(to+1, 0);
        for(auto &counting: nominator_countings){
            for(int k=from+1; k<to+1; k++){
                partial_C[k]=partial_C[k-1] + counting.partial_partial_C[k];
                counting.C[k]=counting.C[k-1] + partial_C[k];
            }
        }
        for(auto &counting: denominator_countings){
            for(int k=from+1; k<to+1; k++){
                partial_C[k]=partial_C[k-1] + counting.partial_partial_C[k];
                counting.C[k]=counting.C[k-1] + partial_C[k];
            }
        }
    }
};

void compute_frequencies_by_datasets(ProkrusteanExtension &ext, BrayCurtisIntermediate &intermediate){
    assert(ext.stratum_incoming_degrees.size()==ext.prokrustean.stratum_count);
    
    Vertex vertex;
    StratumSize overlap;
    stack<StratumId> completed_strata;
    for(SeqId i=0; i<ext.prokrustean.sequence_count; i++){
        ext.prokrustean.get_sequence(i, vertex);
        DatasetId dt_id=intermediate.dataset_ids_by_sequence[i];
        for(CoveringRegionIdx j=0; j<vertex.s_edges.size(); j++){
            auto &edge=vertex.s_edges[j];
            intermediate.stratum_frequencies_by_dataset[dt_id][edge.stratum_id]+=1;
            
            assert(ext.stratum_incoming_degrees[edge.stratum_id]>0);
            ext.stratum_incoming_degrees[edge.stratum_id]--;
            if(ext.stratum_incoming_degrees[edge.stratum_id]==0){
                completed_strata.push(edge.stratum_id);
            }

            overlap=vertex.overlap_length_on_left(j);
            if(overlap<ext.prokrustean.lmin){
                continue;
            }
            // pinpoint the overlap matching stratum
            assert(ext.prokrustean.stratums__size[edge.stratum_id]>overlap);
            StratumId stratum_id_to_be_resolved=_find_leftmost_descendant_of_matching_length(edge.stratum_id, overlap, ext);
            intermediate.stratum_frequencies_by_dataset[dt_id][stratum_id_to_be_resolved]-=1;
        }
    }
    
    StratumId stratum_id;
    while(!completed_strata.empty()){
        stratum_id=completed_strata.top();
        completed_strata.pop();

        ext.prokrustean.get_stratum(stratum_id, vertex);
        for(CoveringRegionIdx j=0; j<vertex.s_edges.size(); j++){
            auto &edge=vertex.s_edges[j];
            for(DatasetId dt_id=0; dt_id<intermediate.dataset_count; dt_id++){
                intermediate.stratum_frequencies_by_dataset[dt_id][edge.stratum_id]+=intermediate.stratum_frequencies_by_dataset[dt_id][stratum_id];
            }

            // incoming check
            assert(ext.stratum_incoming_degrees[edge.stratum_id]>0);
            ext.stratum_incoming_degrees[edge.stratum_id]--;
            if(ext.stratum_incoming_degrees[edge.stratum_id]==0){
                completed_strata.push(edge.stratum_id);
            }

            overlap=vertex.overlap_length_on_left(j);
            if(overlap<ext.prokrustean.lmin){
                continue;
            }
            // pinpoint the overlap matching stratum
            assert(ext.prokrustean.stratums__size[edge.stratum_id]>overlap);
            StratumId stratum_id_to_be_resolved=_find_leftmost_descendant_of_matching_length(edge.stratum_id, overlap, ext);
            for(DatasetId dt_id=0; dt_id<intermediate.dataset_count; dt_id++){
                intermediate.stratum_frequencies_by_dataset[dt_id][stratum_id_to_be_resolved]-=intermediate.stratum_frequencies_by_dataset[dt_id][stratum_id];
            }
        }
    }
}

void _count_k_mers_per_datasets(Vertex &vertex, int k, BrayCurtisIntermediate &intermediate){
    /* _count_k_mers */
    if(vertex.size<k) return;
    assert(vertex.is_stratum);

    int rgn_cnt=vertex.s_edges.size();
    int64_t cnt=0;
    if(rgn_cnt==0){
        cnt=vertex.size-k+1;
    } else {
        for(int i=0; i<rgn_cnt; i++){
            // reflecting at front
            if(i==0 && vertex.s_edges[i].from>0){
                cnt+=vertex.s_edges[i].from;
            }
            // reflecting at the middle (between stratifieds)
            if(i<rgn_cnt-1){
                // important - consider unsigned type for the comparison
                if(vertex.s_edges[i].to >= vertex.s_edges[i+1].from + k-1){
                } else {
                    cnt+=vertex.s_edges[i+1].from-vertex.s_edges[i].to+(k-1);
                }
            }
            // reflecting at back
            if(i==rgn_cnt-1 && vertex.s_edges[i].to<vertex.size){
                cnt+=vertex.size-vertex.s_edges[i].to;
            }
        }
    }
    intermediate.set_C(vertex.id, cnt, k);
}

void _reflecting_contributions_when_no_stratified_exists(Vertex &vertex, int from, BrayCurtisIntermediate &intermediate){
    if(vertex.size<from){
        return;
    }
    if(vertex.s_edges.size()==0){
        intermediate.set_partial_partial_C(vertex.id, from+1, -1);
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
        intermediate.set_partial_partial_C(vertex.id, k_mer_starts_decreasing_at, -1);
    }
    intermediate.set_partial_partial_C(vertex.id, vertex.size+2, +1);
}

void _reflecting_contributions_between_two_strata(Vertex &vertex, int from, BrayCurtisIntermediate &intermediate){
    if(vertex.size<from){
        return;
    }
    int rgn_cnt=vertex.s_edges.size();
    if(rgn_cnt<2){
        return;
    }
    int maximum_length_on_right=0;
    for(int i=rgn_cnt-2; i>=0; i--){
        if(maximum_length_on_right<vertex.s_edges[i+1].size()){
            maximum_length_on_right=vertex.s_edges[i+1].size();
        }
        // find when a reflecting region starts growing in the edge i 
        int reflecting_appears_at=0;
        if(vertex.s_edges[i].to<=vertex.s_edges[i+1].from){
            // no intersection
            reflecting_appears_at=from+1;
        } else {
            // intersection
            reflecting_appears_at=vertex.s_edges[i].to - vertex.s_edges[i+1].from + 2;
            if(reflecting_appears_at<from + 1){
                reflecting_appears_at=from + 1;
            }
        }
        intermediate.set_partial_partial_C(vertex.id, reflecting_appears_at, +1);
        // find when a reflecting region does not grow or does not start in the edge.
        auto reflecting_stops_increasing_at=2 + std::min<SequenceSize>(maximum_length_on_right, vertex.s_edges[i].size());
        intermediate.set_partial_partial_C(vertex.id, reflecting_stops_increasing_at, -1);
    }
}

void compute_braycurtis_k_range(uint64_t from, uint64_t to, ProkrusteanExtension &ext, vector<DatasetId> &dataset_ids, DatasetId dataset_count, vector<BrayCurtisOutput> &outputs){
    // Definitions:
    // - Def partial partial C(k): The contribution to the change of the changing rate of distinct k-mers count.

    assert(from>0 && from<=to);
    assert(ext.stratum_incoming_degrees.size()==ext.prokrustean.stratum_count);
    
    BrayCurtisIntermediate intermediate;
    intermediate.from=from;
    intermediate.to=to;
    intermediate.dataset_count=dataset_count;
    intermediate.dataset_ids_by_sequence=dataset_ids;
    intermediate.initialize(ext.prokrustean);
    
    compute_frequencies_by_datasets(ext, intermediate);
    
    Vertex vertex;
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        ext.prokrustean.get_stratum(i, vertex, from);
        _count_k_mers_per_datasets(vertex, from, intermediate);
        _reflecting_contributions_when_no_stratified_exists(vertex, from, intermediate);
        _reflecting_contributions_between_two_strata(vertex, from, intermediate);
    }
    
    intermediate.complete_countings();
    for(int k=from; k<to+1; k++){
        for(auto counting: intermediate.nominator_countings){
            BrayCurtisOutput output;
            output.k=k;
            output.id1=counting.id1;
            output.id2=counting.id2;
            output.nominator=counting.C[k];
            output.denominator=intermediate.denominator_countings[counting.id1].C[k]+intermediate.denominator_countings[counting.id2].C[k];
            output.value=2*(float64_t)output.nominator/(float64_t)output.denominator;
            outputs.push_back(output);
        }
    }
}

// void count_distinct_kmers_of_range_parallel(uint64_t from, uint64_t to, int thread_cnt, Prokrustean &prokrustean, vector<uint64_t> &output){
//     // Definitions:
//     // - Def partial partial C(k): The contribution to the change of the changing rate of distinct k-mers count.

//     assert(from>0 && from<=to);
//     vector<vector<uint64_t>> outputs(thread_cnt);
//     vector<vector<uint64_t>> partial_partial_Cs(thread_cnt);
//     for(int i=0; i<thread_cnt; i++){
//         outputs[i].resize(to+1,0);
//         partial_partial_Cs[i].resize(to+1,0);
//     }

//     vector<future<void>> futures;
//     atomic<int> seq_idx_gen;
//     atomic<int> stratum_idx_gen;
//     auto func_ = [](Prokrustean &prokrustean, uint64_t from, uint8_t thread_idx, uint8_t thread_cnt, atomic<int> &seq_idx_gen, atomic<int> &stratum_idx_gen, vector<uint64_t> &output, vector<uint64_t> &partial_partial_C) {
//         Vertex vertex;
//         for(uint64_t i=thread_idx; i<prokrustean.sequence_count; i+=thread_cnt){
//             prokrustean.get_sequence(i, vertex, from);
//             output[from]+=_count_k_mers(vertex, from);
//             _reflecting_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
//             _reflecting_contributions_between_two_strata(vertex, from, partial_partial_C);
//         }
//         for(uint64_t i=thread_idx; i<prokrustean.stratum_count; i+=thread_cnt){
//             prokrustean.get_stratum(i, vertex, from);
//             output[from]+=_count_k_mers(vertex, from);
//             _reflecting_contributions_when_no_stratified_exists(vertex, from, partial_partial_C);
//             _reflecting_contributions_between_two_strata(vertex, from, partial_partial_C);
//         }
//     };
//     for(int i=0; i<thread_cnt; i++){futures.push_back(
//         std::async(std::launch::async, func_, ref(prokrustean), from, i, thread_cnt, ref(seq_idx_gen), ref(stratum_idx_gen), ref(outputs[i]), ref(partial_partial_Cs[i]))
//     );}
//     for (auto &f : futures) {f.wait();}
//     output.clear();
//     output.resize(to+1, 0);
//     vector<uint64_t> partial_partial_C(to+1, 0);
//     vector<uint64_t> dOutput(to+1, 0);
//     for(int i=0; i<thread_cnt; i++){
//         for(int k=0; k<to+1; k++){
//             output[k]+=outputs[i][k];    
//             partial_partial_C[k]+=partial_partial_Cs[i][k];
//         }
//     }
//     for(int k=from+1; k<to+1; k++){
//         dOutput[k]=dOutput[k-1] + partial_partial_C[k];
//         output[k]=output[k-1] + dOutput[k];
//     }
// }
#endif