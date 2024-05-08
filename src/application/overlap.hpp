#ifndef APPLICATION_OVERLAP_HPP_
#define APPLICATION_OVERLAP_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"
#include "../util/string.access.hpp"
#include "../util/data.store.hpp"

using namespace std;


void count_overlap_degrees(ProkrusteanExtension &ext, uint32_t min_length, vector<uint8_t> &incoming_degrees, vector<uint32_t> &overlap_degrees_in, vector<uint32_t> &overlap_degrees_out, bool allow_multi_edge){
    overlap_degrees_out.resize(ext.prokrustean.sequence_count, 0);
    overlap_degrees_in.resize(ext.prokrustean.sequence_count, 0);
    vector<uint32_t> strata_prefix_occurrences(ext.prokrustean.stratum_count, 0);
    vector<uint32_t> strata_suffix_occurrences(ext.prokrustean.stratum_count, 0);
    vector<bool> strata_containing_symmetric_prefsuff(ext.prokrustean.stratum_count, false);

    /*
    Compute prefix/suffix occurrences for each stratum
    */ 
    Vertex vertex;
    stack<StratumId> completed_strata;
    for(SeqId i=0; i<ext.prokrustean.sequence_count; i++){
        ext.prokrustean.get_sequence(i, vertex);
        for(CoveringRegionIdx j=0; j<vertex.s_edges.size(); j++){
            auto &s_edge=vertex.s_edges[j];
            assert(incoming_degrees[s_edge.stratum_id]>0);
            incoming_degrees[s_edge.stratum_id]--;
            if(s_edge.size()<min_length){
                continue;
            }

            if(incoming_degrees[s_edge.stratum_id]==0){
                completed_strata.push(s_edge.stratum_id);
            }
            if(j==0 && s_edge.from==0){
                strata_prefix_occurrences[s_edge.stratum_id]+=1;
            }
            if(j+1==vertex.s_edges.size() && s_edge.to==vertex.size){
                strata_suffix_occurrences[s_edge.stratum_id]+=1;
            }
        }
    }

    StratumId stratum_id;
    while(!completed_strata.empty()){
        stratum_id=completed_strata.top();
        completed_strata.pop();

        ext.prokrustean.get_stratum(stratum_id, vertex);

        for(CoveringRegionIdx j=0; j<vertex.s_edges.size(); j++){
            auto &s_edge=vertex.s_edges[j];
            assert(incoming_degrees[s_edge.stratum_id]>0);
            incoming_degrees[s_edge.stratum_id]--;
            if(s_edge.size()<min_length){
                continue;
            }

            if(incoming_degrees[s_edge.stratum_id]==0){
                completed_strata.push(s_edge.stratum_id);
            }
            if(j==0 && s_edge.from==0){
                strata_prefix_occurrences[s_edge.stratum_id]+=strata_prefix_occurrences[stratum_id];
            }
            if(j+1==vertex.s_edges.size() && s_edge.to==vertex.size){
                strata_suffix_occurrences[s_edge.stratum_id]+=strata_suffix_occurrences[stratum_id];
            }
        }
    }

    // compute symmetrics
    if(!allow_multi_edge){
        StratifiedEdge front_edge;
        StratifiedEdge back_edge;
        for(StratumId i=0; i<ext.prokrustean.stratum_count; i++){
            if(ext.prokrustean.get_stratum_size(i)<min_length){
                continue;
            }
            ext.prokrustean.get_stratum(i, vertex);
            stratum_id=i;
            if(ext.get_stratum_first_stratified(stratum_id, front_edge) 
                && front_edge.from==0
                && front_edge.size()>=min_length
                && ext.get_stratum_last_stratified(stratum_id, back_edge)
                && back_edge.to==ext.prokrustean.get_stratum_size(stratum_id)
                && back_edge.size()>=min_length){

                if(front_edge.stratum_id==back_edge.stratum_id){
                    strata_containing_symmetric_prefsuff[stratum_id]=true;
                }
            }
        }
    }
    
    StratifiedEdge edge;
    for(SeqId i=0; i<ext.prokrustean.sequence_count; i++){
        if(ext.get_sequence_first_stratified(i, edge) && edge.from==0){
            stratum_id=edge.stratum_id;
            while(true){
                if(edge.size()<min_length){
                    break;
                }

                if(allow_multi_edge){
                    // this means overlap will be computed multiple times for certain pairs
                    overlap_degrees_in[i]+=strata_suffix_occurrences[stratum_id];
                } else if(strata_containing_symmetric_prefsuff[stratum_id]) {
                   // 
                } else {
                    overlap_degrees_in[i]+=strata_suffix_occurrences[stratum_id];
                }

                if(ext.get_stratum_first_stratified(stratum_id, edge) && edge.from==0){
                    stratum_id=edge.stratum_id;
                } else {
                    break;
                }
            }
        }
        
        if(ext.get_sequence_last_stratified(i, edge) && edge.to==ext.prokrustean.sequences__size[i]){
            stratum_id=edge.stratum_id;
            while(true){
                if(edge.size()<min_length){
                    break;
                }

                if(allow_multi_edge){
                    // this means overlap will be computed multiple times for certain pairs
                    overlap_degrees_out[i]+=strata_prefix_occurrences[stratum_id];
                } else if(strata_containing_symmetric_prefsuff[stratum_id]) {
                   // 
                } else {
                    overlap_degrees_out[i]+=strata_prefix_occurrences[stratum_id];
                }
                if(ext.get_stratum_last_stratified(stratum_id, edge) && edge.to==ext.prokrustean.get_stratum_size(stratum_id)){
                    stratum_id=edge.stratum_id;
                } else {
                    break;
                }
            }
        }

        // cancel out incorrect self-loop occurred when the whole region is stratifying
        if(ext.get_sequence_first_stratified(i, edge) && edge.size()==ext.prokrustean.sequences__size[i]){
            overlap_degrees_in[i]-=1;
            overlap_degrees_out[i]-=1;
        }
    }
}

#endif