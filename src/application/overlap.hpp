#ifndef APPLICATION_OG_HPP_
#define APPLICATION_OG_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"
#include "../util/string.access.hpp"
#include "../util/data.store.hpp"

/* 
*/

struct HOGNode {
    optional<SeqId> node_id_sequence; 
    optional<StratumId> node_id_stratum;
    Pos edge_prefix_length=0;
    optional<StratumId> edge_prefix_id;
    Pos edge_suffix_length=0;
    optional<StratumId> edge_suffix_id;
};

void _store_hog_node(HOGNode &hog_node, AbstractStringDataStore &store){
    string expr;
    if(hog_node.node_id_sequence.has_value()){
        expr=to_string((int)hog_node.node_id_sequence.value());
        expr+=",";
        expr+="T";
    } else {
        expr=to_string((int)hog_node.node_id_stratum.value());
        expr+=",";
        expr+="F";
    }
    expr+=",";
    expr+=to_string((int)hog_node.edge_prefix_length);
    expr+=",";
    expr+=hog_node.edge_prefix_id.has_value()? to_string((int)hog_node.edge_prefix_id.value()) : to_string((int)-1);
    expr+=",";
    expr+=to_string((int)hog_node.edge_suffix_length);
    expr+=",";
    expr+=hog_node.edge_suffix_id.has_value()? to_string((int)hog_node.edge_suffix_id.value()) : to_string((int)-1);
    store.store(expr);
}

void extract_hierarchical_overlap_graph(ProkrusteanExtension &ext, int k, AbstractStringDataStore &store){
    StratifiedEdge stra_edge;
    HOGNode hog_node;
    /* 
    Strata nodes
    */
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]<k){
            continue;
        }
        hog_node.node_id_stratum=i;

        bool has_first=ext.get_stratum_first_stratified(i, stra_edge);
        if(has_first && stra_edge.from==0 && stra_edge.size()>=k){
            hog_node.edge_prefix_id=stra_edge.stratum_id;
            hog_node.edge_prefix_length=ext.prokrustean.stratums__size[i]-stra_edge.size();
        } else {
            hog_node.edge_prefix_id=nullopt;
            hog_node.edge_prefix_length=ext.prokrustean.stratums__size[i];
        }

        bool has_last=ext.get_stratum_last_stratified(i, stra_edge);
        if(has_last && stra_edge.to==ext.prokrustean.stratums__size[i] && stra_edge.size()>=k){
            hog_node.edge_suffix_id=stra_edge.stratum_id;
            hog_node.edge_suffix_length=ext.prokrustean.stratums__size[i]-stra_edge.size();
        } else {
            hog_node.edge_suffix_id=nullopt;
            hog_node.edge_suffix_length=ext.prokrustean.stratums__size[i];
        }

        _store_hog_node(hog_node, store);
    }
    
    hog_node.node_id_stratum=nullopt;

    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        if(ext.prokrustean.sequences__size[i]<k){
            continue;
        }
        hog_node.node_id_sequence=i;

        bool has_first=ext.get_sequence_first_stratified(i, stra_edge);
        if(has_first && stra_edge.from==0 && stra_edge.size()>=k){
            hog_node.edge_prefix_id=stra_edge.stratum_id;
            hog_node.edge_prefix_length=ext.prokrustean.sequences__size[i]-stra_edge.size();
        } else {
            hog_node.edge_prefix_id=nullopt;
            hog_node.edge_prefix_length=0;
        }

        bool has_last=ext.get_sequence_last_stratified(i, stra_edge);
        if(has_last && stra_edge.to==ext.prokrustean.sequences__size[i] && stra_edge.size()>=k){
            hog_node.edge_suffix_id=stra_edge.stratum_id;
            hog_node.edge_suffix_length=ext.prokrustean.sequences__size[i]-stra_edge.size();
        } else {
            hog_node.edge_suffix_id=nullopt;
            hog_node.edge_suffix_length=0;
        }

        _store_hog_node(hog_node, store);
    }
}

#endif