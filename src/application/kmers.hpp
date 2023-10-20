#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "../util/data.store.hpp"
#include "../util/string.access.hpp"

using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_vertex_kmers(int k, Vertex &vertex, vector<Edge> &edges, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store){
    ext.prokrustean.get_spectrum(vertex, k, edges);
    if(vertex.is_sequence){
        for(auto &rgn: edges){
            if(rgn.is_reflected){
                auto str=sequence_access.get_substring(vertex.id, rgn.from, rgn.size());
                // auto str=seq_txts[vertex.id].substr(rgn.from, rgn.size());
                for(int p=0; p<str.size()-(k-1); p++){
                    string mer = str.substr(p, k);
                    string_store.store(mer);
                }
            }
        }
    } else { // stratum
        for(auto &rgn: edges){
            if(rgn.is_reflected){
                auto seq_id = ext.stratum_sample_occ_seq_id[vertex.id];
                auto pos = ext.stratum_sample_occ_pos[vertex.id];
                auto str=sequence_access.get_substring(seq_id, pos+rgn.from, rgn.size());
                // auto str=seq_txts[seq_id].substr(pos+rgn.from, rgn.size());
                for(int p=0; p<str.size()-(k-1); p++){
                    string mer = str.substr(p, k);
                    string_store.store(mer);
                }
            }
        }
    }
}

void get_distinct_kmers_(int k, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store){
    Vertex vertex; 
    vector<Edge> edges;
    // vector<string> unitigs;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        ext.prokrustean.get_sequence(i, vertex);
        if(vertex.size<k){
            continue;
        }
        get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
    }

    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        ext.prokrustean.get_stratum(i, vertex);
        if(vertex.size<k){
            continue;
        }
        get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
    }
}

#endif