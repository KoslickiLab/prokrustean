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
void _get_vertex_kmers(int k, Vertex &vertex, vector<Edge> &edges, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store){
    string str;
    ext.prokrustean.get_spectrum(vertex, k, edges);
    if(vertex.is_sequence){
        for(auto &rgn: edges){
            if(rgn.is_reflected){
               sequence_access.read_seq_substr(vertex.id, rgn.from, rgn.size(), str);
               string_store.chop_and_store(str, k);
                // for(int p=0; p<str.size()-(k-1); p++){
                //     string mer = str.substr(p, k);
                //     string_store.store(mer);
                // }
            }
        }
    } else { // stratum
        for(auto &rgn: edges){
            if(rgn.is_reflected){
                auto seq_id = ext.stratum_sample_occ_seq_id[vertex.id];
                auto pos = ext.stratum_sample_occ_pos[vertex.id];
                sequence_access.read_seq_substr(seq_id, pos+rgn.from, rgn.size(), str);
                string_store.chop_and_store(str, k);
                // for(int p=0; p<str.size()-(k-1); p++){
                //     string mer = str.substr(p, k);
                //     string_store.store(mer);
                // }
            }
        }
    }
}

void get_distinct_kmers(int k, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store){
    Vertex vertex; 
    vector<Edge> edges;
    // vector<string> unitigs;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        ext.prokrustean.get_sequence(i, vertex);
        if(vertex.size<k){
            continue;
        }
        _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
    }

    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        ext.prokrustean.get_stratum(i, vertex);
        if(vertex.size<k){
            continue;
        }
        _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
    }
}

void get_distinct_kmers_parallel(int k, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, DiskStringDataStore &string_store, int thread_cnt){
    vector<future<void>> futures;
    auto func_ = [](ProkrusteanExtension &ext, int k, uint8_t thread_idx, uint8_t thread_cnt, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store) {
        Vertex vertex; 
        vector<Edge> edges;
        // vector<string> unitigs;
        for(int i=thread_idx; i<ext.prokrustean.sequence_count; i+=thread_cnt){
            ext.prokrustean.get_sequence(i, vertex);
            if(vertex.size<k){
                continue;
            }
            _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
        }

        for(int i=thread_idx; i<ext.prokrustean.stratum_count; i+=thread_cnt){
            ext.prokrustean.get_stratum(i, vertex);
            if(vertex.size<k){
                continue;
            }
            _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(
        std::async(std::launch::async, func_, ref(ext), k, i, thread_cnt, ref(sequence_access), ref(string_store))
    );}
    for (auto &f : futures) {f.wait();}
}


void get_distinct_kmers_parallel_multi_file(int k, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, vector<DiskStringDataStore*> &string_store_list, int thread_cnt){
    vector<future<void>> futures;
    auto func_ = [](ProkrusteanExtension &ext, int k, uint8_t thread_idx, uint8_t thread_cnt, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store) {
        Vertex vertex; 
        vector<Edge> edges;
        // vector<string> unitigs;
        for(int i=thread_idx; i<ext.prokrustean.sequence_count; i+=thread_cnt){
            ext.prokrustean.get_sequence(i, vertex);
            if(vertex.size<k){
                continue;
            }
            _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
        }

        for(int i=thread_idx; i<ext.prokrustean.stratum_count; i+=thread_cnt){
            ext.prokrustean.get_stratum(i, vertex);
            if(vertex.size<k){
                continue;
            }
            _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(
        std::async(std::launch::async, func_, ref(ext), k, i, thread_cnt, ref(sequence_access), ref(*string_store_list[i]))
    );}
    for (auto &f : futures) {f.wait();}
}

// void get_distinct_kmers_parallel(int k, ProkrusteanExtension &ext, vector<DiskSequenceAccess> &sequence_access_list, AbstractStringDataStore &string_store, int thread_cnt){
//     vector<future<void>> futures;
//     auto func_ = [](ProkrusteanExtension &ext, int k, uint8_t thread_idx, uint8_t thread_cnt, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store) {
//         Vertex vertex; 
//         vector<Edge> edges;
//         // vector<string> unitigs;
//         for(int i=thread_idx; i<ext.prokrustean.sequence_count; i+=thread_cnt){
//             ext.prokrustean.get_sequence(i, vertex);
//             if(vertex.size<k){
//                 continue;
//             }
//             _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
//         }

//         for(int i=thread_idx; i<ext.prokrustean.stratum_count; i+=thread_cnt){
//             ext.prokrustean.get_stratum(i, vertex);
//             if(vertex.size<k){
//                 continue;
//             }
//             _get_vertex_kmers(k, vertex, edges, ext, sequence_access, string_store);
//         }
//     };
//     for(int i=0; i<thread_cnt; i++){futures.push_back(
//         std::async(std::launch::async, func_, ref(ext), k, i, thread_cnt, ref(sequence_access_list[i]), ref(string_store))
//     );}
//     for (auto &f : futures) {f.wait();}
// }

#endif