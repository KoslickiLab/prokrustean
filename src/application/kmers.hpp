#ifndef APPLICATION_KMER_HPP_
#define APPLICATION_KMER_HPP_
#include <algorithm>
#include "../prokrustean.hpp"
#include "../util/data.store.hpp"
#include "../util/string.access.hpp"

using namespace std;


void get_reflectums(int k, ProkrusteanExtension &ext, vector<string> &seq_texts, vector<string> &output){
    auto &prokrustean=ext.prokrustean;

    vector<Edge> spectrum;
    //prokrustean
    for(int i=0; i<prokrustean.sequence_count; i++){
        auto seq=prokrustean.get_sequence(i);
        if(seq.size<k){
            continue;
        }
        prokrustean.get_spectrum(seq, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
            } else {
                output.push_back(seq_texts[i].substr(rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[i].substr(rgn.from, rgn.size()) << endl; 
            }
        }
    }

    for(int i=0; i<prokrustean.stratum_count; i++){
        auto stratum=prokrustean.get_stratum(i);
        if(stratum.size<k){
            continue;
        }
        prokrustean.get_spectrum(stratum, k, spectrum);
        for(auto &rgn: spectrum){
            if(rgn.is_stratified){
                
            } else {
                auto seq_id = ext.stratum_sample_occ_seq_id[i];
                auto pos = ext.stratum_sample_occ_pos[i];
                output.push_back(seq_texts[seq_id].substr(pos+rgn.from, rgn.size()));
                // cout << "from: " << rgn.from << " to: " << rgn.to << " " << seq_texts[seq_id].substr(pos+rgn.from, rgn.size()) << endl; 
            }
        }
    }
}

void get_distinct_kmers(int k, ProkrusteanExtension &ext, vector<string> &seq_texts, vector<string> &output){
    output.clear();
    // vector<string> unitigs;
    get_reflectums(k, ext, seq_texts, output);
    int uniform_cnt=output.size();
    
    // output.reserve(unitigs.size());
    for(int i=0; i<output.size(); i++){
        if(output[i].size()==k){
            // preserve the existing kmer

        } else if(output[i].size()>k) {
            // split the uniform unitig and add from 2nds
            string s=output[i];
            // swtich mer1
            string mer1 = s.substr(0, k);
            output[i]=mer1;
            for(int p=1; p<s.size()-(k-1); p++){
                string mer = s.substr(p, k);
                output.push_back(mer);
            }
        } else {
            // reflectums of degree k cannot be shorter than k
            assert(false);
        }
    }
    // cout << "k:" << k  << " uniform: " << uniform_cnt << " kmers:" << output.size() << endl;
}

void get_distinct_kmers_parallel(int k, ProkrusteanExtension &ext, vector<string> &seq_texts, vector<string> &output){
    output.clear();
    // vector<string> unitigs;
    get_reflectums(k, ext, seq_texts, output);
    int uniform_cnt=output.size();
    
    // output.reserve(unitigs.size());
    for(int i=0; i<output.size(); i++){
        if(output[i].size()==k){
            // preserve the existing kmer

        } else if(output[i].size()>k) {
            // split the uniform unitig and add from 2nds
            string s=output[i];
            // swtich mer1
            string mer1 = s.substr(0, k);
            output[i]=mer1;
            for(int p=1; p<s.size()-(k-1); p++){
                string mer = s.substr(p, k);
                output.push_back(mer);
            }
        } else {
            // reflectums of degree k cannot be shorter than k
            assert(false);
        }
    }
    // cout << "k:" << k  << " uniform: " << uniform_cnt << " kmers:" << output.size() << endl;
}

void store_kmers(const vector<string>& data, const std::string& filename) {
    std::ofstream outputFile(filename);
    
    for(auto &mer: data){
        outputFile << mer << endl;
    }
    
    outputFile.close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_vertex_kmers(int k, Vertex &vertex, vector<Edge> &edges, ProkrusteanExtension &ext, AbstractSequenceAccess &sequence_access, AbstractStringDataStore &string_store){
    ext.prokrustean.get_spectrum(vertex, k, edges);
    if(vertex.is_sequence){
        for(auto &rgn: edges){
            if(rgn.is_reflected){
                auto str=sequence_access.get_substring(vertex.id, rgn.from, rgn.size());
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
                for(int p=0; p<str.size()-(k-1); p++){
                    string mer = str.substr(p, k);
                    string_store.store(mer);
                }
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