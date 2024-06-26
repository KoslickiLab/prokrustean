#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <stack>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include "data_types.hpp"

using namespace std;

// for succinct left/right counts
int HALF_BITS_IN_CHAR = 4*sizeof(CharCount);
CharCount MASK_HALF_BITS = (1 << HALF_BITS_IN_CHAR) - 1;

struct Edge {
    Pos from;
    Pos to;

    bool is_stratified;
    bool is_reflected;
    // stratified only 
    StratumId stratum_id;

    Edge(){}
    Edge(Pos from, Pos to): Edge(from, to, false, 0){}
    Edge(Pos from, Pos to, StratumId stratum_id): Edge(from, to, true, stratum_id){}
    Edge(Pos from, Pos to, bool is_stratified, StratumId stratum_id): from(from), to(to), is_stratified(is_stratified), is_reflected(!is_stratified), stratum_id(stratum_id) {}

    uint64_t size(){
        assert(from<to);
        return to-from;
    }

    void update(Pos from, Pos to, bool is_stratified, StratumId stratum_id){
        this->from=from;
        this->to=to;
        this->is_stratified=is_stratified;
        this->is_reflected=!is_stratified;
        this->stratum_id=stratum_id;
    }

    void print(){
        if(is_stratified){
            cout<<"stratified("<< stratum_id <<")";
        } else {
            cout<<"refracted";
        }
        cout << " from:to " << from << ":" << to << endl;
    }
};

struct StratifiedEdge: Edge{
    StratifiedEdge(){}
    StratifiedEdge(Pos from, Pos to, StratumId stratum_id): Edge(from, to, stratum_id) {}
};

struct RefractedEdge: Edge{
    RefractedEdge(){}
    RefractedEdge(Pos from, Pos to): Edge(from, to) {}
};


struct StratifiedData {
    StratumId stratum_id=0;
    Pos pos=0;
    
    StratifiedData(){}
    StratifiedData(Pos pos, StratumId stratum_id): pos(pos), stratum_id(stratum_id)
    {}
};

struct Vertex {
    uint32_t id;
    uint32_t size;
    vector<StratifiedEdge> s_edges;

    bool is_sequence;
    bool is_stratum;

    Vertex(){}
    Vertex(uint32_t id, uint32_t size, vector<StratifiedEdge> &regions, bool is_stratum) :id(id),size(size),s_edges(regions),is_stratum(is_stratum), is_sequence(!is_stratum){}
    Vertex(uint32_t id, uint32_t size, bool is_stratum, StratifiedData* data, StratifyingRegionCount rgn_cnt, vector<StratumSize> &stratum_sizes): id(id), size(size), is_stratum(is_stratum), is_sequence(!is_stratum){
        s_edges.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData d=data[rgn_cnt];
            auto rgn=StratifiedEdge(d.pos, d.pos+stratum_sizes[d.stratum_id], d.stratum_id);
            s_edges[rgn_cnt]=rgn;
            assert(0<=rgn.from && rgn.from < size && 0<rgn.to && rgn.to <= size);
        } 
    }
    Vertex(uint32_t id, uint32_t size, bool is_stratum, StratifiedData* data, StratifyingRegionCount rgn_cnt, int k, vector<StratumSize> &stratum_sizes): id(id), size(size), is_stratum(is_stratum), is_sequence(!is_stratum){
        int valid_rgn_cnt=0;
        for(int i=0; i<rgn_cnt; i++){
            if(stratum_sizes[data[i].stratum_id]>=k){
                valid_rgn_cnt++;
            }
        }
        s_edges.resize(valid_rgn_cnt);
        while(valid_rgn_cnt>0){
            rgn_cnt--;
            if(stratum_sizes[data[rgn_cnt].stratum_id]<k){
                continue;
            }
            valid_rgn_cnt--;
            StratifiedData d=data[rgn_cnt];
            auto rgn=StratifiedEdge(d.pos, d.pos+stratum_sizes[d.stratum_id], d.stratum_id);
            s_edges[valid_rgn_cnt]=rgn;
            // assert(0<=rgn.from && rgn.from < size && 0<rgn.to && rgn.to <= size && rgn.from < rgn.to);
        } 
    }
    void get_valid_indices(int k, vector<StratifyingRegionIdx>& indices){
        indices.clear();
        for(StratifyingRegionIdx i=0; i<s_edges.size(); i++){
            if(s_edges[i].size()>=k){
                indices.push_back(i);
            }
        }
    }
    StratumSize overlap_length_on_left(StratifyingRegionIdx idx){
        if(idx==0) {
            return 0;
        } else {
            return s_edges[idx-1].to > s_edges[idx].from? s_edges[idx-1].to-s_edges[idx].from : 0; 
        }
    }

    void print(){
        if(is_sequence) cout << "sequence("<< id <<"): " << endl;
        if(is_stratum)  cout << "stratum("<< id <<"): " << endl; 
        for(auto &r: s_edges){
            cout << "stratified: (" << r.from << ", " << r.to << ") - " << r.stratum_id << endl;
        }
    }
};

struct StratumVertex: Vertex{
    StratumVertex(){}
    StratumVertex(uint32_t id, uint32_t size, vector<StratifiedEdge> &regions): Vertex(id, size, regions, true){}
    StratumVertex(uint32_t id, uint32_t size, StratifiedData* data, StratifyingRegionCount rgn_cnt, vector<StratumSize> &stratum_sizes): Vertex(id, size, true, data, rgn_cnt, stratum_sizes){}
    StratumVertex(uint32_t id, uint32_t size, StratifiedData* data, StratifyingRegionCount rgn_cnt, int k, vector<StratumSize> &stratum_sizes): Vertex(id, size, true, data, rgn_cnt, k, stratum_sizes){}
};

struct SequenceVertex: Vertex{
    SequenceVertex(){}
    SequenceVertex(uint32_t id, uint32_t size, vector<StratifiedEdge> &regions): Vertex(id, size, regions, false){}
    SequenceVertex(uint32_t id, uint32_t size, StratifiedData* data, StratifyingRegionCount rgn_cnt, vector<StratumSize> &stratum_sizes): Vertex(id, size, false, data, rgn_cnt, stratum_sizes){}
    SequenceVertex(uint32_t id, uint32_t size, StratifiedData* data, StratifyingRegionCount rgn_cnt, int k, vector<StratumSize> &stratum_sizes): Vertex(id, size, false, data, rgn_cnt, k, stratum_sizes){}
};

// int CHECKSUM=stoi(string("I am a prokrustean file"));

struct Prokrustean {
    string file_name;
    int kmin;
    int version=1; // for file loading compatibility 
    bool contains_stratum_extension_count;
    bool contains_stratum_frequency;

    uint64_t sequence_count;
    uint64_t total_sequence_region_count=0;
    vector<SequenceSize> sequences__size;
    vector<StratifiedData*> sequences__region;
    vector<StratifyingRegionCount> sequences__region_cnt;
    
    uint64_t stratum_count;
    uint64_t total_strata_region_count=0;
    vector<StratumSize> stratums__size;
    vector<StratifiedData*> stratums__region;
    vector<StratifyingRegionCount> stratums__region_cnt;

    // extension
    vector<CharCount> stratums__character_extensions_cnt;
    vector<FrequencyCount> stratums__frequency_cnt;

    SequenceVertex get_sequence(SeqId id){
        auto sequence=SequenceVertex(id, sequences__size[id], sequences__region[id], sequences__region_cnt[id], stratums__size);
        assert(sequence.is_sequence);
        return sequence;
    }

    StratumVertex get_stratum(StratumId id){
        auto stratum=StratumVertex(id, stratums__size[id], stratums__region[id], stratums__region_cnt[id], stratums__size);
        assert(stratum.is_stratum);
        return stratum;
    }

    SequenceVertex get_sequence(SeqId id, int k){
        auto sequence=SequenceVertex(id, sequences__size[id], sequences__region[id], sequences__region_cnt[id], k, stratums__size);
        assert(sequence.is_sequence);
        return sequence;
    }

    StratumVertex get_stratum(StratumId id, int k){
        auto stratum=StratumVertex(id, stratums__size[id], stratums__region[id], stratums__region_cnt[id], k, stratums__size);
        assert(stratum.is_stratum);
        return stratum;
    }

    void get_sequence(SeqId id, Vertex &vertex, int k){
        // memory efficient
        vertex.id=id;
        vertex.size= sequences__size[id];
        vertex.is_sequence=true;
        vertex.is_stratum=false;
        auto rgn_cnt=sequences__region_cnt[id];
        vertex.s_edges.resize(rgn_cnt);
        auto valid_rgn_idx=0;
        for(int i=0; i<rgn_cnt; i++){
            StratifiedData &d=sequences__region[id][i];
            if(k<=stratums__size[d.stratum_id]){
                StratifiedEdge &edge=vertex.s_edges[valid_rgn_idx];
                edge.from=d.pos;
                edge.to=d.pos+stratums__size[d.stratum_id];
                edge.stratum_id=d.stratum_id;
                edge.is_stratified=true;
                edge.is_reflected=false;
                valid_rgn_idx++;
            }
        }
        if(valid_rgn_idx<rgn_cnt){
            vertex.s_edges.resize(valid_rgn_idx);
        }
    }

    void get_stratum(StratumId id, Vertex &vertex, int k){
        // memory efficient
        vertex.id=id;
        vertex.size= stratums__size[id];
        vertex.is_sequence=false;
        vertex.is_stratum=true;
        auto rgn_cnt=stratums__region_cnt[id];
        vertex.s_edges.resize(rgn_cnt);
        auto valid_rgn_idx=0;
        for(int i=0; i<rgn_cnt; i++){
            StratifiedData &d=stratums__region[id][i];
            if(k<=stratums__size[d.stratum_id]){
                StratifiedEdge &edge=vertex.s_edges[valid_rgn_idx];
                edge.from=d.pos;
                edge.to=d.pos+stratums__size[d.stratum_id];
                edge.stratum_id=d.stratum_id;
                edge.is_stratified=true;
                edge.is_reflected=false;
                valid_rgn_idx++;
            }
        }
        if(valid_rgn_idx<rgn_cnt){
            vertex.s_edges.resize(valid_rgn_idx);
        }
    }

    void get_sequence(SeqId id, Vertex &vertex){
        // memory efficient
        vertex.id=id;
        vertex.size= sequences__size[id];
        vertex.is_sequence=true;
        vertex.is_stratum=false;
        auto rgn_cnt=sequences__region_cnt[id];
        vertex.s_edges.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData &d=sequences__region[id][rgn_cnt];
            vertex.s_edges[rgn_cnt].from=d.pos;
            vertex.s_edges[rgn_cnt].to=d.pos+stratums__size[d.stratum_id];
            vertex.s_edges[rgn_cnt].stratum_id=d.stratum_id;
            vertex.s_edges[rgn_cnt].is_stratified=true;
            vertex.s_edges[rgn_cnt].is_reflected=false;
        }
    }

    void get_stratum(StratumId id, Vertex &vertex){
        // memory efficient
        vertex.id=id;
        vertex.size= stratums__size[id];
        vertex.is_sequence=false;
        vertex.is_stratum=true;
        auto rgn_cnt=stratums__region_cnt[id];
        vertex.s_edges.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData &d=stratums__region[id][rgn_cnt];
            vertex.s_edges[rgn_cnt].from=d.pos;
            vertex.s_edges[rgn_cnt].to=d.pos+stratums__size[d.stratum_id];
            vertex.s_edges[rgn_cnt].stratum_id=d.stratum_id;
            vertex.s_edges[rgn_cnt].is_stratified=true;
            vertex.s_edges[rgn_cnt].is_reflected=false;
        }
    }

    void get_spectrum(Vertex &v, int k,vector<Edge> &output){
        output.clear();
        assert(v.size>=k);
        // cout << "---- spectrum of -----" << endl;
        // v.print();
        vector<StratifiedEdge> regions_least_k;
        for(auto &r:v.s_edges){
            if(r.size()>=k){
                regions_least_k.push_back(r);
            }
        }
        int cnt=regions_least_k.size();
        // single reflected 
        if(cnt==0){
            // cout << "single reflected" << endl;
            auto rgn=RefractedEdge(0, v.size);
            output.push_back(rgn);
            assert(rgn.is_reflected&&!rgn.is_stratified);
            return;
        }

        // regions are already sorted
        for(int i=0; i<cnt; i++){
            if(i==0 && regions_least_k[i].from>0){
                output.push_back(RefractedEdge(0, regions_least_k[i].from+(k-1)));  
            }

            output.push_back(regions_least_k[i]);

            if(i<cnt-1){
                // careful - unsigned integers can flip the result
                if(regions_least_k[i+1].from <= regions_least_k[i].to && regions_least_k[i].to - regions_least_k[i+1].from >= k-1){
                    // large intersection
                } else {
                    output.push_back(RefractedEdge(regions_least_k[i].to-(k-1), regions_least_k[i+1].from+(k-1)));          
                }
            }

            if(i==cnt-1 && regions_least_k[i].to<v.size){
                output.push_back(RefractedEdge(regions_least_k[i].to-(k-1), v.size));  
            }
        }

        // for(auto &rgn: output){
        //     rgn.print();
        // }
    } 
    
    StratumSize get_stratum_size(StratumId id){
        return this->stratums__size[id];
    }
    CharCount get_left_cnt(StratumId id) {
        // Extract the first 'numBits' bits (leftmost)
        return (this->stratums__character_extensions_cnt[id] >> HALF_BITS_IN_CHAR) & MASK_HALF_BITS;
    }
    CharCount get_right_cnt(StratumId id) {
        // the next 'numBits' bits (rightmost)
        return this->stratums__character_extensions_cnt[id] & MASK_HALF_BITS;
    }

    void set_stratum_left_right_extensions(StratumId id, CharCount left_cnt, CharCount right_cnt) {
        // Ensure that both integers are within the valid range
        if (left_cnt >= (1 << HALF_BITS_IN_CHAR) || right_cnt >= (1 << HALF_BITS_IN_CHAR)) {
            // Handle invalid input
            throw std::invalid_argument("left cnt, right cnt too large. The data space is designed to include both left and right extension counts in one byte");
        }
        // Shift the first integer to the left and combine it with the second integer using bitwise OR
        this->stratums__character_extensions_cnt[id]=(left_cnt << HALF_BITS_IN_CHAR) | (right_cnt & MASK_HALF_BITS);
    }

    void set_stratum_frequency(StratumId id, FrequencyCount frequency) {
        assert(frequency>0);
        this->stratums__frequency_cnt[id]=frequency;
    }

    void set_seq_count(uint64_t seq_cnt){
        this->sequences__size.resize(seq_cnt);
        this->sequences__region.resize(seq_cnt);
        this->sequences__region_cnt.resize(seq_cnt, 0);
        this->sequence_count=seq_cnt;
    }

    void set_stratum_count(uint64_t stratum_cnt){
        this->stratums__region.resize(stratum_cnt);
        this->stratums__region_cnt.resize(stratum_cnt, 0);
        this->stratums__size.resize(stratum_cnt); 
        this->stratum_count=stratum_cnt;
        if(contains_stratum_extension_count){
            this->stratums__character_extensions_cnt.resize(stratum_cnt);
        }
        if(contains_stratum_frequency){
            this->stratums__frequency_cnt.resize(stratum_cnt);
        }
    }

    void set_stratum_size(StratumId id, StratumSize stratum_size){
        this->stratums__size[id]=stratum_size;
    }

    void set_seq_size(SeqId id, uint64_t seq_size){
        this->sequences__size[id]=seq_size;
    }

    void set_seq_regions(SeqId id, StratifiedData* data, uint64_t rgn_cnt){
        assert(rgn_cnt<=numeric_limits<StratifyingRegionCount>::max());
        this->sequences__region[id]=data;
        this->sequences__region_cnt[id]=rgn_cnt;
    }

    void set_stratum_regions(StratumId id, StratifiedData* data, uint64_t rgn_cnt){
        assert(rgn_cnt<=numeric_limits<StratifyingRegionCount>::max());
        this->stratums__region[id]=data;
        this->stratums__region_cnt[id]=rgn_cnt;
    }

    void print_abstract(){
        cout << "------------------------------------------------------" << endl;
        cout << "---------------    prokrustean    --------------------" << endl;
        cout << "------------------------------------------------------" << endl;
        cout << "-- " << "sequences: " << sequence_count << endl;
        cout << "-- " << "stratified regions in sequences: " << total_sequence_region_count << endl;
        cout << "-- " << "strata: " << stratum_count << endl;
        cout << "-- " << "stratified regions in strata: " << total_strata_region_count << endl;
        cout << "-- " << "kmin used: " << kmin << endl;
        cout << "------------------------------------------------------" << endl;
    }
};


#endif