#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <stack>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "data_types.hpp"

using namespace std;

struct Region {
    Pos from;
    Pos to;

    bool is_stratified;
    bool is_reflected;
    // stratified only 
    StratumId stratum_id;

    Region(){}
    Region(Pos from, Pos to): Region(from, to, false, 0){}
    Region(Pos from, Pos to, StratumId stratum_id): Region(from, to, true, stratum_id){}
    Region(Pos from, Pos to, bool is_stratified, StratumId stratum_id): from(from), to(to), is_stratified(is_stratified), is_reflected(!is_stratified), stratum_id(stratum_id) {}

    uint64_t size(){
        assert(from<to);
        return to-from;
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

struct StratifiedRegion: Region{
    StratifiedRegion(){}
    StratifiedRegion(Pos from, Pos to, StratumId stratum_id): Region(from, to, stratum_id) {}
};

struct ReflectedRegion: Region{
    ReflectedRegion(){}
    ReflectedRegion(Pos from, Pos to): Region(from, to) {}
};


struct StratifiedData {
    StratumId stratum_id;
    Pos pos;
    
    StratifiedData(){}
    StratifiedData(Pos pos, StratumId stratum_id): pos(pos), stratum_id(stratum_id)
    {}
};

struct Vertex {
    uint32_t id;
    uint32_t size;
    vector<StratifiedRegion> s_edges;

    bool is_sequence;
    bool is_stratum;

    Vertex(){}
    Vertex(uint32_t id, uint32_t size, vector<StratifiedRegion> &regions, bool is_stratum) :id(id),size(size),s_edges(regions),is_stratum(is_stratum), is_sequence(!is_stratum){}
    Vertex(uint32_t id, uint32_t size, bool is_stratum, StratifiedData* data, uint8_t rgn_cnt, vector<StratumSize> &stratum_sizes): id(id), size(size), is_stratum(is_stratum), is_sequence(!is_stratum){
        s_edges.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData d=data[rgn_cnt];
            auto rgn=StratifiedRegion(d.pos, d.pos+stratum_sizes[d.stratum_id], d.stratum_id);
            s_edges[rgn_cnt]=rgn;
            assert(0<=rgn.from && rgn.from < size && 0<rgn.to && rgn.to <= size);
        } 
    }
    void get_valid_indices(int k, vector<uint8_t>& indices){
        indices.clear();
        for(uint8_t i=0; i<s_edges.size(); i++){
            if(s_edges[i].size()>=k){
                indices.push_back(i);
            }
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

struct Stratum: Vertex{
    Stratum(){}
    Stratum(uint32_t id, uint32_t size, vector<StratifiedRegion> &regions): Vertex(id, size, regions, true){}
    Stratum(uint32_t id, uint32_t size, StratifiedData* data, uint8_t rgn_cnt, vector<StratumSize> &stratum_sizes): Vertex(id, size, true, data, rgn_cnt, stratum_sizes){}
};

struct Sequence: Vertex{
    Sequence(){}
    Sequence(uint32_t id, uint32_t size, vector<StratifiedRegion> &regions): Vertex(id, size, regions, false){}
    Sequence(uint32_t id, uint32_t size, StratifiedData* data, uint8_t rgn_cnt, vector<StratumSize> &stratum_sizes): Vertex(id, size, false, data, rgn_cnt, stratum_sizes){}
};


class SpinLock {
    std::atomic_flag flag = ATOMIC_FLAG_INIT;

public:
    SpinLock(){}
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire)) {}
    }

    void unlock() {
        flag.clear(std::memory_order_release);
    }
};


struct Prokrustean {
    uint64_t lmin;
    uint64_t sequence_count;
    vector<SequenceSize> sequences__size;
    vector<StratifiedData*> sequences__region;
    vector<uint8_t> sequences__region_cnt;
    
    uint64_t stratum_count;
    vector<StratumSize> stratums__size;
    vector<StratifiedData*> stratums__region;
    vector<uint8_t> stratums__region_cnt;


    Sequence get_sequence(SeqId id){
        auto sequence=Sequence(id, sequences__size[id], sequences__region[id], sequences__region_cnt[id], stratums__size);
        assert(sequence.is_sequence);
        return sequence;
    }

    Stratum get_stratum(StratumId id){
        auto stratum=Stratum(id, stratums__size[id], stratums__region[id], stratums__region_cnt[id], stratums__size);
        assert(stratum.is_stratum);
        return stratum;
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

    void get_stratum(SeqId id, Vertex &vertex){
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

    void get_spectrum(Vertex &v, int k,vector<Region> &output){
        output.clear();
        assert(v.size>=k);
        // cout << "---- spectrum of -----" << endl;
        // v.print();
        vector<StratifiedRegion> regions_least_k;
        for(auto &r:v.s_edges){
            if(r.size()>=k){
                regions_least_k.push_back(r);
            }
        }
        int cnt=regions_least_k.size();
        // single reflected 
        if(cnt==0){
            // cout << "single reflected" << endl;
            auto rgn=ReflectedRegion(0, v.size);
            output.push_back(rgn);
            assert(rgn.is_reflected&&!rgn.is_stratified);
            return;
        }

        // regions are already sorted
        for(int i=0; i<cnt; i++){
            if(i==0 && regions_least_k[i].from>0){
                output.push_back(ReflectedRegion(0, regions_least_k[i].from+(k-1)));  
            }

            output.push_back(regions_least_k[i]);

            if(i<cnt-1){
                // careful - unsigned integers can flip the result
                if(regions_least_k[i+1].from <= regions_least_k[i].to && regions_least_k[i].to - regions_least_k[i+1].from >= k-1){
                    // large intersection
                } else {
                    output.push_back(ReflectedRegion(regions_least_k[i].to-(k-1), regions_least_k[i+1].from+(k-1)));          
                }
            }

            if(i==cnt-1 && regions_least_k[i].to<v.size){
                output.push_back(ReflectedRegion(regions_least_k[i].to-(k-1), v.size));  
            }
        }

        // for(auto &rgn: output){
        //     rgn.print();
        // }
    } 
    
    StratumSize get_stratum_size(StratumId id){
        return this->stratums__size[id];
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
        // if in the middle of construction, stratum sizes are already set though
        this->stratums__size.resize(stratum_cnt);
        this->stratum_count=stratum_cnt;
    }

    void set_stratum_size(StratumSize stratum_size){
        this->stratums__size.push_back(stratum_size);
    }

    void set_seq_regions(SeqId id, uint64_t seq_size, StratifiedData* data, uint8_t rgn_cnt){
        this->sequences__size[id]=seq_size;
        this->sequences__region[id]=data;
        this->sequences__region_cnt[id]=rgn_cnt;
    }

    void set_stratum_regions(StratumId id, StratifiedData* data, uint8_t rgn_cnt){
        this->stratums__region[id]=data;
        this->stratums__region_cnt[id]=rgn_cnt;
    }
};


#endif