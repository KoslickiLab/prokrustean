#ifndef BASE_HPP_
#define BASE_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

/* These types influence the size of outputs A LOT */

typedef uint32_t StratumId; // 4,294,967,295
typedef uint32_t SeqId; // 4,294,967,295
typedef uint16_t Pos; // 4,294,967,295 == Sequence max length
typedef uint16_t StratumSize; //65,535
typedef uint32_t SequenceSize; //65,535


struct Region {
    Pos from;
    Pos to;

    bool is_stratified;
    bool is_reflected;
    // stratified only 
    StratumId stratum_id;

    Region(){}
    Region(Pos from, Pos to, bool is_stratified): Region(from, to, is_stratified, 0){}
    Region(Pos from, Pos to, bool is_stratified, StratumId stratum_id): from(from), to(to), is_stratified(is_stratified), is_reflected(!is_stratified), stratum_id(stratum_id) {}

    uint64_t size(){
        assert(from<to);
        return to-from;
    }
};

struct StratifiedRegion: Region{
    StratifiedRegion(){}
    StratifiedRegion(Pos from, Pos to, StratumId stratum_id): Region(from, to, true, stratum_id) {}
};

struct ReflectedRegion: Region{
    ReflectedRegion(){}
    ReflectedRegion(Pos from, Pos to): Region(from, to, false) {}
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
    vector<StratifiedRegion> s_regions;

    bool is_sequence;
    bool is_stratum;
    
    // stratum only
    optional<tuple<SeqId, Pos>> example_occ;

    void set_occ(SeqId seq_id, Pos pos){
        example_occ=make_tuple(seq_id, pos);
    }
    Vertex(){}
    Vertex(uint32_t id, uint32_t size, vector<StratifiedRegion> &regions, bool is_stratum) :Vertex(id, size, regions, is_stratum, nullopt){}
    Vertex(uint32_t id, uint32_t size, vector<StratifiedRegion> &regions, bool is_stratum, optional<tuple<SeqId, Pos>> example_occ) :id(id),size(size),s_regions(regions),is_stratum(is_stratum), is_sequence(!is_stratum), example_occ(example_occ) {}    
    Vertex(uint32_t id, uint32_t size, bool is_stratum, StratifiedData* data, uint8_t rgn_cnt, vector<StratumSize> &stratum_sizes): id(id), size(size), is_stratum(is_stratum), is_sequence(!is_stratum){
        s_regions.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData d=data[rgn_cnt];
            auto rgn=StratifiedRegion(d.pos, d.pos+stratum_sizes[d.stratum_id], d.stratum_id);
            s_regions[rgn_cnt]=rgn;
            assert(0<=rgn.from && rgn.from < size && 0<rgn.to && rgn.to <= size);
        } 
    }
    void get_valid_indices(int k, vector<uint8_t>& indices){
        indices.clear();
        for(uint8_t i=0; i<s_regions.size(); i++){
            if(s_regions[i].size()>=k){
                indices.push_back(i);
            }
        }
    }
    void print(){
        if(is_sequence) cout << "sequence("<< id <<"): " << endl;
        if(is_stratum)  cout << "stratum("<< id <<"): " << endl; 
        for(auto &r: s_regions){
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
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire)) {}
    }

    void unlock() {
        flag.clear(std::memory_order_release);
    }
};

#endif