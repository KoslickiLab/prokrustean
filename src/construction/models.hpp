#ifndef CONSTRUCTION_MODEL_HPP_
#define CONSTRUCTION_MODEL_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../prokrustean.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

struct MaximalRepeatAnnotation {
    // 
    uint64_t size;

    vector<SuffixArrayIdx> repr_indexes;

    SuffixArrayIdx first_repr_idx;
};

struct PositionAnnotation {
    Pos pos; 
    
    SuffixArrayIdx sa_idx; 

    // order desc for rep sizes
    vector<MaximalRepeatAnnotation> reps;

    // corresponding ids
    vector<RepId> rep_ids;
};

struct SequenceAnnotation {
    //
    SeqId id;

    // 
    uint64_t size;

    //
    vector<PositionAnnotation> positions;

    //
    optional<string> sequence; 
};

struct ReprSuffixRank {
    //
    sdsl::bit_vector bv;
    //
    sdsl::rank_support_v<> rb;

    ReprSuffixRank(){}

    void initialize(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        bv.resize(seq_length);
        for(auto rep: repeats){
            for(auto sa_idx: rep.repr_indexes){
                bv[sa_idx]=true;
            }
        }
        
        rank_support_v<> _rb(&bv);
        this->rb = _rb;
    }
    
    uint64_t get_repr_size(){
        return rb.rank(bv.size());
    }

    // 
    bool exists(uint64_t sa_idx){
        return bv[sa_idx];
    }

    uint64_t rank(uint64_t sa_idx){
        return rb.rank(sa_idx);
    }
};

struct ReprSuffixAnnotation {
private:
    ReprSuffixRank sa_rank;
    vector<vector<RepId>> repr_suffixes;

public:
    void initialize_rank(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        sa_rank.initialize(seq_length, repeats);
    }

    void initialize_repr_sa(vector<MaximalRepeatAnnotation> &repeats){
        repr_suffixes.resize(sa_rank.get_repr_size());
        for(uint64_t rep_id=0; rep_id< repeats.size(); rep_id++){
            for(auto sa_idx: repeats[rep_id].repr_indexes){
                uint64_t r = sa_rank.rank(sa_idx);
                repr_suffixes[r].push_back(rep_id);
            }
        }
    }

    optional<vector<RepId>> get_repeats(uint64_t sa_idx){
        if(!sa_rank.exists(sa_idx)){
            return nullopt;
        } else {
            auto rank = sa_rank.rank(sa_idx);
            return repr_suffixes[rank];
        }
    }
};



class ReprSuffixAnnotationParallel {
    /* thread safe imple*/
    vector<vector<RepId>> repr_suffixes;

    ReprSuffixAnnotationParallel(uint64_t repr_size){
        repr_suffixes.resize(repr_size);
    }

    // 
    void set_repeats(vector<MaximalRepeatAnnotation> &repeats, uint64_t from, uint64_t to){

    }

    // 
    vector<RepId> get_repeats(uint64_t sa_idx){

    }
};

#endif