#ifndef CONSTRUCTION_MODEL_HPP_
#define CONSTRUCTION_MODEL_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../prokrustean.hpp"
#include "../BitMagic/bm.h"

using namespace std;

struct MaximalRepeatAnnotation {
    // 
    uint64_t size;

    //
    vector<uint64_t> repr_sa_indexes;
};

struct SequenceAnnotation {
    // 
    uint64_t size;

    //
    vector<tuple<Pos, vector<RepId>>> annotations;
};

struct ReprSuffixRank {
    //
    bm::bvector<> bv;
    //
    bm::bvector<>::rs_index_type* rs_idx;

    ReprSuffixRank(){}

    ReprSuffixRank(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        this->bv.resize(seq_length);
        auto _rs_idx(new bm::bvector<>::rs_index_type());
        this->rs_idx = _rs_idx;

        for(auto rep: repeats)
        for(auto sa_idx: rep.repr_sa_indexes)
        bv.set_bit(sa_idx);
            
        bv.build_rs_index(rs_idx);
    }
    
    uint64_t get_repr_size(){
        return bv.rank(bv.size(), *rs_idx);
    }

    // 
    bool exists(uint64_t sa_idx){
        return bv.get_bit(sa_idx);
    }

    uint64_t rank(uint64_t sa_idx){
        return bv.rank(sa_idx, *rs_idx);
    }
};

struct ReprSuffixAnnotation {
private:
    ReprSuffixRank sa_rank;
    vector<vector<RepId>> repr_suffixes;

public:
    void initialize_rank(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        sa_rank=ReprSuffixRank(seq_length, repeats);
    }

    void initialize_repr_sa(vector<MaximalRepeatAnnotation> &repeats){
        repr_suffixes.resize(sa_rank.get_repr_size());

        for(uint64_t rep_id; rep_id< repeats.size(); rep_id++){
            for(auto sa_idx: repeats[rep_id].repr_sa_indexes){
                repr_suffixes[sa_idx].push_back(rep_id);
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