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


class ReprSuffixAnnotation {
    bm::bvector<> bv;
    std::unique_ptr<bm::bvector<>::rs_index_type> rs_idx;
    uint64_t repr_suff_cnt;
    vector<vector<RepId>> repr_suffixes;
public:
    ReprSuffixAnnotation(uint64_t seq_length){
        bv.resize(seq_length);
        auto rs_idx(new bm::bvector<>::rs_index_type());
    }

    // 
    void set_repr_sa_rank(vector<MaximalRepeatAnnotation> &repeats){
        for(auto rep: repeats){
            for(auto sa_idx: rep.repr_sa_indexes){
                bv.set_bit(sa_idx);
            }
        }
        bv.build_rs_index(rs_idx.get());

        //set repr suffix cardinality
        repr_suff_cnt = bv.rank(bv.size(), *rs_idx);
        repr_suffixes.resize(repr_suff_cnt);
    }

    void set_repr_sa_annotation(vector<MaximalRepeatAnnotation> &repeats){
        uint64_t rep_id = 0;
        for(auto rep: repeats){
            for(auto sa_idx: rep.repr_sa_indexes){
                repr_suffixes[sa_idx].push_back(rep_id);
            }
            rep_id++;
        }
    }
    // 
    optional<vector<RepId>> get_repeats(uint64_t sa_idx){
        if(~bv.get_bit(sa_idx)){
            return nullopt;
        }

        uint64_t rank = bv.rank(sa_idx, *rs_idx);
        return repr_suffixes[rank];
    }
};

class ReprSuffixAnnotationParallel {
    /* thread safe imple*/
public:
    ReprSuffixAnnotationParallel(uint64_t sa_count){

    }

    // 
    void set_repeats(vector<MaximalRepeatAnnotation> &repeats, uint64_t from, uint64_t to){

    }

    // 
    vector<RepId> get_repeats(uint64_t sa_idx){

    }
};

#endif