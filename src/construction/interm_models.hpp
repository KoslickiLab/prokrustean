#ifndef CONSTRUCTION_REPR_SA_HPP_
#define CONSTRUCTION_REPR_SA_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

class ReprSASimple: public ReprSABase {
    std::vector<bool> sa_mark;
    std::unordered_map<uint64_t, vector<uint64_t>> dictionary;

public:
    ReprSASimple(uint64_t n){
        // bit vector of size the total sequence length
        sa_mark = vector<bool>(n);
    }
    bool is_repr(uint64_t sa_idx){
        return sa_mark[sa_idx];
    };
    // returns repeat id list
    vector<uint64_t> pop_repr(uint64_t sa_idx){
        if(~sa_mark[sa_idx]){
            throw std::invalid_argument("the sa is not marked. Is not a reprsentative SA or has already been popped." );    
        }
        vector<uint64_t> repeats = dictionary[sa_idx];
        sa_mark[sa_idx] = false;
        return repeats;
    };

    void insert_repr(uint64_t sa_idx, uint64_t repeat_idx){
        sa_mark[sa_idx] = true;
        // mutex required.
        if (dictionary.find(sa_idx) == dictionary.end()){
            dictionary[sa_idx] = {repeat_idx};
        }
    };
};

#endif;