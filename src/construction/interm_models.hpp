#ifndef CONSTRUCTION_REPR_SA_HPP_
#define CONSTRUCTION_REPR_SA_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "build_prokrustean.hpp"
#include "../prokrustean.hpp"

using namespace std;

class ReprSASimple: public AbstractReprSA {
    std::vector<bool> sa_mark;
    std::unordered_map<uint64_t, vector<RepId>> sa_dict;

public:
    ReprSASimple(uint64_t n){
        // bit vector of size the total sequence length
        sa_mark = vector<bool>(n);
    }

    bool exists(uint64_t sa_idx){
        return sa_mark[sa_idx];
    };
    // returns repeat id list
    vector<uint64_t> pop(uint64_t sa_idx){
        if(~sa_mark[sa_idx]){
            throw std::invalid_argument("the sa is not marked. Is not a reprsentative SA or has already been popped." );    
        }
        vector<uint64_t> repeats = sa_dict[sa_idx];
        sa_mark[sa_idx] = false;
        return repeats;
    };

    void push(RepId rep_idx, vector<uint64_t> sa_indexes){
        for(auto sa_idx: sa_indexes){
            sa_mark[sa_idx] = true;
            // mutex required.
            if (sa_dict.find(sa_idx) == sa_dict.end()){
                sa_dict[sa_idx] = {rep_idx};
            } else {
                sa_dict[sa_idx].push_back(rep_idx);
            }
        }
    };
};

#endif;