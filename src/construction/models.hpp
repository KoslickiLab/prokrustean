#ifndef CONSTRUCTION_MODEL_HPP_
#define CONSTRUCTION_MODEL_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "build_prokrustean.hpp"
#include "../prokrustean.hpp"

using namespace std;

class SimpleReprSuffPositions: public AbsReprSuffPositions {
    std::vector<bool> sa_mark;

public:
    SimpleReprSuffPositions(uint64_t n){
        // bit vector of size the total sequence length
        sa_mark = vector<bool>(n);
    }

    // 
    uint64_t count(){

    }

    // 
    void set(uint64_t sa_idx){

    };

    // 
    bool exists(uint64_t sa_idx){

    };

    // 
    uint64_t rank(uint64_t sa_idx){

    };
};

class SimpleReprSuff: public AbsReprSuff {
    void initialize(uint64_t repr_sa_cnt){

    }

    //
    void set(uint64_t sa_index, RepId rep_id){

    }

    // 
    vector<RepId> get(uint64_t sa_index, bool remove=true){

    }
};

#endif