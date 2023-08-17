#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

typedef uint64_t RepId;
typedef uint64_t SeqId;
typedef uint64_t Pos;  

struct MinCover {
    uint64_t id;
    bool is_rep;
    uint64_t size;
    vector<tuple<Pos, RepId>> mc_reps;
};

struct Prokrustean {
    vector<MinCover> seq_mcs;
    vector<MinCover> rep_mcs;
    
    optional<vector<string>> sequences;

    void set_sizes(uint64_t seq_no, uint64_t rep_no){
        seq_mcs.resize(seq_no);
        rep_mcs.resize(rep_no);
    }
};

#endif