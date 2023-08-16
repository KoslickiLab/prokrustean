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
    uint64_t size;
    vector<tuple<Pos, RepId>> mc_rep;
};

struct Prokrustean {
    vector<MinCover> seq_mc;
    vector<MinCover> rep_mc;

    void set_sizes(uint64_t seq_no, uint64_t rep_no){
        seq_mc.resize(seq_no);
        rep_mc.resize(rep_no);
    }
};

#endif