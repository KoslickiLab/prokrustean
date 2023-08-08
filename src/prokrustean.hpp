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

class MinCover {
    bool is_rep;
    uint64_t id;
    uint64_t size;
    vector<tuple<uint64_t, RepId>> mc_rep;
public:
    MinCover(){}

    void addCovOcc(RepId id, Pos pos){
    }
};

class Prokrustean {
    vector<MinCover> seq_mc;
    vector<MinCover> rep_mc;
public:
    Prokrustean(){}
};

#endif