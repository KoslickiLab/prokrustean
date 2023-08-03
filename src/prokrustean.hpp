#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

typedef uint64_t RepId;
typedef uint64_t SeqId;
typedef uint64_t Pos;  

class MinCover {
    uint64_t size;
    vector<pair<uint64_t, RepId>> mc_rep;
public:
    MinCover(){}

    void AddCovOcc(RepId id, Pos pos){
    }
};

class Prokrustean {
    vector<MinCover> seq_mc;
    vector<MinCover> rep_mc;
public:
    Prokrustean(){}

    MinCover GetMC(SeqId id){
    }
    MinCover GetMC(RepId id){
    }
};