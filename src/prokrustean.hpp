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

    MinCover getMC(SeqId id){
    }
    MinCover getMC(RepId id){
    }
    RepId registerRep(uint64_t size){
        rep_mc.push_back(MinCover());
        return rep_mc.size();
    }
};