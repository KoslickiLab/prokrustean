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
    
    //optional
    optional<vector<string>> sequences;

    void set_sizes(uint64_t seq_no, uint64_t rep_no, bool recover_sequences = false){
        seq_mcs.resize(seq_no);
        rep_mcs.resize(rep_no);
        if(recover_sequences)
        sequences = vector<string>(seq_no);
    }
};

void _print_rep(RepId rid, string str, int depth, Prokrustean &pk, vector<bool> &printed_rep){
    if(printed_rep[rid]){
        cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
        if(pk.rep_mcs[rid].mc_reps.size()>0){
            cout << string(2*(depth+1), '-')<< "...dup skip..." << endl;
        }
        return;
    }
    printed_rep[rid]=true;
    cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
    for(auto r: pk.rep_mcs[rid].mc_reps){
        Pos _pos = get<0>(r);
        RepId _rid = get<1>(r);
        _print_rep(_rid, str.substr(_pos, pk.rep_mcs[_rid].size), depth+1, pk, printed_rep);
    }
}

void print_prokrustean(Prokrustean pk){
    assert(pk.sequences.has_value());
    vector<bool> printed_rep(pk.rep_mcs.size());
    for(uint64_t i=0; i<pk.seq_mcs.size(); i++){
        auto mc_reps = pk.seq_mcs[i].mc_reps;
        cout << "seq" << i << ", size " << pk.seq_mcs[i].size << ", reps " << mc_reps.size() << endl;
        cout << pk.sequences.value()[i] << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            auto str = pk.sequences.value()[i].substr(pos, pk.rep_mcs[rid].size);
            _print_rep(rid, str, 1, pk, printed_rep);
        }
        cout << endl;
    }
}

void print_bare_prokrustean(Prokrustean pk){
    for(uint64_t i=0; i<pk.seq_mcs.size(); i++){
        auto mc_reps = pk.seq_mcs[i].mc_reps;
        cout << "seq" << i << ", size " << pk.seq_mcs[i].size << ", reps " << mc_reps.size() << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            cout << "(" << "pos: " << pos << ", " << "size: " << pk.rep_mcs[rid].size<< ", R" << rid << ")" << " ";
        }
        cout << endl;
    }
    for(uint64_t i=0; i<pk.rep_mcs.size(); i++){
        auto mc_reps = pk.rep_mcs[i].mc_reps;
        cout << "rep" << i << ", size " << pk.rep_mcs[i].size << ", reps " << mc_reps.size() << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            cout << "(" << pos << ", " << pk.rep_mcs[rid].size<< ", R" << rid << ")" << " ";
        }
        cout << endl;
    }
}
#endif