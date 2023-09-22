#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

/* Important: the types were fixed due to the purpose of the work. Need to be adjusted based on needs */

typedef uint32_t RepId; // 4,294,967,295
typedef uint32_t SeqId; // 4,294,967,295
typedef uint32_t Id; // 4,294,967,295 Both RepId SeqId
typedef uint64_t Pos; // 4,294,967,295 == Sequence max length
typedef uint16_t Size; //65,535 maximal repeat size

struct Stratification {
    uint64_t id;
    bool is_stratum;
    uint64_t size;
    vector<tuple<Pos, RepId>> regions;
};

struct Prokrustean {
    vector<Stratification> seqs;
    vector<Stratification> stratums;
    
    //optional
    optional<vector<string>> sequences;

    void set_sizes(uint64_t seq_no, uint64_t rep_no, bool recover_sequences = false){
        seqs.resize(seq_no);
        stratums.resize(rep_no);
        if(recover_sequences)
        sequences = vector<string>(seq_no);
    }
};

void print_prokrustean_statistics(Prokrustean pk){
    cout << "--- pk stat ---" << endl;
    cout << "seq mc: " << pk.seqs.size() << endl;
    uint64_t occurrence=0;
    for(auto mc: pk.seqs){
        occurrence+=mc.regions.size();
    }
    cout << "seq mc occurrences: " << occurrence << endl;

    cout << "rep mc: " << pk.stratums.size() << endl;
    occurrence=0;
    for(auto mc: pk.stratums){
        occurrence+=mc.regions.size();
    }
    cout << "rep mc occurrences: " << occurrence << endl;
    
}

void _print_rep(RepId rid, string str, int depth, Prokrustean &pk, vector<bool> &printed_rep){
    if(printed_rep[rid]){
        cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
        if(pk.stratums[rid].regions.size()>0){
            cout << string(2*(depth+1), '-')<< "...dup skip..." << endl;
        }
        return;
    }
    printed_rep[rid]=true;
    cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
    for(auto r: pk.stratums[rid].regions){
        Pos _pos = get<0>(r);
        RepId _rid = get<1>(r);
        _print_rep(_rid, str.substr(_pos, pk.stratums[_rid].size), depth+1, pk, printed_rep);
    }
}

void print_prokrustean(Prokrustean pk){
    assert(pk.sequences.has_value());
    cout << "-- print prokrustean --" << endl;
    cout << "rep no.: " << pk.stratums.size() << endl;
    vector<bool> printed_rep(pk.stratums.size());
    for(uint64_t i=0; i<pk.seqs.size(); i++){
        auto mc_reps = pk.seqs[i].regions;
        cout << "seq" << i << ", size " << pk.seqs[i].size << ", reps " << mc_reps.size() << endl;
        cout << pk.sequences.value()[i] << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            auto str = pk.sequences.value()[i].substr(pos, pk.stratums[rid].size);
            _print_rep(rid, str, 1, pk, printed_rep);
        }
        cout << endl;
    }
}

void print_bare_prokrustean(Prokrustean pk){
    for(uint64_t i=0; i<pk.seqs.size(); i++){
        auto mc_reps = pk.seqs[i].regions;
        cout << "seq" << i << ", size " << pk.seqs[i].size << ", reps " << mc_reps.size() << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            cout << "(" << "pos: " << pos << ", " << "size: " << pk.stratums[rid].size<< ", R" << rid << ")" << " ";
        }
        cout << endl;
    }
    for(uint64_t i=0; i<pk.stratums.size(); i++){
        auto mc_reps = pk.stratums[i].regions;
        cout << "rep" << i << ", size " << pk.stratums[i].size << ", reps " << mc_reps.size() << endl;
        for(auto r: mc_reps){
            Pos pos = get<0>(r);
            RepId rid = get<1>(r);
            cout << "(" << pos << ", " << pk.stratums[rid].size<< ", R" << rid << ")" << " ";
        }
        cout << endl;
    }
}

void print_single_rep_to_single_rep_relationships(Prokrustean pk){
    int cnt = 0;
    for(auto mc: pk.stratums){
        if(mc.regions.size()==1){
            auto unique_repid = get<1>(mc.regions[0]);
            if(pk.stratums[unique_repid].regions.size()>0){
                cnt++;
            }
        }
    }
    cout << "-- print_single_rep_to_single_rep_relationships --" << endl;
    cout << "total: " << cnt << endl;
}
#endif