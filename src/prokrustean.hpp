#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "base.hpp"

using namespace std;

struct Region{
    StratumId stratum_id;
    Pos pos;
    Region(){}
    Region(Pos pos, StratumId stratum_id): pos(pos), stratum_id(stratum_id)
    {}
};

struct Stratum {
    StratumSize size;
    vector<Region> regions;
};

struct Sequence {
    SequenceSize size;
    vector<Region> regions;
};

struct Prokrustean {
    /* data structure is succinct as possible*/
    vector<SequenceSize> sequences__size;
    vector<Region*> sequences__region;
    vector<uint8_t> sequences__region_cnt;
    
    vector<StratumSize> stratums__size;
    vector<Region*> stratums__region;
    vector<uint8_t> stratums__region_cnt;
    
    //optional
    optional<vector<string>> sequences;

    uint64_t sequence_count(){
        return this->sequences__size.size();
    }

    uint64_t stratum_count(){
        return this->stratums__size.size();
    }

    Sequence get_sequence(SeqId id){
        auto sequence=Sequence();
        auto rgn_cnt=sequences__region_cnt[id];
        sequence.size=sequences__size[id];
        sequence.regions.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            sequence.regions[rgn_cnt]=*sequences__region[rgn_cnt];
        }
        return sequence;
    }

    Stratum get_stratum(StratumId id){
        auto stratum=Stratum();
        auto rgn_cnt=stratums__region_cnt[id];
        stratum.size=stratums__size[id];
        stratum.regions.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            stratum.regions[rgn_cnt]=*stratums__region[rgn_cnt];
        }
        return stratum;
    }
    
    uint64_t get_cardinality(){
        uint64_t cnt = 0;
        for(auto rcnt: sequences__region_cnt){
            cnt+=rcnt;
        }
        for(auto rcnt: stratums__region_cnt){
            cnt+=rcnt;
        }
        return cnt;
    }
};

// void print_prokrustean_statistics(Prokrustean pk){
//     cout << "--- pk stat ---" << endl;
//     cout << "seq mc: " << pk.seqs.size() << endl;
//     uint64_t occurrence=0;
//     for(auto mc: pk.seqs){
//         occurrence+=mc.regions.size();
//     }
//     cout << "seq mc occurrences: " << occurrence << endl;

//     cout << "rep mc: " << pk.stratums.size() << endl;
//     occurrence=0;
//     for(auto mc: pk.stratums){
//         occurrence+=mc.regions.size();
//     }
//     cout << "rep mc occurrences: " << occurrence << endl;
    
// }

// void _print_rep(StratumId rid, string str, int depth, Prokrustean &pk, vector<bool> &printed_rep){
//     if(printed_rep[rid]){
//         cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
//         if(pk.stratums[rid].regions.size()>0){
//             cout << string(2*(depth+1), '-')<< "...dup skip..." << endl;
//         }
//         return;
//     }
//     printed_rep[rid]=true;
//     cout << string(2*depth, '-')<< str << " (" << "R" << rid << ")"<< endl;
//     for(auto r: pk.stratums[rid].regions){
//         Pos _pos = r.pos;
//         StratumId _rid = r.stratum_id;
//         _print_rep(_rid, str.substr(_pos, pk.stratums[_rid].size), depth+1, pk, printed_rep);
//     }
// }

// void print_prokrustean(Prokrustean pk){
//     assert(pk.sequences.has_value());
//     cout << "-- print prokrustean --" << endl;
//     cout << "rep no.: " << pk.stratums.size() << endl;
//     vector<bool> printed_rep(pk.stratums.size());
//     for(uint64_t i=0; i<pk.seqs.size(); i++){
//         auto mc_reps = pk.seqs[i].regions;
//         cout << "seq" << i << ", size " << pk.seqs[i].size << ", reps " << mc_reps.size() << endl;
//         cout << pk.sequences.value()[i] << endl;
//         for(auto r: mc_reps){
//             Pos pos = r.pos;
//             StratumId rid = r.stratum_id;
//             auto str = pk.sequences.value()[i].substr(pos, pk.stratums[rid].size);
//             _print_rep(rid, str, 1, pk, printed_rep);
//         }
//         cout << endl;
//     }
// }

// void print_bare_prokrustean(Prokrustean pk){
//     for(uint64_t i=0; i<pk.seqs.size(); i++){
//         auto mc_reps = pk.seqs[i].regions;
//         cout << "seq" << i << ", size " << pk.seqs[i].size << ", reps " << mc_reps.size() << endl;
//         for(auto r: mc_reps){
//             Pos pos = r.pos;
//             StratumId rid = r.stratum_id;
//             cout << "(" << "pos: " << pos << ", " << "size: " << pk.stratums[rid].size<< ", R" << rid << ")" << " ";
//         }
//         cout << endl;
//     }
//     for(uint64_t i=0; i<pk.stratums.size(); i++){
//         auto mc_reps = pk.stratums[i].regions;
//         cout << "rep" << i << ", size " << pk.stratums[i].size << ", reps " << mc_reps.size() << endl;
//         for(auto r: mc_reps){
//             Pos pos = r.pos;
//             StratumId rid = r.stratum_id;
//             cout << "(" << pos << ", " << pk.stratums[rid].size<< ", R" << rid << ")" << " ";
//         }
//         cout << endl;
//     }
// }

// void print_single_rep_to_single_rep_relationships(Prokrustean pk){
//     int cnt = 0;
//     for(auto mc: pk.stratums){
//         if(mc.regions.size()==1){
//             auto unique_repid = mc.regions[0].stratum_id;
//             if(pk.stratums[unique_repid].regions.size()>0){
//                 cnt++;
//             }
//         }
//     }
//     cout << "-- print_single_rep_to_single_rep_relationships --" << endl;
//     cout << "total: " << cnt << endl;
// }
#endif