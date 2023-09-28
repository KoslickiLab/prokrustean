#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "base.hpp"

using namespace std;

struct Region {
    Pos from;
    Pos to;
    
    Region(){}
    Region(Pos from, Pos to): from(from), to(to){}

    uint64_t size(){
        assert(from<to);
        return to-from;
    }
};

struct StratifiedRegion: Region {
    StratumId stratum_id;

    StratifiedRegion(){}
    StratifiedRegion(Pos from, Pos to, StratumId stratum_id): Region(from, to), stratum_id(stratum_id){}
};

struct ReflectedRegion: Region {

    ReflectedRegion(){}
    ReflectedRegion(Pos from, Pos to): Region(from, to){}
};

struct StratifiedData {
    StratumId stratum_id;
    Pos pos;
    
    StratifiedData(){}
    StratifiedData(Pos pos, StratumId stratum_id): pos(pos), stratum_id(stratum_id)
    {}
};

struct Stratum {
    StratumSize size;
    vector<StratifiedRegion> regions;
    tuple<SeqId, Pos> example_occ;

    void set_occ(SeqId seq_id, Pos pos){
        example_occ=make_tuple(seq_id, pos);
    }
};

struct Sequence {
    SequenceSize size;
    vector<StratifiedRegion> regions;
};

struct Prokrustean {
    /* data structure is succinct as possible*/
    vector<SequenceSize> sequences__size;
    vector<StratifiedData*> sequences__region;
    vector<uint8_t> sequences__region_cnt;
    
    vector<StratumSize> stratums__size;
    vector<StratifiedData*> stratums__region;
    vector<uint8_t> stratums__region_cnt;

    uint64_t sequence_count(){
        return this->sequences__size.size();
    }

    uint64_t stratum_count(){
        return this->stratums__size.size();
    }

    void get_stratum_example_occ(vector<tuple<SeqId, Pos>> &stratum_pos){
        stratum_pos.resize(stratum_count());
        vector<bool> visits(stratum_count());
        stack<Stratum> stratum_stack;
        for(int i=0; i<sequence_count(); i++){
            for(auto &rgn: get_sequence(i).regions){
                auto stratum = get_stratum(rgn.stratum_id);
                stratum.set_occ(i, rgn.from);
                stratum_pos[rgn.stratum_id]=stratum.example_occ;
                visits[rgn.stratum_id]=true;

                stratum_stack.push(stratum);
                while(!stratum_stack.empty()){
                    auto stratum=stratum_stack.top();
                    stratum_stack.pop();
                    for(auto &c_rgn: stratum.regions){
                        if(visits[c_rgn.stratum_id]) continue;
                        SeqId seq_id = get<0>(stratum.example_occ);
                        Pos rel_pos = get<1>(stratum.example_occ)+c_rgn.from;
                        auto c_stratum = get_stratum(c_rgn.stratum_id);
                        c_stratum.set_occ(seq_id, rel_pos);
                        stratum_pos[c_rgn.stratum_id]=c_stratum.example_occ;
                        visits[c_rgn.stratum_id]=true;
                        stratum_stack.push(c_stratum);
                    }
                }
            }
        }
        for(auto &el: stratum_pos){
            cout << get<0>(el) << ", " << get<1>(el) << endl;
        }

    }

    void get_spectrum(std::variant<Sequence, Stratum> v, int k,vector<std::variant<StratifiedRegion, ReflectedRegion>> &output){
        uint64_t size;
        vector<StratifiedRegion> regions;
        if (std::holds_alternative<Sequence>(v)) {
            size=std::get<Sequence>(v).size;
            regions=std::get<Sequence>(v).regions;
        } else {
            size=std::get<Stratum>(v).size;
            regions=std::get<Stratum>(v).regions;
        }
        
        int cnt=regions.size();
        // single reflected 
        if(cnt==0){
            output.push_back(ReflectedRegion(0, size));
            return;
        }

        // already sorted
        for(int i=0; i<cnt; i++){
            if(i==0 && regions[i].from>0){
                output.push_back(ReflectedRegion(0, regions[i].from+(k-1)));  
            }

            output.push_back(regions[i]);

            if(i<cnt-1){
                if(regions[i].to - regions[i+1].from >= k-1){
                    // stratified are overlapped too much that no reflected region exists
                } else {
                    output.push_back(ReflectedRegion(regions[i].to-(k-1), regions[i+1].from+(k-1)));      
                }
            }

            if(i==cnt-1 && regions[i].to<size){
                output.push_back(ReflectedRegion(regions[i].to-(k-1), size));  
            }
        }
    } 

    Sequence get_sequence(SeqId id){
        auto sequence=Sequence();
        auto rgn_cnt=sequences__region_cnt[id];
        sequence.size=sequences__size[id];
        sequence.regions.resize(rgn_cnt);
        while(rgn_cnt>0){
            rgn_cnt--;
            StratifiedData data=sequences__region[id][rgn_cnt];
            auto rgn=StratifiedRegion(data.pos, data.pos +this->stratums__size[data.stratum_id], data.stratum_id);
            sequence.regions[rgn_cnt]=rgn;
            assert(0<=rgn.from && rgn.from < sequence.size && 0<rgn.to && rgn.to <= sequence.size);
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
            auto data=stratums__region[id][rgn_cnt];
            auto rgn=StratifiedRegion(data.pos, data.pos +this->stratums__size[data.stratum_id], data.stratum_id);
            stratum.regions[rgn_cnt]=rgn;
            assert(0<=rgn.from && rgn.from < stratum.size && 0<rgn.to && rgn.to <= stratum.size);
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