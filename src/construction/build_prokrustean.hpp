
#include "../fm_index/tree.hpp"
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include <algorithm>

using namespace std;

#ifndef CONSTRUCTION_PROKRUSTEAN_HPP_
#define CONSTRUCTION_PROKRUSTEAN_HPP_

class AbstractProkrusteanPrep {
public:
    virtual void initialize(uint64_t seq_cnt) = 0;

    virtual void set_seq_size(SeqId id, uint64_t size) = 0;

    virtual RepId register_rep(uint64_t size) = 0;

    virtual void set_mc_rep(SeqId id, vector<tuple<uint64_t, RepId>> mc_rep) = 0;

    virtual void set_mc_rep(RepId id, vector<tuple<uint64_t, RepId>> mc_rep) = 0;

    virtual vector<MinCover> get_seq_MCs() = 0;

    virtual vector<MinCover> get_rep_MCs() = 0;
};

class AbstractReprSA {
public:
    // 
    virtual bool exists(uint64_t sa_idx) = 0;

    // returns repeat id list
    virtual vector<uint64_t> pop(uint64_t sa_idx) = 0;

    virtual void push(RepId rep_idx, vector<uint64_t> sa_indexes) = 0;
};

void build_prokrustean(FmIndex &fm_idx, AbstractProkrusteanPrep &pk, AbstractReprSA &repr_sa){
    Interval root = get_root(fm_idx);
    _collect_representative_suffix_array(root, fm_idx, pk, repr_sa);

    // for(auto idx: fm_idx.seq_cnt()){
    //     _build_mc_rep(idx, fm_idx, pk, repr_sa);
    // }
    
}

void _collect_representative_suffix_array(Interval &curr_interval, FmIndex &fm_idx, AbstractProkrusteanPrep &pk, AbstractReprSA &repr_sa){
    left_extension extension = left_extend(fm_idx, curr_interval);
    vector<tuple<char,char>> pairs = extension.collect_extentions();
    bool is_rep = false;
    uint64_t repeat_idx = 0;
    vector<tuple<char,char>> repr_pairs;
    /* find rep */
    // pk.register_rep()
    // repr_sa.push(sa_idx, repeat_idx);
    
    /* find repr sa */
    // for(auto pair: repr_pairs){
    //     uint64_t sa_idx = extension.first(get<0>(pair), get<1>(pair));
    //     repr_sa.insert_repr(sa_idx, repeat_idx);
        
    // }
}

void _build_mc_rep(SeqId id, FmIndex &fm_idx, AbstractProkrusteanPrep &pk, AbstractReprSA &repr_sa){
    
    uint64_t L = id;
    uint64_t F = fm_idx.LF(L);
    /* plot occurrences */
    vector<pair<uint64_t, RepId>> occurrences;
    while(F >= fm_idx.seq_cnt()){
        // sa.push_back(make_tuple(F, seq));
        // if(repr_sa.exists(F)){
        //   repr_sa.pop()
        // }
        L = F;
        F = fm_idx.LF(L);
    }
    
    /* build mc rep */
    // stack
    // while occurrences
    // pk.set_seq_info()
    // pk.set_seq_size
}
#endif