
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
    virtual vector<RepId> pop(uint64_t sa_idx) = 0;

    virtual void push(RepId rep_idx, vector<uint64_t> sa_indexes) = 0;
};

void build_prokrustean(FmIndex &fm_idx, AbstractProkrusteanPrep &pk_prep, AbstractReprSA &repr_sa){
    Interval root = get_root(fm_idx);
    _collect_representative_suffix_array(root, fm_idx, pk_prep, repr_sa);

    // for(auto idx: fm_idx.seq_cnt()){
    //     _build_mc_rep(idx, fm_idx, pk, repr_sa);
    // }
    
}

void _collect_representative_suffix_array(Interval &curr_interval, FmIndex &fm_idx, AbstractProkrusteanPrep &pk, AbstractReprSA &repr_sa){
    left_extension extension = left_extend(fm_idx, curr_interval);
    /* find rep */
    bool is_rep = extension.left_maximal();
    if(is_rep){
        /* find repr sa */
        vector<uint64_t> sa_indexes;
        // search not extensible
        int char_cnt = fm_idx.characters.size();
        vector<CharId> left_paired_example(char_cnt);
        vector<CharId> right_paired_example(char_cnt);
        vector<int> left_paired_cnt(char_cnt);
        vector<int> right_paired_cnt(char_cnt);
        for(auto pair:extension.distinct_extensions()){
            left_paired_cnt[get<0>(pair)]++;
            left_paired_example[get<0>(pair)]=get<1>(pair);
            right_paired_cnt[get<1>(pair)]++;
            right_paired_example[get<1>(pair)]=get<0>(pair);
        }
        for(int i=1/*skip term*/; i<char_cnt-1; i++){
            if(left_paired_cnt[i]==0) 
            continue;
            //left exclusive but not bi-exclusive
            if(left_paired_cnt[i]==1
            && left_paired_example[i]!=0
            && right_paired_cnt[left_paired_example[i]]>1) 
            continue;
            sa_indexes.push_back(extension.first_l(i));
        }
        for(int i=1/*skip term*/; i<char_cnt-1; i++){
            if(right_paired_cnt[i]==0) 
            continue;
            //even if bi-exclusive, it must have been covered above.
            if(right_paired_cnt[i]==1
            && right_paired_example[i]!=0) 
            continue;
            sa_indexes.push_back(extension.first_r(i));
        }

        uint64_t rep_idx = pk.register_rep(curr_interval.depth);
        repr_sa.push(rep_idx, sa_indexes);
    }
    // navigate
    for(auto c_interval:extension.intervals){
        if(c_interval.right_maximal()){
            _collect_representative_suffix_array(c_interval, fm_idx, pk, repr_sa);
        }
    }
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