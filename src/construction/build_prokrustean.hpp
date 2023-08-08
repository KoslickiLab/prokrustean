
#include "../fm_index/tree.hpp"
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include <algorithm>

using namespace std;

#ifndef CONSTRUCTION_PROKRUSTEAN_HPP_
#define CONSTRUCTION_PROKRUSTEAN_HPP_

struct MaxRepeat {
public:
    // 
    uint64_t size;
    vector<uint64_t> repr_sa_indexes;
};

class AbsReprSuffPositions {
public:
    // 
    virtual uint64_t count() = 0;

    // 
    virtual void set(uint64_t sa_idx) = 0;

    // 
    virtual bool get(uint64_t sa_idx) = 0;

    // 
    virtual uint64_t rank(uint64_t sa_idx) = 0;
};

class AbsReprSuff {
public:
    virtual void initialize(uint64_t repr_sa_cnt) = 0;

    //
    virtual void set(uint64_t sa_index, RepId rep_id) = 0;

    // 
    virtual vector<RepId> get(uint64_t sa_index, bool remove=true) = 0;
};

class AbsProkrusteanPrep {
public:
    virtual void initialize(uint64_t seq_cnt, uint64_t rep_cnt) = 0;

    virtual void set_mc_rep(SeqId id, uint64_t size, vector<tuple<uint64_t, RepId>> mc_rep) = 0;

    virtual void set_mc_rep(RepId id, uint64_t size, vector<tuple<uint64_t, RepId>> mc_rep) = 0;

    virtual vector<MinCover> get_seq_MCs() = 0;

    virtual vector<MinCover> get_rep_MCs() = 0;
};

void build_prokrustean(FmIndex &fm_idx, 
                       AbsReprSuffPositions &repr_sa_pos,
                       AbsReprSuff &repr_sa,
                       AbsProkrusteanPrep &pk_prep){
    // step1 collect representative suffix array
    Interval root = get_root(fm_idx);
    // parallelize by subtree
    vector<MaxRepeat> repeats = _collect_max_repeat_and_repr_sa(root, fm_idx, repr_sa_pos);

    // step2 flip structure
    repr_sa.initialize(repr_sa_pos.count());
    // parallelize by repeats
    _build_repr_sa(repeats, repr_sa);

    // step3 build structure
    // parallelize by sequences
    for(uint64_t i; i < fm_idx.seq_cnt(); i++){
        _build_mc_rep(i, fm_idx, repr_sa_pos, repr_sa, pk_prep);
    }
}

vector<MaxRepeat> _collect_max_repeat_and_repr_sa(Interval &curr_interval, FmIndex &fm_idx, AbsReprSuffPositions &repr_sa_pos){
    vector<MaxRepeat> reps;
    left_extension extension = left_extend(fm_idx, curr_interval);
    /* find rep */
    bool is_rep = extension.left_maximal() && curr_interval.right_maximal();
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
        reps.push_back({curr_interval.depth, sa_indexes});
        for(auto sa_idx :sa_indexes){
            repr_sa_pos.set(sa_idx);
        }
    }
    // navigate
    for(auto c_interval:extension.intervals){
        if(c_interval.count()>1){
            vector<MaxRepeat> c_reps = _collect_max_repeat_and_repr_sa(c_interval, fm_idx, repr_sa_pos);
            reps.insert(reps.end(), c_reps.begin(), c_reps.end());
        }
    }
    return reps;
}

void _build_repr_sa(vector<MaxRepeat> repeats, AbsReprSuff &repr_sa){
    uint64_t rep_id = 0;
    for(auto rep: repeats){
        for(auto sa_idx: rep.repr_sa_indexes){
            repr_sa.set(sa_idx, rep_id);
            rep_id++;
        }
    }
}

void _build_mc_rep(SeqId id, 
                   FmIndex &fm_idx, 
                   AbsReprSuffPositions &repr_sa_pos,
                   AbsReprSuff &repr_sa,
                   AbsProkrusteanPrep &pk_prep){
    
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