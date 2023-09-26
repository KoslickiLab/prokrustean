#ifndef CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#include "models.hpp"
#include "../fm_index/tree_new.hpp"
#include <algorithm>
#include <stack>
#include <tuple>

using namespace std;

void decide_repr_sa_new(SuffixArrayNodeExtension_NEW &ext){
    // initialize cnts
    for(int i=0; i<ext.characters_cnt; i++){
        ext.repr_space.left_cnts[i]=0;
        ext.repr_space.right_cnts[i]=0;
    }
    // for each letter
    for(int c=0; c < ext.characters_cnt; c++){
        // c_node not exists
        if(!ext.c_nodes_open[c]) continue;

        //for each right a
        for(int a=0; a<ext.characters_cnt; a++){
            if(ext.c_nodes[c].firsts[a]<ext.c_nodes[c].firsts[a+1]){
                ext.repr_space.left_cnts[c]++;
                ext.repr_space.left_paired_a_char[c]=a;
                ext.repr_space.right_cnts[a]++;
                ext.repr_space.right_paired_c_char[a]=c;
            }
        }
    }
    /* find repr characters */
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        // cout << "left: " << i << ", left cnt" << left_cnt << ", right info: " << (int)r << endl;
        if(ext.repr_space.left_cnts[i]==0){
            ext.repr_space.left_repr[i]=false;
        } //left exclusive but not bi-exclusive
        else if(ext.repr_space.left_cnts[i]==1 
        && ext.repr_space.left_paired_a_char[i]!=0 
        && ext.repr_space.right_cnts[ext.repr_space.left_paired_a_char[i]]>1){
            ext.repr_space.left_repr[i]=false;
        } else {
            ext.repr_space.left_repr[i]=true;
        }
    }
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        if(ext.repr_space.right_cnts[i]==0)
        {
            ext.repr_space.right_repr[i]=false;
        } //right exclusive but not bi-exclusive
        else if(ext.repr_space.right_cnts[i]==1 
        && ext.repr_space.right_paired_c_char[i]!=0 
        && ext.repr_space.left_cnts[ext.repr_space.right_paired_c_char[i]]>1) {
            ext.repr_space.right_repr[i]=false;
        } else {
            ext.repr_space.right_repr[i]=true;
        }
    }
}


void report_repr_sa(FmIndex &index, SuffixArrayNodeExtension_NEW &ext, vector<MaximalRepeatAnnotation> &outs){
    decide_repr_sa_new(ext);

    ext.repr_space.uniq_repr_sa.clear();
    for(int c=0; c<ext.characters_cnt; c++){
        // the first suffix index of Wa form is explored by default.
        if(ext.repr_space.right_repr[c]){
            ext.repr_space.uniq_repr_sa.push_back(ext.first_r(c));
        }
        // the first suffix index of cW should be found and checked for duplication
        if(ext.repr_space.left_repr[c]){
            //revert the rank process
            SuffixArrayIdx sa_idx=0;
            CharId a;
            for(int i=0; i<ext.characters_cnt; i++){
                if(ext.c_nodes[c].firsts[i]<ext.c_nodes[c].firsts[i+1]) {
                    // first suffix array index of cW 
                    sa_idx=ext.c_nodes[c].firsts[i];
                    a=i;
                    break;
                }
            }
            assert(sa_idx!=0);
            // auto start2 = std::chrono::steady_clock::now();
            uint64_t rank = sa_idx - index.C[c] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
            // ext.any_measure[6]+=(std::chrono::steady_clock::now()-start2).count();
            
            // here, check if the suffix array picked by the right letter a of cWa is duplicated.
            // this cannot be checked in advance..
            if(ext.repr_space.right_repr[a] && prev_sa_idx==ext.first_r(a)){
            } else {
                ext.repr_space.uniq_repr_sa.push_back(prev_sa_idx);
            }
        }
    }
    
    /*important. If termination is placed in both ends, all are representative. 
    decide_repr_sa_extensions marks it by simply including 0 in left character. 
    Then in here we should collect all suffixes where the form is #W#.
    */
    for(auto sa_idx: ext.both_ext_terms){
        ext.repr_space.uniq_repr_sa.push_back(sa_idx);
    }

    assert(ext.repr_space.uniq_repr_sa.size()>0);
    // SuffixArrayIdx min_idx2 = ext.repr_sa_workspace.uniq_repr_sa[0];
    SuffixArrayIdx min_idx = ext.repr_space.uniq_repr_sa[0];
    for(auto idx: ext.repr_space.uniq_repr_sa){
        // min_idx2 = ext.repr_sa_workspace.uniq_repr_sa[i] < min_idx2? ext.repr_sa_workspace.uniq_repr_sa[i]: min_idx2;
        min_idx = idx < min_idx? idx: min_idx;
    }

    outs.push_back({ext.node.depth, ext.repr_space.uniq_repr_sa, min_idx});
}


SequenceAnnotation get_sequence_annotations_new(SeqId seq_id, 
                                            FmIndex &fm_idx, 
                                            ReprSuffixAnnotation &repr_sa, 
                                            vector<MaximalRepeatAnnotation> &rep_annots, 
                                            bool recover_sequences=false){
    // 
    uint64_t L = seq_id;
    uint64_t F = fm_idx.LF(L);
    string seq;
    uint64_t seq_length=0;
    uint64_t repr_cnt=0;
    /* get seq length, recover seq */
    while(F >= fm_idx.seq_cnt()){
        if(repr_sa.exists(F)){
            repr_cnt++;
        }
        if(recover_sequences){
            seq = fm_idx.get_character(L) + seq;
        }
        L = F;
        F = fm_idx.LF(L);
        seq_length++;
    }

    uint32_t pos = seq_length-1;
    uint32_t loc_idx = repr_cnt-1;
    vector<PositionAnnotation> locs(repr_cnt);
    L = seq_id;
    F = fm_idx.LF(L);
    /* plot occurrences */
    while(F >= fm_idx.seq_cnt()){
        auto v = repr_sa.get_repeats(F);
        if (v.has_value()){
            vector<StratumId> rep_ids = v.value();
            vector<MaximalRepeatAnnotation> reps;
            for(auto id: rep_ids){
                reps.push_back(rep_annots[id]);
            }
            sort(rep_ids.begin(), rep_ids.end(),
            [&rep_annots](uint64_t r1, uint64_t r2) {return rep_annots[r1].size < rep_annots[r2].size; });
            sort(reps.begin(), reps.end(),
            [](MaximalRepeatAnnotation r1, MaximalRepeatAnnotation r2) {return r1.size < r2.size; });

            locs[loc_idx]={pos, F, reps, rep_ids};
            loc_idx--;
        }
        L = F;
        F = fm_idx.LF(L);
        pos--;
    }
    return {seq_id, seq_length, locs, seq};
}


#endif