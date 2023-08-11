#include "models.hpp"
#include "../fm_index/tree.hpp"
#include <algorithm>
#include <stack>

using namespace std;

#ifndef CONSTRUCTION_ALGO_PROCEDURES_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_HPP_

vector<uint64_t> decide_repr_sa_extensions(int char_cnt, vector<tuple<CharId, CharId, uint64_t>> distinct_extensions){
    vector<uint64_t> sa_indexes;
    vector<tuple<CharId, uint64_t>> left_paired_1st(char_cnt);
    vector<tuple<CharId, uint64_t>> right_paired_1st(char_cnt);
    vector<int> left_paired_cnt(char_cnt);
    vector<int> right_paired_cnt(char_cnt);
    for(auto pair:distinct_extensions){
        auto left = get<0>(pair);
        auto right = get<1>(pair);
        auto sa_idx = get<2>(pair);
        
        if(left_paired_cnt[left]==0 || sa_idx < get<1>(left_paired_1st[left])){
            left_paired_1st[left]=make_tuple(right, sa_idx);
        }
        left_paired_cnt[left]++;

        if(right_paired_cnt[right]==0 || sa_idx < get<1>(right_paired_1st[right])){
            right_paired_1st[right]=make_tuple(left, sa_idx);
        }
        right_paired_cnt[right]++;
    }
    for(int i=1/*skip term*/; i<char_cnt-1; i++){ 
        if(left_paired_cnt[i]==0) 
        continue;
        //left exclusive but not bi-exclusive
        CharId r = get<0>(left_paired_1st[i]);
        if(left_paired_cnt[i]==1 && r!=0 && right_paired_cnt[r]>1) 
        continue;
        uint64_t sa_idx = get<1>(left_paired_1st[i]);
        sa_indexes.push_back(sa_idx);
    }
    for(int i=1/*skip term*/; i<char_cnt-1; i++){ 
        if(right_paired_cnt[i]==0) 
        continue;
        //left exclusive but not bi-exclusive
        CharId l = get<0>(right_paired_1st[i]);
        if(left_paired_cnt[i]==1 && l!=0) 
        continue;
        uint64_t sa_idx = get<1>(right_paired_1st[i]);
        sa_indexes.push_back(sa_idx);
    }
}

optional<MaximalRepeatAnnotation> get_rep_annot(SuffixArrayNode &node, FmIndex &fm_idx){
    if(node.left_maximal() && node.interval.right_maximal()){
        vector<uint64_t> repr_sa = decide_repr_sa_extensions(fm_idx.characters.size(), node.distinct_extensions());
        MaximalRepeatAnnotation rep = {node.interval.depth, repr_sa};
        return rep;
    } else {
        return nullopt;
    }
}

// vector<MaximalRepeat> collect_repeats_in_subtree(SuffixArrayNode &root, FmIndex &fm_idx){
//     stack<SuffixArrayNode> stack;
//     vector<MaximalRepeat> reps;
    
//     stack.push(root);
//     while(~stack.empty()){
//         SuffixArrayNode node = stack.top();
//         stack.pop();
//         /* find rep */
//         SuffixArrayNode extension = extend_left(fm_idx, node);
//         bool is_rep = extension.left_maximal() && node.right_maximal();
//         if(is_rep){
//             /* find repr sa */
//             vector<uint64_t> sa_indexes = decide_repr_sa_extensions(fm_idx.characters.size(), extension.distinct_extensions());
//             reps.push_back({node.depth, sa_indexes});
//         }
//         /* navigate */
//         for(auto c_node:extension.c_nodes){
//             if(c_node.count()>1){
//                 stack.push(c_node);
//             }
//         }
//     }
//     return reps;
// }

SequenceAnnotation get_seq_annot(SeqId id, FmIndex &fm_idx, ReprSuffixAnnotation &repr_sa){
    vector<tuple<Pos, vector<RepId>>> annotations;

    uint64_t L = id;
    uint64_t F = fm_idx.LF(L);
    uint64_t rev_pos = 0;
    /* plot occurrences */
    while(F >= fm_idx.seq_cnt()){
        auto reps = repr_sa.get_repeats(F);
        if (reps.has_value()){
            annotations.push_back(make_tuple(rev_pos, reps.value()));
        }
        L = F;
        F = fm_idx.LF(L);
        rev_pos++;
    }
    SequenceAnnotation annot = {rev_pos, annotations};
    return annot;
}

vector<MinCover> get_min_covers(SequenceAnnotation annotation){

}

#endif