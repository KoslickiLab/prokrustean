#include "models.hpp"
#include "../fm_index/tree.hpp"
#include <algorithm>
#include <stack>

using namespace std;

#ifndef CONSTRUCTION_ALGO_PROCEDURES_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_HPP_

tuple<vector<CharId>, vector<CharId>> decide_repr_sa_extensions(int char_max, vector<tuple<CharId, CharId>> distinct_extensions){
    vector<tuple<int, CharId>> left_pair_info(char_max);
    vector<tuple<int, CharId>> right_pair_info(char_max);
    vector<CharId> left_repr_extensions;
    vector<CharId> right_repr_extensions;
    for(auto pair:distinct_extensions){
        auto left = get<0>(pair);
        auto right = get<1>(pair);
        left_pair_info[left] = make_tuple(get<0>(left_pair_info[left])+1, right);
        right_pair_info[right] = make_tuple(get<0>(right_pair_info[right])+1, left);
    }
    for(int i=1/*skip term*/; i<char_max; i++){ 
        int left_cnt = get<0>(left_pair_info[i]);
        CharId r = get<1>(left_pair_info[i]);
        if(left_cnt==0) 
        continue;
        //left exclusive but not bi-exclusive
        if(left_cnt==1 && r!=0 && get<0>(right_pair_info[r])>1)
        continue;
        left_repr_extensions.push_back(i);
    }
    for(int i=1/*skip term*/; i<char_max; i++){ 
        int right_cnt = get<0>(right_pair_info[i]);
        CharId l = get<1>(right_pair_info[i]);
        if(right_cnt==0) 
        continue;
        //right exclusive but not bi-exclusive
        if(right_cnt==1 && l!=0 && get<0>(left_pair_info[l])>1) 
        continue;
        right_repr_extensions.push_back(i);
    }
    return make_tuple(left_repr_extensions, right_repr_extensions);
}


// vector<SuffixArrayIdx> _decide_repr_sa_extensions(int char_max, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
//     vector<SuffixArrayIdx> sa_indexes;
//     vector<tuple<CharId, SuffixArrayIdx>> left_paired_1st(char_max);
//     vector<tuple<CharId, SuffixArrayIdx>> right_paired_1st(char_max);
//     vector<int> left_paired_cnt(char_max);
//     vector<int> right_paired_cnt(char_max);
//     for(auto pair:distinct_extensions){
//         auto left = get<0>(pair);
//         auto right = get<1>(pair);
//         auto sa_idx = get<2>(pair);
        
//         if(left_paired_cnt[left]==0 || sa_idx < get<1>(left_paired_1st[left])){
//             left_paired_1st[left]=make_tuple(right, sa_idx);
//         }
//         left_paired_cnt[left]++;

//         if(right_paired_cnt[right]==0 || sa_idx < get<1>(right_paired_1st[right])){
//             right_paired_1st[right]=make_tuple(left, sa_idx);
//         }
//         right_paired_cnt[right]++;
//     }
//     for(int i=1/*skip term*/; i<char_max; i++){ 
//         if(left_paired_cnt[i]==0) 
//         continue;
//         //left exclusive but not bi-exclusive
//         CharId r = get<0>(left_paired_1st[i]);
//         if(left_paired_cnt[i]==1 && r!=0 && right_paired_cnt[r]>1) 
//         continue;
//         uint64_t sa_idx = get<1>(left_paired_1st[i]);
//         sa_indexes.push_back(sa_idx);
        
//         cout << "left char: ";
//         switch (i)
//         {
//         case 1: cout << "A"; break;
//         case 2: cout << "C"; break;
//         case 3: cout << "G"; break;
//         case 4: cout << "T"; break;
//         default: break;
//         }
//         cout << ", " << sa_idx << endl;
//     }
//     for(int i=1/*skip term*/; i<char_max; i++){ 
//         if(right_paired_cnt[i]==0) 
//         continue;
//         //left exclusive but not bi-exclusive
//         CharId l = get<0>(right_paired_1st[i]);
//         if(left_paired_cnt[i]==1 && l!=0) 
//         continue;
//         uint64_t sa_idx = get<1>(right_paired_1st[i]);
//         sa_indexes.push_back(sa_idx);

//         cout << "right char: ";
//         switch (i)
//         {
//         case 1: cout << "A"; break;
//         case 2: cout << "C"; break;
//         case 3: cout << "G"; break;
//         case 4: cout << "T"; break;
//         default: break;
//         }
//         cout << ", " << sa_idx << endl;
//     }
//     return sa_indexes;
// }

optional<MaximalRepeatAnnotation> get_rep_annot(SuffixArrayNodeExtension &ext){
    if(ext.left_maximal() && ext.node.right_maximal()){
        tuple<vector<CharId>, vector<CharId>> repr_extensions = decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions());
        // Remove possible duplications by set. Duplications cannot be predicted in advance.
        set<SuffixArrayIdx> uniq_repr_sa;
        for(auto l: get<0>(repr_extensions)){
            uniq_repr_sa.insert(ext.first_l(l));
        }
        for(auto r: get<1>(repr_extensions)){
            uniq_repr_sa.insert(ext.first_r(r));
        }
        vector<SuffixArrayIdx> repr_sa(uniq_repr_sa.begin(), uniq_repr_sa.end());
        MaximalRepeatAnnotation rep = {ext.node.depth, repr_sa};
        return rep;
    } else {
        return nullopt;
    }
}

SequenceAnnotation get_seq_annot(SeqId id, FmIndex &fm_idx, ReprSuffixAnnotation &repr_sa){
    stack<tuple<Pos, vector<RepId>>> temp_annotations;

    uint64_t L = id;
    uint64_t F = fm_idx.LF(L);
    uint64_t rev_pos = 0;
    /* plot occurrences */
    while(F >= fm_idx.seq_cnt()){
        optional<vector<RepId>> reps = repr_sa.get_repeats(F);
        if (reps.has_value()){
            temp_annotations.push(make_tuple(rev_pos, reps.value()));
        }
        L = F;
        F = fm_idx.LF(L);
        rev_pos++;
    }

    uint64_t seq_len = rev_pos;
    vector<tuple<Pos, vector<RepId>>> annotations;
    while(!temp_annotations.empty()){
        auto item = temp_annotations.top();
        Pos rpos = get<0>(item);
        vector<RepId> reps = get<1>(item);
        annotations.push_back(make_tuple(seq_len-1-rpos, reps));
    }
    return {seq_len, annotations};
}

vector<MinCover> get_min_covers(SequenceAnnotation annotation){

}

#endif