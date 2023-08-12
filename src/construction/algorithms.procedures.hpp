#include "models.hpp"
#include "../fm_index/tree.hpp"
#include <algorithm>
#include <stack>

using namespace std;

#ifndef CONSTRUCTION_ALGO_PROCEDURES_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_HPP_

vector<SuffixArrayIdx> decide_repr_sa_extensions(int char_max, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
    vector<SuffixArrayIdx> sa_indexes;
    vector<tuple<CharId, SuffixArrayIdx>> left_paired_1st(char_max);
    vector<tuple<CharId, SuffixArrayIdx>> right_paired_1st(char_max);
    vector<int> left_paired_cnt(char_max);
    vector<int> right_paired_cnt(char_max);
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
    for(int i=1/*skip term*/; i<char_max; i++){ 
        if(left_paired_cnt[i]==0) 
        continue;
        //left exclusive but not bi-exclusive
        CharId r = get<0>(left_paired_1st[i]);
        if(left_paired_cnt[i]==1 && r!=0 && right_paired_cnt[r]>1) 
        continue;
        uint64_t sa_idx = get<1>(left_paired_1st[i]);
        sa_indexes.push_back(sa_idx);
        
        cout << "left char: ";
        switch (i)
        {
        case 1: cout << "A"; break;
        case 2: cout << "C"; break;
        case 3: cout << "G"; break;
        case 4: cout << "T"; break;
        default: break;
        }
        cout << ", " << sa_idx << endl;
    }
    for(int i=1/*skip term*/; i<char_max; i++){ 
        if(right_paired_cnt[i]==0) 
        continue;
        //left exclusive but not bi-exclusive
        CharId l = get<0>(right_paired_1st[i]);
        if(left_paired_cnt[i]==1 && l!=0) 
        continue;
        uint64_t sa_idx = get<1>(right_paired_1st[i]);
        sa_indexes.push_back(sa_idx);

        cout << "right char: ";
        switch (i)
        {
        case 1: cout << "A"; break;
        case 2: cout << "C"; break;
        case 3: cout << "G"; break;
        case 4: cout << "T"; break;
        default: break;
        }
        cout << ", " << sa_idx << endl;
    }
    return sa_indexes;
}

optional<MaximalRepeatAnnotation> get_rep_annot(SuffixArrayNodeExtension &ext){
    if(ext.left_maximal() && ext.node.right_maximal()){
        cout << "max repeat of " << ext.node.depth << endl;
        vector<uint64_t> repr_sa = decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions());
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