#include "build_prokrustean.hpp"
#include <algorithm>

using namespace std;

#ifndef CONSTRUCTION_ABSTRACT_HPP_
#define CONSTRUCTION_ABSTRACT_HPP_

class ProkrusteanProcedures: AbsProkrusteanProcedures{
public:
    vector<MaxRepeat> collect_repeats_from_subtree(Interval &curr_interval, FmIndex &fm_idx, AbsReprSuffPositions &repr_sa_pos){
        vector<MaxRepeat> reps;
        left_extension extension = left_extend(fm_idx, curr_interval);
        /* find rep */
        bool is_rep = extension.left_maximal() && curr_interval.right_maximal();
        if(is_rep){
            /* find repr sa */
            vector<uint64_t> sa_indexes = decide_repr_sa(fm_idx.characters.size(), extension.distinct_extensions());
            reps.push_back({curr_interval.depth, sa_indexes});
            for(auto sa_idx :sa_indexes){
                repr_sa_pos.set(sa_idx);
            }
        }
        // navigate
        for(auto c_interval:extension.intervals){
            if(c_interval.count()>1){
                vector<MaxRepeat> c_reps = collect_repeats_from_subtree(c_interval, fm_idx, repr_sa_pos);
                reps.insert(reps.end(), c_reps.begin(), c_reps.end());
            }
        }
        return reps;
    }

    vector<uint64_t> decide_repr_sa(int char_cnt, vector<tuple<CharId, CharId, uint64_t>> distinct_extensions){
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

    void index_representative_suffix_array(vector<MaxRepeat> repeats, AbsReprSuff &repr_sa){
        uint64_t rep_id = 0;
        for(auto rep: repeats){
            for(auto sa_idx: rep.repr_sa_indexes){
                repr_sa.set(sa_idx, rep_id);
                rep_id++;
            }
        }
    }

    SeqAnnot get_annotation_structure(SeqId id, FmIndex &fm_idx, AbsReprSuffPositions &repr_sa_pos, AbsReprSuff &repr_sa){
        
        vector<tuple<Pos, vector<RepId>>> annotations;

        uint64_t L = id;
        uint64_t F = fm_idx.LF(L);
        uint64_t rev_pos = 0;
        /* plot occurrences */
        while(F >= fm_idx.seq_cnt()){
            if(repr_sa_pos.exists(F)){
                vector<RepId> reps = repr_sa.get(F);
                annotations.push_back(make_tuple(rev_pos, reps));
            }
            L = F;
            F = fm_idx.LF(L);
            rev_pos++;
        }
        SeqAnnot annot = {rev_pos, annotations};
        return annot;
    }

    vector<MinCover> get_min_covers(SeqAnnot annotation){

    }
};

#endif