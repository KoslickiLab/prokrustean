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

optional<MaximalRepeatAnnotation> get_repeat_annotations(SuffixArrayNodeExtension &ext){
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
        SuffixArrayIdx min_idx = repr_sa[0];
        for(auto i: repr_sa) 
        min_idx = i < min_idx? i: min_idx;
        
        MaximalRepeatAnnotation rep = {ext.node.depth, repr_sa, min_idx};
        return rep;
    } else {
        return nullopt;
    }
}

SequenceAnnotation get_sequence_annotations(SeqId seq_id, 
                                            FmIndex &fm_idx, 
                                            ReprSuffixAnnotation &repr_sa, 
                                            vector<MaximalRepeatAnnotation> &rep_annots, 
                                            bool recover_sequences=false){
    // 
    uint64_t L = seq_id;
    uint64_t F = fm_idx.LF(L);
    string seq;
    uint64_t seq_length;
    /* get seq length, recover seq */
    while(F >= fm_idx.seq_cnt()){
        if(recover_sequences){
            seq = fm_idx.get_character(L) + seq;
        }
        L = F;
        F = fm_idx.LF(L);
        seq_length++;
    }

    vector<PositionAnnotation> locs;
    L = seq_id;
    F = fm_idx.LF(L);
    uint64_t pos = seq_length-1;
    /* plot occurrences */
    while(pos >= 0){
        auto v = repr_sa.get_repeats(F);
        if (v.has_value()){
            vector<RepId> rep_ids = v.value();
            vector<MaximalRepeatAnnotation> reps;
            for(auto id: rep_ids){
                reps.push_back(rep_annots[id]);
            }
            locs.push_back({pos, F, reps, rep_ids});
        }
        L = F;
        F = fm_idx.LF(L);
        pos--;
    }
    return {seq_id, seq_length, locs, seq};
}


vector<MinCover> get_min_covers(SequenceAnnotation &seq_annot){
    /* 
    input: annotations (grouped by positions, sorted by asc size.)
    variables:
    * current rep reference for each pos: CURR(pos) -> rep. 
    * an ancestor stack and its topmost item (parent pos, parent rep).
    * active pos
    operation: 
    * left(pos),right(pos)
    * right(pos, rep) = (right(pos), CURR(right(pos)))
    * left(pos, rep) = (left(pos), CURR(left(pos)))
    algo:
    make mc_rep(seq)
    set pos <- first pos, first rep
    set curr <- all first reps
    iterate below until break
    if moved_to_new_pos,
        make mc_rep(rep)
        if(rep not first){
            add prev(rep) to mc_reps(rep)
        }
    else
        moved_to_new_pos = true
    (pos, rep) <- (active pos, curr(active pos))
    if (pos, rep) includes right(pos, rep),
        stack.push(pos, rep)
        (parent pos, parent rep) <- (pos, rep) 
        active pos <- right(pos).
        go to next iteration.
    else if (parent pos, parent rep) includes right(pos, rep),
        active pos <- right(pos).
        go to next iteration.
    
    if (parent pos, parent rep) != null,
        if not last rep,
            CURR(pos) <- next(rep)
            if (parent pos, parent rep) includes (pos, next(rep)),
                go to next iteration
        add rep to (parent rep).mc_reps
        active pos <- left(pos)
        moved_to_new_pos = false
        if (pos, rep) == (parent pos, parent rep)
            parent pos <- stack.pop() (null if empty)
    else
        if not last rep,
            CURR(pos) <- next(rep)
        else
            CURR(pos) <- null
            add rep to (seq).mc_reps
            if pos is last
                break
            else
                active pos <- right(pos).
    go to next iteration
    */
    
    MinCover seq_mc;
    seq_mc.id = seq_annot.id;
    seq_mc.size = seq_annot.size;
    seq_mc.is_rep = false;
    vector<MinCover> mcs({seq_mc});

    if(seq_annot.positions.size()==0)
    return mcs;

    // position is already sorted
    for(auto p: seq_annot.positions){
        sort(p.reps.begin(), p.reps.end(), 
        [](MaximalRepeatAnnotation rep1, MaximalRepeatAnnotation rep2) {return rep1.size < rep2.size; });
    }

    // index is the position's index, value is the reps's index
    vector<uint64_t> curr_rep_indexes(seq_annot.positions.size());
    stack<tuple<uint64_t, MinCover>> ancestors;
    optional<tuple<uint64_t, MinCover>> parent;
    uint64_t active_pos_idx = 0;
    bool moved_to_new_pos_at_last_iter = true;
    while(true){
        PositionAnnotation active_pos = seq_annot.positions[active_pos_idx];
        uint64_t active_rep_idx = curr_rep_indexes[active_pos_idx];
        MaximalRepeatAnnotation active_rep = active_pos.reps[active_rep_idx];
        optional<MinCover> new_mc;
        if(moved_to_new_pos_at_last_iter){
            //if the suffix array index is the same as the first repr suffix array index of the repeat.
            bool rep_is_the_primary = active_rep.first_repr_idx == active_pos.sa_idx;
            if(rep_is_the_primary){
                MinCover rep_mc;
                rep_mc.id = active_pos.rep_ids[active_rep_idx];
                rep_mc.size = active_rep.size;
                rep_mc.is_rep = true;
                if(active_rep_idx>0){
                    auto prev_rep_id = active_pos.rep_ids[active_rep_idx-1];
                    Pos relative_pos = 0;
                    rep_mc.mc_reps.push_back(make_tuple(relative_pos, prev_rep_id));
                }
                new_mc = rep_mc;
                mcs.push_back(rep_mc);
            }
        } else {
            moved_to_new_pos_at_last_iter = true;
        }

        if(active_pos_idx+1 < seq_annot.positions.size()){
            auto right_pos = seq_annot.positions[active_pos_idx+1];
            auto right_rep = right_pos.reps[curr_rep_indexes[active_pos_idx+1]];
            // active includes right
            if(right_pos.pos+right_rep.size <= active_pos.pos+active_rep.size){
                //if the repeat becomes a parent, it should be primary.
                assert(new_mc.has_value());
                parent = make_tuple(active_pos_idx, new_mc.value());
                ancestors.push(parent.value());
                active_pos_idx++;
                continue;
            } // parent includes right
            else if(parent.has_value()) {
                uint64_t parent_pos_idx = get<0>(parent.value());
                MinCover parent_mc = get<1>(parent.value());
                PositionAnnotation parent_pos = seq_annot.positions[parent_pos_idx];
                MaximalRepeatAnnotation parent_rep = parent_pos.reps[curr_rep_indexes[parent_pos_idx]];
                if(right_pos.pos+right_rep.size < parent_pos.pos+parent_rep.size){
                    active_pos_idx++;
                    continue;
                }
            }
        }

        if(parent.has_value()){
            uint64_t parent_pos_idx = get<0>(parent.value());
            MinCover parent_mc = get<1>(parent.value());
            PositionAnnotation parent_pos = seq_annot.positions[parent_pos_idx];
            MaximalRepeatAnnotation parent_rep = parent_pos.reps[curr_rep_indexes[parent_pos_idx]];
            //not last rep
            if(active_rep_idx+1 < active_pos.reps.size()){
                // move to the next rep
                curr_rep_indexes[active_pos_idx] = active_rep_idx+1;
                auto next_rep = active_pos.reps[active_rep_idx+1];
                // parent includes the next larger rep
                if(active_pos.pos+next_rep.size < parent_pos.pos+parent_rep.size){
                    continue;
                }
            }
            // insert into parent
            Pos relative_pos = active_pos.pos-parent_pos.pos;
            parent_mc.mc_reps.push_back(make_tuple(relative_pos, active_pos.rep_ids[active_rep_idx]));

            //active_pos_idx cannot be 0 because parent would not exists
            assert(active_pos_idx>0);
            active_pos_idx--;
            
            moved_to_new_pos_at_last_iter = false;
            if(active_pos_idx == parent_pos_idx){
                if(ancestors.empty()) 
                parent = nullopt;
                else 
                parent = ancestors.top();
            }
        } else { // no parent
            curr_rep_indexes[active_pos_idx] = active_rep_idx + 1;

            if(active_rep_idx == active_pos.reps.size()){ 
                Pos relative_pos = active_pos.pos;
                seq_mc.mc_reps.push_back(make_tuple(relative_pos, active_pos.rep_ids[active_rep_idx]));
                if(active_pos_idx == seq_annot.positions.size()){
                    break;
                } else {
                    active_pos_idx++;
                }
            }
        }
    }
}

#endif