#include "models.hpp"
#include "../fm_index/tree.hpp"
#include <algorithm>
#include <stack>

using namespace std;

#ifndef CONSTRUCTION_ALGO_PROCEDURES_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_HPP_

tuple<vector<CharId>, vector<CharId>> decide_repr_sa_extensions(int char_max, vector<tuple<CharId, CharId>> distinct_extensions){
    // cout << "-- distinct -- " << endl;
    // for(auto ext: distinct_extensions){
    //     cout << (int)get<0>(ext) << ", " << (int)get<1>(ext) << endl;
    // }
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
    // cout << "left: ";
    // for(auto ext: left_repr_extensions){
    //     cout << (int)(ext);
    // }
    // cout<<endl;
    // cout << "right: ";
    // for(auto ext: left_repr_extensions){
    //     cout << (int)(ext);
    // }
    // cout<<endl;
    return make_tuple(left_repr_extensions, right_repr_extensions);
}

optional<MaximalRepeatAnnotation> get_repeat_annotations(SuffixArrayNodeExtension &ext){
    if(ext.left_maximal() && ext.node.right_maximal()){
        tuple<vector<CharId>, vector<CharId>> repr_extensions = decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions());
        // Remove possible duplications by set. Duplications cannot be predicted in advance.
        set<SuffixArrayIdx> uniq_repr_sa;
        vector<SuffixArrayIdx> left_ext_indexes;
        for(auto l: get<0>(repr_extensions)){
            uniq_repr_sa.insert(ext.first_l(l));
            left_ext_indexes.push_back(l);
        }
        for(auto r: get<1>(repr_extensions)){
            uniq_repr_sa.insert(ext.first_r(r));
        }
        vector<SuffixArrayIdx> repr_sa(uniq_repr_sa.begin(), uniq_repr_sa.end()); 
        SuffixArrayIdx min_idx = repr_sa[0];
        for(auto i: repr_sa) 
        min_idx = i < min_idx? i: min_idx;

        MaximalRepeatAnnotation rep = {ext.node.depth, repr_sa, min_idx, left_ext_indexes};
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

    uint64_t pos = seq_length-1;
    uint64_t loc_idx = repr_cnt-1;
    vector<PositionAnnotation> locs(repr_cnt);
    L = seq_id;
    F = fm_idx.LF(L);
    /* plot occurrences */
    while(F >= fm_idx.seq_cnt()){
        auto v = repr_sa.get_repeats(F);
        if (v.has_value()){
            vector<RepId> rep_ids = v.value();
            vector<MaximalRepeatAnnotation> reps;
            for(auto id: rep_ids){
                reps.push_back(rep_annots[id]);
            }
            locs[loc_idx]={pos, F, reps, rep_ids};
            loc_idx--;
        }
        L = F;
        F = fm_idx.LF(L);
        pos--;
    }
    return {seq_id, seq_length, locs, seq};
}

void _print__get_min_cover__progress(vector<uint64_t> &curr_rep_indexes, 
                                    stack<tuple<uint64_t, uint64_t>> &ancestors, 
                                    uint64_t active_pos_idx,
                                    SequenceAnnotation &seq_annot){
    PositionAnnotation active_pos = seq_annot.positions[active_pos_idx];
    uint64_t active_rep_idx = curr_rep_indexes[active_pos_idx];
    MaximalRepeatAnnotation active_rep = active_pos.reps[active_rep_idx];
    optional<tuple<uint64_t, uint64_t>> parent;
    if(!ancestors.empty()){
        parent=ancestors.top();
    } 
    int txt_block_size = 5;
    uint64_t max_rep=0;
    cout<< "---- algo status -----"<<endl;
    cout<< "ancestor("<<ancestors.size()<<")"<<endl;
    if(parent.has_value())
    cout<< "parent pos: "<< seq_annot.positions[get<0>(parent.value())].pos<<endl;
    else
    cout<<"parent empty"<<endl;
    if(active_rep.first_repr_idx != active_pos.sa_idx){
        cout << "active is not primary" << endl;
    }
    cout<< "pos   ";
    for(auto p: seq_annot.positions){
        if(max_rep<p.reps.size())
        max_rep = p.reps.size();
        string str = to_string(p.pos);
        cout << str;
        if(str.size()<txt_block_size)
        cout << std::string(txt_block_size-str.size(), ' ');
    }
    cout<<endl;
    for(int i=0; i<max_rep; i++){
        string label = " - ";
        cout<< label;
        if(label.size()<txt_block_size)
        cout << std::string(txt_block_size-label.size(), ' ');
        for(int p=0; p<seq_annot.positions.size(); p++){
            auto pos = seq_annot.positions[p];
            string str;
            if(i<pos.reps.size()){
                str += to_string(pos.reps[i].size);
                if(i==curr_rep_indexes[p]){
                    str = ">"+str;
                    if(p==active_pos_idx){
                       str = str+"v"; 
                    }
                }
                cout << str;
                if(str.size()<txt_block_size)
                cout << std::string(txt_block_size-str.size(), ' ');
            } else {
                cout << std::string(txt_block_size, ' ');
            }
            
        }
        cout<<endl;
    }
    cout<<endl;
}



// vector<MinCover> get_min_covers(SequenceAnnotation &seq_annot){
//     /* 
//     input: annotations (grouped by positions, sorted by asc size.)
//     variables:
//     * current rep reference for each pos: CurrRep(pos) -> rep. 
//     * first not visited rep reference for each pos: NonVisitedRep(pos) -> rep. 
//     * an ancestor stack and its topmost item (parent pos, parent rep).
//     * active pos
//     * was_last_move_left
//     operation: 
//     * left(pos),right(pos)
//     * right(pos, rep) = (right(pos), CurrRep(right(pos)))
//     * left(pos, rep) = (left(pos), CurrRep(left(pos)))
//     algo:
//     make mc_rep(seq)
//     set pos <- first pos, first rep
//     set curr <- all first reps
//     set was_last_move_left <- false

//     **** iterate below until break ******
//     (pos, rep) <- (active pos, curr(active pos))
//     (pos, rep)_is_first_visit = NonVisitedRep(active pos) == rep
//     if (pos, rep)_is_first_visit,
//         NonVisitedRep(active pos) <- 1_larger(rep)
        

//     (pos, rep)_is_first_visit = NOT was_last_move_left
//     was_last_move_left = false;

//     if NOT (pos, rep)_is_first_visit AND (pos, rep) is primary,
//         make mc_rep(rep)
//         if rep is not lowest(smallest at the pos)
//             ADD 1_smaller(rep) to (rep)->mc_reps
    
//     **** active ⊃ right? ******
//     if (pos, rep)_is_first_visit
//             AND (pos, rep) primary 
//             AND (pos, rep) is NOT rightmost
//             AND (pos, rep) includes right(pos, rep)
//         stack.push(pos, rep)
//         (parent pos, parent rep) <- (pos, rep) 
//         active pos <- right(pos).
//         go to next iteration.

//     **** parent ⊃ right? ******
//     if parent exists 
//             AND (pos, rep)_is_first_visit
//             AND (parent pos, parent rep) includes right(pos, rep),
//         active pos <- right(pos).
//         go to next iteration.

//     **** parent ⊃ not only current interval but also the larger next rep? ******
//     if parent exists
//             AND (pos, rep) is NOT highest(NOT largest at the pos) 
//             AND (parent pos, parent rep) includes (pos, 1_larger(rep))
//         CurrRep(pos) <- 1_larger(rep)
//         go to next iteration
    
//     **** parent ⊃ current interval (already determined) ******
//     if parent exists 
//         Add rep to (parent rep)->mc_reps
//         active pos <- left(pos)
//         was_last_move_left = true
//         if active pos == parent pos
//             parent pos <- stack.pop() (null if empty)
//         go to next iteration
    
//     **** at this place no parent exists  ******
//     check parent does not exist

//     **** next interval ⊃ current interval ******
//     if active rep is NOT highest(NOT largests at the pos),
//         CurrRep(pos) <- next(rep)
//         go to next iteration
    
//     **** current interval is highest ******
//     CurrRep(pos) <- null
//     Add rep to (seq)->mc_reps
//     if pos is rightmost
//         break
//     else
//         active pos <- right(pos).
//     go to next iteration
//     */
//     MinCover seq_mc;
//     seq_mc.id = seq_annot.id;
//     seq_mc.size = seq_annot.size;
//     seq_mc.is_rep = false;
//     vector<MinCover> mcs;

//     if(seq_annot.positions.size()==0){
//         mcs.push_back(seq_mc);
//         return mcs;
//     }
    

//     // position is already sorted
//     for(auto p: seq_annot.positions){
//         sort(p.reps.begin(), p.reps.end(), 
//         [](MaximalRepeatAnnotation rep1, MaximalRepeatAnnotation rep2) {return rep1.size < rep2.size; });
//     }

//     // index is the position's index, value is the reps's index
//     vector<uint64_t> curr_rep_indexes(seq_annot.positions.size());
//     vector<uint64_t> non_visited_rep_indexes(seq_annot.positions.size());
//     stack<tuple<uint64_t, MinCover>> ancestors;
//     optional<tuple<uint64_t, MinCover>> parent;
//     uint64_t active_pos_idx = 0;
//     while(true){
//         _print__get_min_cover__progress(curr_rep_indexes, ancestors, parent, active_pos_idx, seq_annot);
//         PositionAnnotation active_pos = seq_annot.positions[active_pos_idx];
//         uint64_t active_rep_idx = curr_rep_indexes[active_pos_idx];
//         MaximalRepeatAnnotation active_rep = active_pos.reps[active_rep_idx];
//         bool is_active_repr_sa_PRIMARY_in_the_rep = active_rep.first_repr_idx == active_pos.sa_idx;
//         bool is_active_position_rightmost = active_pos_idx+1 == seq_annot.positions.size();
//         bool is_active_rep_largest = active_rep_idx+1 == active_pos.reps.size();
//         bool is_active_first_visit = non_visited_rep_indexes[active_pos_idx] == active_rep_idx;
//         //visited
//         if(is_active_first_visit){
//             non_visited_rep_indexes[active_pos_idx] = active_rep_idx+1;
//         }
//         /* new mc */
//         optional<MinCover> new_mc;
//         if(is_active_first_visit && is_active_repr_sa_PRIMARY_in_the_rep){
//             //if the suffix array index is the same as the first repr suffix array index of the repeat.
//             MinCover rep_mc;
//             rep_mc.id = active_pos.rep_ids[active_rep_idx];
//             rep_mc.size = active_rep.size;
//             rep_mc.is_rep = true;
//             if(active_rep_idx>0){
//                 auto prev_rep_id = active_pos.rep_ids[active_rep_idx-1];
//                 Pos relative_pos = 0;
//                 rep_mc.mc_reps.push_back(make_tuple(relative_pos, prev_rep_id));
//             }
//             new_mc = rep_mc;
//             mcs.push_back(rep_mc);
//         }

//         /**** active ⊃ right? ******/
//         if(is_active_first_visit && !is_active_position_rightmost && is_active_repr_sa_PRIMARY_in_the_rep){
//             auto right_pos = seq_annot.positions[active_pos_idx+1];
//             auto right_rep = right_pos.reps[curr_rep_indexes[active_pos_idx+1]];
//             // cout << "active inclusion?: a pos " << active_pos.pos << ", a size " << active_rep.size << ", r pos " << right_pos.pos << ", r size " << right_rep.size   << endl;
//             bool does_active_interval_include_right_interval = right_pos.pos + right_rep.size <= active_pos.pos + active_rep.size;
//             if(does_active_interval_include_right_interval){
//                 //if the repeat becomes a parent, it should be primary.
//                 assert(new_mc.has_value());
//                 if(parent.has_value()){
//                     ancestors.push(parent.value());
//                 }
//                 parent = make_tuple(active_pos_idx, new_mc.value());
//                 active_pos_idx++;
//                 continue;
//             }
//         }

//         if(parent.has_value()) {
//             uint64_t parent_pos_idx = get<0>(parent.value());
//             MinCover parent_mc = get<1>(parent.value());
//             PositionAnnotation parent_pos = seq_annot.positions[parent_pos_idx];
//             MaximalRepeatAnnotation parent_rep = parent_pos.reps[curr_rep_indexes[parent_pos_idx]];

//             /**** parent ⊃ right? ******/
//             if(is_active_first_visit && !is_active_position_rightmost){
//                 auto right_pos = seq_annot.positions[active_pos_idx+1];
//                 auto right_rep = right_pos.reps[curr_rep_indexes[active_pos_idx+1]];
//                 bool does_parent_interval_include_right_interval = right_pos.pos+right_rep.size <= parent_pos.pos+parent_rep.size;
//                 if(does_parent_interval_include_right_interval){
//                     active_pos_idx++;
//                     continue;
//                 }
//             }
        
//             /**** parent ⊃ not only current interval but also the larger next rep? ******/
//             if(!is_active_rep_largest){
//                 // move to the next rep
//                 curr_rep_indexes[active_pos_idx] = active_rep_idx+1;
//                 auto next_rep = active_pos.reps[active_rep_idx+1];
//                 // parent includes the next larger rep
//                 bool does_parent_interval_include_higher_interval = active_pos.pos+next_rep.size < parent_pos.pos+parent_rep.size;
//                 if(does_parent_interval_include_higher_interval){
//                     continue;
//                 }
//             }
        
//             /**** parent ⊃ current interval (already determined) ******/
//             Pos relative_pos = active_pos.pos-parent_pos.pos;
//             parent_mc.mc_reps.push_back(make_tuple(relative_pos, active_pos.rep_ids[active_rep_idx]));
//             //active_pos_idx cannot be 0 because parent would not exist
//             assert(active_pos_idx>0);
//             active_pos_idx--;
//             if(active_pos_idx == parent_pos_idx){
//                 if(!ancestors.empty()) 
//                 parent = ancestors.top();
//                 else 
//                 parent = nullopt;
//             }
//             continue;
//         }

//         /**** at this place no parent exists  ******/
//         assert(!parent.has_value());
        
//         /**** next interval ⊃ current interval ******/
//         if(!is_active_rep_largest){
//             curr_rep_indexes[active_pos_idx] = active_rep_idx + 1;
//             continue;
//         }
    
//         /**** current interval is highest ******/
//         if(is_active_first_visit){
//             Pos relative_pos = active_pos.pos;
//             seq_mc.mc_reps.push_back(make_tuple(relative_pos, active_pos.rep_ids[active_rep_idx]));
//         }
        
//         if(is_active_position_rightmost){
//             break;
//         } else {
//             active_pos_idx++;
//         }
//     }
//     mcs.push_back(seq_mc);
//     return mcs;
// }

vector<MinCover> get_min_covers(SequenceAnnotation &seq_annot){
/* Summary
    set pos=first pos

    Iterate intervals 
    curr = (pos, Layer(pos))
    right = (pos+1, Layer(pos+1))
    upper = (pos, Layer(pos)+1)

    0. 1st visit(curr), primary(curr)?
    make mc
    ADD (rep)⊃(0, child rep) (∃child: layer(pos)>0)

    1. 1st visit(curr), primary(curr), curr ⊃ right? (∃right) 
    curr <- right, possible-parent <- curr (ancestor push)
    jump iter

    *** ∃possible-parent ****
    2. possible-parent ⊃ upper? (∃possible-parent, ∃upper)
    curr <- upper
    jump iter

    3. possible-parent NOT ⊃ upper? (∃possible-parent, ∃upper)
    // the UNIQUE left-extensible parent found
    ADD (parent rep)⊃(rel pos, rep)
    layer(pos)++

    4. possible-parent ⊃ right? (∃possible-parent, ∃right)
    curr <- right
    jump iter

    5. Came here: ALL mc_rep found for possible-parent.
    curr <- possible-parent, possible-parent <- (ancestor pop or NULL)
    jump iter

    *** NOT ∃possible-parent ****

    6. NOT ∃possible-parent, ∃upper?
    layer(pos)++
    jump iter

    7. NOT ∃possible-parent, NOT ∃upper, ∃right?
    ADD (seq)⊃(pos, rep)
    curr <- right
    jump iter

    8. NOT ∃possible-parent, NOT ∃upper, ∃right, NOT 1st visit(curr)?
    curr <- right
    jump iter

    9. NOT ∃possible-parent, NOT ∃upper, NOT ∃right?
    finish.
*/
    vector<MinCover> mcs;
    MinCover seq_mc;
    seq_mc.id = seq_annot.id;
    seq_mc.size = seq_annot.size;
    seq_mc.is_rep = false;

    if(seq_annot.positions.size()==0){
        mcs.push_back(seq_mc);
        return mcs;
    }
    
    // sort rep. position has been sorted already
    for(auto p: seq_annot.positions){
        sort(p.reps.begin(), p.reps.end(), 
        [](MaximalRepeatAnnotation rep1, MaximalRepeatAnnotation rep2) {return rep1.size < rep2.size; });
    }
    
    // set theoretical limit on iterations: each pair (pos, rep) does not visit the iteration as the 'curr' more than twice.
    uint64_t iter_lim = 0;
    for(auto p: seq_annot.positions){
        iter_lim += 2*p.reps.size();
    }

    // index is the position's index, value is the reps's index
    vector<uint64_t> rep_layers(seq_annot.positions.size(),0);
    vector<optional<uint64_t>> rep_layers_visited(seq_annot.positions.size(),nullopt);
    vector<optional<uint64_t>> rep_layers_completed(seq_annot.positions.size(),nullopt);
    stack<tuple<uint64_t, uint64_t>> ancestors;
    uint64_t curr_pos_idx = 0;
    uint64_t iter = 0;
    while(true){
        if(iter==iter_lim){
            cout << "Iteration limit: the loop may contain a flaw resulting in an infinite loop." << endl;
            assert(false);
        }
        iter++;
        _print__get_min_cover__progress(rep_layers, ancestors, curr_pos_idx, seq_annot);
        /* setup variables */
        PositionAnnotation curr_pos = seq_annot.positions[curr_pos_idx];
        uint64_t curr_rep_idx = rep_layers[curr_pos_idx];
        MaximalRepeatAnnotation curr_rep = curr_pos.reps[curr_rep_idx];
        optional<tuple<uint64_t, uint64_t>> parent;
        // 1st visit(curr)
        // rep_layers_n_visited[curr_pos_idx] == curr_rep_idx;
        bool is_curr_first_visit = !(rep_layers_visited[curr_pos_idx].has_value() && rep_layers_visited[curr_pos_idx].value() == curr_rep_idx);
        bool is_curr_completed = rep_layers_completed[curr_pos_idx].has_value() && rep_layers_completed[curr_pos_idx].value() == curr_rep_idx;
        // primary(curr)
        bool is_curr_primary = curr_rep.first_repr_idx == curr_pos.sa_idx;
        // ∃right
        bool right_exists = curr_pos_idx+1 < seq_annot.positions.size();
        // ∃lower
        bool lower_exists = curr_rep_idx>0;
        // ∃upper
        bool upper_exists = curr_rep_idx+1 < curr_pos.reps.size();
        // ∃possible-parent
        if(!ancestors.empty()){
            parent=ancestors.top();
        }
        
        bool possible_parent_exists = parent.has_value();
        optional<Pos> curr_relative_pos_to_parent_opt;
        // optional<MinCover> parent_mc_opt;
        optional<uint64_t> parent_pos_idx_opt;
        optional<uint64_t> parent_mc_idx_opt;
        // curr ⊃ right?
        bool curr_include_right_opt=false;
        // possible-parent ⊃ right?
        bool possible_parent_include_right_opt=false;
        // possible-parent ⊃ upper?
        bool possible_parent_include_upper_opt=false;
        // check inclusions
        Pos closing_pos_curr = curr_pos.pos + curr_rep.size;
        optional<Pos> closing_pos_right_opt;
        optional<Pos> closing_pos_upper_opt;
        optional<Pos> closing_pos_parent_opt;
        
        if(right_exists){
            assert(curr_pos_idx+1<seq_annot.positions.size());
            auto right_pos = seq_annot.positions[curr_pos_idx+1];
            assert(curr_pos_idx+1<rep_layers.size());
            auto right_rep = right_pos.reps[rep_layers[curr_pos_idx+1]];
            closing_pos_right_opt = right_pos.pos + right_rep.size;
        }
        
        if(upper_exists){
            assert(curr_rep_idx+1<curr_pos.reps.size());
            auto upper_rep = curr_pos.reps[curr_rep_idx+1];
            closing_pos_upper_opt = curr_pos.pos + upper_rep.size;
        }
        
        if(possible_parent_exists){
            parent_pos_idx_opt = get<0>(parent.value());
            parent_mc_idx_opt = get<1>(parent.value());
            // parent_mc_opt = mcs[mc_idx];
            PositionAnnotation parent_pos = seq_annot.positions[parent_pos_idx_opt.value()];
            MaximalRepeatAnnotation parent_rep = parent_pos.reps[rep_layers[parent_pos_idx_opt.value()]];
            closing_pos_parent_opt = parent_pos.pos + parent_rep.size;
            curr_relative_pos_to_parent_opt = curr_pos.pos-parent_pos.pos;
        }

        if(right_exists){
            curr_include_right_opt = closing_pos_right_opt.value() <= closing_pos_curr;
        }
        
        if(possible_parent_exists && right_exists)
        possible_parent_include_right_opt = closing_pos_right_opt.value() <= closing_pos_parent_opt.value();

        if(possible_parent_exists && upper_exists)
        possible_parent_include_upper_opt = closing_pos_upper_opt.value() <= closing_pos_parent_opt.value();
        
        /***********************/
        /* start state updates */
        /***********************/
        if(is_curr_first_visit){
            //push non visited
            if(rep_layers_visited[curr_pos_idx].has_value()){
                if(rep_layers_visited[curr_pos_idx].value()+1 < curr_pos.reps.size()){
                    rep_layers_visited[curr_pos_idx]=rep_layers_visited[curr_pos_idx].value()+1;
                }
            } else {
                rep_layers_visited[curr_pos_idx]=0;
            }
        }
        
        /* new mc */
        optional<uint64_t> new_mc_idx_opt;
        if(is_curr_first_visit && is_curr_primary){
            //if the suffix array index is the same as the first repr suffix array index of the repeat.
            MinCover rep_mc;
            rep_mc.id = curr_pos.rep_ids[curr_rep_idx];
            rep_mc.size = curr_rep.size;
            rep_mc.is_rep = true;
            if(lower_exists){
                cout << "Add lower" << endl;
                auto prev_rep_id = curr_pos.rep_ids[curr_rep_idx-1];
                Pos relative_pos = 0;
                rep_mc.mc_reps.push_back(make_tuple(relative_pos, prev_rep_id));
            }
            mcs.push_back(rep_mc);
            // new_mc_opt = &mcs[mcs.size()-1];
            new_mc_idx_opt = mcs.size()-1;
        }
        
        /* 1. 1st visit(curr), primary(curr), curr ⊃ right? (∃right) */
        if(is_curr_first_visit && is_curr_primary && curr_include_right_opt){
            //if the repeat becomes a parent, it should be primary.
            assert(new_mc_idx_opt.has_value());
            ancestors.push(make_tuple(curr_pos_idx, new_mc_idx_opt.value()));
            curr_pos_idx++;
            cout << "condition 1" << ": parent now has" << mcs[new_mc_idx_opt.value()].mc_reps.size() << endl;
            continue;
        }

        /* 2. possible-parent ⊃ upper? (∃possible-parent, ∃upper) */
        if(possible_parent_include_upper_opt){
            rep_layers[curr_pos_idx]++;
            cout << "condition 2" << endl;
            continue;
        }
        
        /* 3. possible-parent NOT ⊃ upper? (∃possible-parent) */
        if(possible_parent_exists && !possible_parent_include_upper_opt){
            mcs[parent_mc_idx_opt.value()].mc_reps.push_back(make_tuple(curr_relative_pos_to_parent_opt.value(), curr_pos.rep_ids[curr_rep_idx]));
            // parent_mc_opt.value().mc_reps.push_back(make_tuple(curr_relative_pos_to_parent_opt.value(), curr_pos.rep_ids[curr_rep_idx]));
            rep_layers_completed[curr_pos_idx]=curr_rep_idx;
            // the unique left-extensible parent found so no need to give focus 
            if(rep_layers[curr_pos_idx]+1<curr_pos.reps.size()){
                rep_layers[curr_pos_idx]++;
            }
            cout << "condition 3" << ": parent now has" << mcs[parent_mc_idx_opt.value()].mc_reps.size() << endl;
        }

        /* 4. possible-parent ⊃ right? (∃possible-parent, ∃right) */
        if(possible_parent_include_right_opt){
            curr_pos_idx++;
            cout << "condition 4" << endl;
            continue;
        }

        /* 5. possible-parent is now obsolete(∃possible-parent): ALL possible mc_rep found for possible-parent. */
        if(possible_parent_exists){
            //active_pos_idx cannot be 0 because parent would not exist
            curr_pos_idx = parent_pos_idx_opt.value();
            ancestors.pop();
            cout << "condition 5" << endl;
            continue;
        }

        //NOTE: write all conditions for clarity even if unnecessary. e.g. parent cannot exist from here
        /* 6. NOT ∃possible-parent, ∃upper? */
        if(upper_exists && !possible_parent_exists){
            rep_layers[curr_pos_idx]++;
            cout << "condition 6" << endl;
            continue;
        }

        if(!is_curr_completed && !upper_exists && !possible_parent_exists){
            seq_mc.mc_reps.push_back(make_tuple(curr_pos.pos, curr_pos.rep_ids[curr_rep_idx]));
        }

        /* 7. NOT ∃possible-parent, NOT ∃upper, ∃right?*/
        if(right_exists && !upper_exists && !possible_parent_exists){
            curr_pos_idx++;
            cout << "condition 7" << endl;
            continue;
        }
    
        /* 8. NOT ∃possible-parent, NOT ∃upper, ∃right, NOT 1st visit(curr)? */
        if(!is_curr_first_visit && !possible_parent_exists && !upper_exists && right_exists){
            curr_pos_idx++;
            cout << "condition 8" << endl;
            continue;
        }

        /* 9. NOT ∃possible-parent, NOT ∃upper, NOT ∃right? */
        if(!possible_parent_exists && !upper_exists && !right_exists){
            break;
        }
        assert(false);
    }
    mcs.push_back(seq_mc);
    return mcs;
}

#endif