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
        ext.repr_sa_workspace.left_cnts[i]=0;
        ext.repr_sa_workspace.right_cnts[i]=0;
    }
    // for each letter
    for(int c=0; c < ext.characters_cnt; c++){
        // c_node not exists
        if(!ext.c_nodes_open[c]) continue;

        //for each right a
        for(int a=0; a<ext.characters_cnt; a++){
            if(ext.c_nodes[c].firsts[a]<ext.c_nodes[c].firsts[a+1]){
                ext.repr_sa_workspace.left_cnts[c]++;
                ext.repr_sa_workspace.left_paired_example[c]=a;
                ext.repr_sa_workspace.right_cnts[a]++;
                ext.repr_sa_workspace.right_paired_example[a]=c;
            }
        }
    }
    /* find repr characters */
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        // cout << "left: " << i << ", left cnt" << left_cnt << ", right info: " << (int)r << endl;
        if(ext.repr_sa_workspace.left_cnts[i]==0){
            ext.repr_sa_workspace.left_repr[i]=false;
        } //left exclusive but not bi-exclusive
        else if(ext.repr_sa_workspace.left_cnts[i]==1 
        && ext.repr_sa_workspace.left_paired_example[i]!=0 
        && ext.repr_sa_workspace.right_cnts[ext.repr_sa_workspace.left_paired_example[i]]>1){
            ext.repr_sa_workspace.left_repr[i]=false;
        } else {
            ext.repr_sa_workspace.left_repr[i]=true;
        }
    }
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        if(ext.repr_sa_workspace.right_cnts[i]==0)
        {
            ext.repr_sa_workspace.right_repr[i]=false;
        } //right exclusive but not bi-exclusive
        else if(ext.repr_sa_workspace.right_cnts[i]==1 
        && ext.repr_sa_workspace.right_paired_example[i]!=0 
        && ext.repr_sa_workspace.left_cnts[ext.repr_sa_workspace.right_paired_example[i]]>1) {
            ext.repr_sa_workspace.right_repr[i]=false;
        } else {
            ext.repr_sa_workspace.right_repr[i]=true;
        }
    }
}

void report_repr_sa(FmIndex &index, SuffixArrayNodeExtension_NEW &ext, vector<MaximalRepeatAnnotation> &outs){
    auto start = std::chrono::steady_clock::now();
    decide_repr_sa_new(ext);
    ext.any_measure[6]+=(std::chrono::steady_clock::now()-start).count();

    start = std::chrono::steady_clock::now();
    set<SuffixArrayIdx> uniq_repr_sa;
    for(int c=0; c<ext.characters_cnt; c++){
        if(ext.repr_sa_workspace.left_repr[c]){
            //revert the rank process
            SuffixArrayIdx sa_idx=0;
            for(int i=0; i<ext.characters_cnt; i++){
                if(ext.c_nodes[c].firsts[i]<ext.c_nodes[c].firsts[i+1]) {
                    sa_idx=ext.c_nodes[c].firsts[i];
                    break;
                }
            }
            assert(sa_idx!=0);
            uint64_t rank = sa_idx - index.C[c] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
            uniq_repr_sa.insert(prev_sa_idx);
            // uniq_repr_sa.insert(ext.first_l(c));
        }
        if(ext.repr_sa_workspace.right_repr[c]){
            uniq_repr_sa.insert(ext.first_r(c));
        }
    }
    
    /*important. If termination is placed in both ends, all are representative. 
    decide_repr_sa_extensions marks it by simply including 0 in left character. 
    Then in here we should collect all suffixes where the form is #W#.
    */
    for(auto sa_idx: ext.both_ext_terms){
        uniq_repr_sa.insert(sa_idx);
    }

    vector<SuffixArrayIdx> repr_sa(uniq_repr_sa.begin(), uniq_repr_sa.end()); 

    ext.any_measure[7]+=(std::chrono::steady_clock::now()-start).count();
    // if(repr_sa.size()==0){
    //     cout << "repr_sa 0 case" << endl;
    //     for(auto pair: ext.distinct_extensions()){
    //         cout << (int)get<0>(pair) << ", " << (int)get<1>(pair) << endl;
    //     }
    //     cout << "repr ex left" << endl; 
    //     for(auto left: get<0>(decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions()))){
    //         cout << (int)left << endl;
    //     }
    //     cout << "repr ex right" << endl; 
    //     for(auto left: get<1>(decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions()))){
    //         cout << (int)left << endl;
    //     }
    //     assert(false);
    // };
    start = std::chrono::steady_clock::now();
    assert(repr_sa.size()>0);
    SuffixArrayIdx min_idx = repr_sa[0];
    for(auto i: repr_sa){
        min_idx = i < min_idx? i: min_idx;
    }
    
    // cout << "rep size: " <<ext.node.depth << " sa: "; 
    // for(auto i: repr_sa){
    //     cout << i << ", ";
    // }
    // cout << endl;

    outs.push_back({ext.node.depth, repr_sa, min_idx});
    ext.any_measure[8]+=(std::chrono::steady_clock::now()-start).count();
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

void _print__in_box_new(string txt, int box_size, char pre1=' ', char pre2=' '){
    cout<< pre1;
    cout<< pre2;
    cout<< txt;
    if(txt.size()<box_size)
        cout << std::string(box_size-txt.size(), ' ');
}

void _print__get_min_cover__progress_new(vector<uint64_t> &curr_rep_indexes, 
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
    int txt_block_size = 3;
    uint64_t max_rep=0;
    cout<< "---- algo status -----"<<endl;
    cout<< "ancestor("<<ancestors.size()<<")"<<endl;
    if(parent.has_value()){
        cout<< "parent pos: "<< seq_annot.positions[get<0>(parent.value())].pos<<endl;
    }
    if(active_rep.first_repr_idx != active_pos.sa_idx){
        cout << "curr is not primary" << endl;
    } else {
        cout << "curr is primary" << endl;
    }
    _print__in_box_new("pos", txt_block_size);
    for(auto p: seq_annot.positions){
        if(max_rep<p.reps.size())
        max_rep = p.reps.size();
        _print__in_box_new(to_string(p.pos), txt_block_size);
    }
    cout<<endl;
    for(int i=0; i<max_rep; i++){
        _print__in_box_new(" - ", txt_block_size);
        for(int p=0; p<seq_annot.positions.size(); p++){
            auto pos = seq_annot.positions[p];
            string str;
            if(i<pos.reps.size()){
                str += to_string(pos.reps[i].size);
                if(i==curr_rep_indexes[p]){
                    if(p==active_pos_idx){
                        _print__in_box_new(str, txt_block_size, '>', '*');
                    } else {
                        _print__in_box_new(str, txt_block_size, ' ', '*');
                    }
                } else {
                    _print__in_box_new(str, txt_block_size);
                }
            } else {
                _print__in_box_new(str, txt_block_size);
            }
        }
        cout<<endl;
    }
    cout<<endl;
}

void _print__seq_annot_new(SequenceAnnotation &seq_annot){   
    // sample seq annot with error.
    // vector<uint64_t> p = {0, 74 };
    // vector<vector<uint64_t>> r = {{84}, {11}, {11, 123}};
    // vector<vector<uint64_t>> f = {{84}, {11}, {11, 123}};
    // vector<PositionAnnotation> a;
    cout << "uint64_t s = " << seq_annot.size << ";" << endl;
    cout << "vector<uint64_t> p = {";
    for(int i=0; i<seq_annot.positions.size(); i++){
        cout << seq_annot.positions[i].pos;
        if(i+1< seq_annot.positions.size())
        cout << ", ";
    }
    cout<< "};" <<endl;
    cout << "vector<uint64_t> sa = {";
    for(int i=0; i<seq_annot.positions.size(); i++){
        cout << seq_annot.positions[i].sa_idx;
        if(i+1< seq_annot.positions.size())
        cout << ", ";
    }
    cout<< "};" <<endl;
    cout << "vector<vector<uint64_t>> r = {";
    for(int i=0; i<seq_annot.positions.size(); i++){
        cout << "{";
        for(int j=0; j<seq_annot.positions[i].reps.size(); j++){
            cout << seq_annot.positions[i].reps[j].size;
            if(j+1< seq_annot.positions[i].reps.size())
            cout << ", ";
        }
        cout<< "}";
        if(i+1< seq_annot.positions.size())
        cout << ", ";
    }
    cout<< "};" <<endl;
    cout << "vector<vector<uint64_t>> f = {";
    for(int i=0; i<seq_annot.positions.size(); i++){
        cout << "{";
        for(int j=0; j<seq_annot.positions[i].reps.size(); j++){
            cout << seq_annot.positions[i].reps[j].first_repr_idx;
            if(j+1< seq_annot.positions[i].reps.size())
            cout << ", ";
        }
        cout<< "}";
        if(i+1< seq_annot.positions.size())
        cout << ", ";
    }
    cout<< "};" <<endl;
}

void _print__mc_rep_addition_new(SequenceAnnotation &seq_annot, 
                            PositionAnnotation curr_pos, 
                            MaximalRepeatAnnotation curr_rep, 
                            tuple<uint64_t, uint64_t> parent,
                            string type){
    if(seq_annot.sequence.has_value()){
        auto repstr = seq_annot.sequence.value().substr(curr_pos.pos, curr_rep.size);
        if(type=="seq"){
            cout << "Added " << repstr << " to seq " << seq_annot.sequence.value() << endl;
        } else if(type=="right") {
            auto parentstr = seq_annot.sequence.value().substr(seq_annot.positions[get<0>(parent)].pos, get<1>(parent));
            cout << "Added " << repstr << " to parent rep " << parentstr << endl;
        } else if(type=="upper") {
            auto childstr = seq_annot.sequence.value().substr(seq_annot.positions[get<0>(parent)].pos, get<1>(parent));
            cout << "Added child " << childstr << " to rep " << repstr << endl;
        }
    } else {
        cout << "cannot print seq" << endl;
    }
                                
}
vector<MinCover> get_min_covers_new(SequenceAnnotation &seq_annot){
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
    
    // set theoretical limit on iterations: each pair (pos, rep) does not visit the iteration as the 'curr' more than twice.
    uint64_t iter_lim = 0;
    for(auto p: seq_annot.positions){
        iter_lim += 3*p.reps.size(); //fix later
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
            _print__get_min_cover__progress_new(rep_layers, ancestors, curr_pos_idx, seq_annot);
            _print__seq_annot_new(seq_annot);
            cout << "Iteration limit: the loop may contain a flaw resulting in an infinite loop." << endl;
            assert(false);
        }
        iter++;
        // _print__get_min_cover__progress_new(rep_layers, ancestors, curr_pos_idx, seq_annot);
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
        bool is_right_completed=false;
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
            is_right_completed = rep_layers_completed[curr_pos_idx+1].has_value() && rep_layers_completed[curr_pos_idx+1].value() == rep_layers[curr_pos_idx+1];
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
                auto prev_rep_id = curr_pos.rep_ids[curr_rep_idx-1];
                Pos relative_pos = 0;
                rep_mc.mc_reps.push_back(make_tuple(relative_pos, prev_rep_id));
                // _print__mc_rep_addition_new(seq_annot, curr_pos, curr_rep, make_tuple(curr_pos_idx, curr_pos.reps[curr_rep_idx-1].size), "upper");
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
            // cout << "condition 1" << ": parent now has" << mcs[new_mc_idx_opt.value()].mc_reps.size() << endl;
            continue;
        }

        /* 2. possible-parent ⊃ upper? (∃possible-parent, ∃upper) */
        if(possible_parent_include_upper_opt){
            rep_layers[curr_pos_idx]++;
            // cout << "condition 2 : upper " << closing_pos_upper_opt.value() << "<= parent" << closing_pos_parent_opt.value() << endl;
            continue;
        }
        
        /* 3. possible-parent NOT ⊃ upper? (∃possible-parent) */
        if(possible_parent_exists && !possible_parent_include_upper_opt){
            mcs[parent_mc_idx_opt.value()].mc_reps.push_back(make_tuple(curr_relative_pos_to_parent_opt.value(), curr_pos.rep_ids[curr_rep_idx]));
            // _print__mc_rep_addition_new(seq_annot, curr_pos, curr_rep, make_tuple(get<0>(parent.value()), mcs[get<1>(parent.value())].size), "right");
            
            rep_layers_completed[curr_pos_idx]=curr_rep_idx;
            // the unique left-extensible parent found so no need to give focus 
            if(rep_layers[curr_pos_idx]+1<curr_pos.reps.size()){
                rep_layers[curr_pos_idx]++;
            }
            // cout << "condition 3" << ": parent now has" << mcs[parent_mc_idx_opt.value()].mc_reps.size() << endl;
        }

        /* 4. possible-parent ⊃ right? (∃possible-parent, ∃right, **!is_right_not_completed) */
        if(possible_parent_include_right_opt && !is_right_completed){
            curr_pos_idx++;
            // cout << "condition 4 right is included" << endl;
            continue;
        }

        /* 5. possible-parent is now obsolete(∃possible-parent): ALL possible mc_rep found for possible-parent. */
        if(possible_parent_exists){
            //active_pos_idx cannot be 0 because parent would not exist
            curr_pos_idx = parent_pos_idx_opt.value();
            ancestors.pop();
            // cout << "condition 5 parent is now obsolete " << endl;
            continue;
        }

        //NOTE: write all conditions for clarity even if unnecessary. e.g. parent cannot exist from here
        /* 6. NOT ∃possible-parent, ∃upper? */
        if(upper_exists && !possible_parent_exists){
            rep_layers[curr_pos_idx]++;
            // cout << "condition 6 go to upper" << endl;
            continue;
        }

        if(!is_curr_completed && !upper_exists && !possible_parent_exists){
            // cout << "condition seq Add to seq" << endl;
            seq_mc.mc_reps.push_back(make_tuple(curr_pos.pos, curr_pos.rep_ids[curr_rep_idx]));
            // _print__mc_rep_addition_new(seq_annot, curr_pos, curr_rep, make_tuple(0, 0), "seq");
        }

        /* 7. NOT ∃possible-parent, NOT ∃upper, ∃right?*/
        if(right_exists && !upper_exists && !possible_parent_exists){
            curr_pos_idx++;
            // cout << "condition 7 no upper go right" << endl;
            continue;
        }
    
        /* 8. NOT ∃possible-parent, NOT ∃upper, ∃right, NOT 1st visit(curr)? */
        if(!is_curr_first_visit && !possible_parent_exists && !upper_exists && right_exists){
            curr_pos_idx++;
            // cout << "condition 8 no upper go right but no first visit" << endl;
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