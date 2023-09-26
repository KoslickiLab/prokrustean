#ifndef CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#include "models.hpp"
#include "../fm_index/tree_new.hpp"
#include "../fm_index/ssa.hpp"
#include <algorithm>
#include <stack>
#include <tuple>

using namespace std;

void decide_repr_sa_new(TreeWorkspace &ext){
    // initialize cnts
    for(int i=0; i<ext.characters_cnt; i++){
        ext.repr_work.left_cnts[i]=0;
        ext.repr_work.right_cnts[i]=0;
    }
    // for each letter
    for(int c=0; c < ext.characters_cnt; c++){
        // c_node not exists
        if(!ext.c_nodes_open[c]) continue;

        //for each right a
        for(int a=0; a<ext.characters_cnt; a++){
            if(ext.c_nodes[c].firsts[a]<ext.c_nodes[c].firsts[a+1]){
                ext.repr_work.left_cnts[c]++;
                ext.repr_work.left_paired_a_char[c]=a;
                ext.repr_work.right_cnts[a]++;
                ext.repr_work.right_paired_c_char[a]=c;
            }
        }
    }
    /* find repr characters */
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        // cout << "left: " << i << ", left cnt" << left_cnt << ", right info: " << (int)r << endl;
        if(ext.repr_work.left_cnts[i]==0){
            ext.repr_work.left_repr[i]=false;
        } //left exclusive but not bi-exclusive
        else if(ext.repr_work.left_cnts[i]==1 
        && ext.repr_work.left_paired_a_char[i]!=0 
        && ext.repr_work.right_cnts[ext.repr_work.left_paired_a_char[i]]>1){
            ext.repr_work.left_repr[i]=false;
        } else {
            ext.repr_work.left_repr[i]=true;
        }
    }
    for(int i=1/*not termination*/; i<ext.characters_cnt; i++){ 
        if(ext.repr_work.right_cnts[i]==0)
        {
            ext.repr_work.right_repr[i]=false;
        } //right exclusive but not bi-exclusive
        else if(ext.repr_work.right_cnts[i]==1 
        && ext.repr_work.right_paired_c_char[i]!=0 
        && ext.repr_work.left_cnts[ext.repr_work.right_paired_c_char[i]]>1) {
            ext.repr_work.right_repr[i]=false;
        } else {
            ext.repr_work.right_repr[i]=true;
        }
    }
}


void report_repr_sa(FmIndex &index, TreeWorkspace &work, StratificationOutput &out){

    decide_repr_sa_new(work);
    
    // work.stratum_id=out.make_stratum(work.node.depth);
    work.repr_work.locations.clear();
    for(int c=0; c<work.characters_cnt; c++){
        // the first suffix index of Wa form is explored by default.
        if(work.repr_work.right_repr[c]){
            work.repr_work.locations.push_back(index.locator->get_location(work.first_r(c)));
        }
        // the first suffix index of cW should be found and checked for duplication
        if(work.repr_work.left_repr[c]){
            //revert the rank process
            SuffixArrayIdx sa_idx=0;
            CharId a;
            for(int i=0; i<work.characters_cnt; i++){
                if(work.c_nodes[c].firsts[i]<work.c_nodes[c].firsts[i+1]) {
                    // first suffix array index of cW 
                    sa_idx=work.c_nodes[c].firsts[i];
                    a=i;
                    break;
                }
            }
            assert(sa_idx!=0);
            // auto start2 = std::chrono::steady_clock::now();
            uint64_t rank = sa_idx - index.C[c] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
            // ext.any_measure[6]+=(std::chrono::steady_clock::now()-start2).count();
            
            // Check duplication, i.e. if the suffix array picked by the right letter a of cWa.
            if(work.repr_work.right_repr[a] && prev_sa_idx==work.first_r(a)){
            } else {
                work.repr_work.locations.push_back(index.locator->get_location(prev_sa_idx));
            }
        }
    }
    
    /*important. If termination is placed in both ends, all are representative. 
    decide_repr_sa_extensions marks it by simply including 0 in left character. 
    Then in here we should collect all suffixes where the form is #W#.
    */
    for(auto sa_idx: work.both_ext_terms){
        work.repr_work.locations.push_back(index.locator->get_location(sa_idx));
    }

    // find primary
    auto loc_cnt=work.repr_work.locations.size();
    auto primary_idx = 0;
    for(int i=1; i<loc_cnt; i++){
        if(work.repr_work.locations[primary_idx] < work.repr_work.locations[i]){
            primary_idx=i;
        }
    }
    // for(int i=1; i<loc_cnt; i++){
    //     out.add_stratified_regions(work.repr_work.locations[i], work.stratum_id, primary_idx==i);
    // }
}


#endif