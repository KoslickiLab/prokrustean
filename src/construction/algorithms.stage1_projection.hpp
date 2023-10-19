#ifndef CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#include <algorithm>
#include <stack>
#include <tuple>
#include "algorithms.stage1_tree_navigation.hpp"
#include "../prokrustean.support.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"
using namespace std;
using namespace sdsl;


struct RawStratifiedRegionBlock{
    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
    // raw
    int content_count=0;
    vector<SuffixArrayIdx_InBlock> raw_sa_indices;
    vector<StratumId> raw_stratum_ids;
    vector<bool> raw_is_primarys;

    void add_projected_region(SuffixArrayIdx_InBlock local_sa_idx, StratumId stratum_id, bool is_primary){
        this->raw_stratum_ids.push_back(stratum_id);
        this->raw_sa_indices.push_back(local_sa_idx);
        this->raw_is_primarys.push_back(is_primary);
        content_count++;
    }

    void dispose_block(){
        this->raw_stratum_ids.clear();
        this->raw_stratum_ids.shrink_to_fit();
        this->raw_sa_indices.clear();
        this->raw_sa_indices.shrink_to_fit();
        this->raw_is_primarys.clear();
        this->raw_is_primarys.shrink_to_fit();
    }
};



struct StratumProjectionWorkspace{
    /* 
    * 1. collect stratums by each thread
    * 2. collect projected stratified regions
    */
    // manage block to occupy less space for stratum size
    // uint8_t stratum_block_unit=numeric_limits<uint8_t>::max();
    StratumId new_stratum_id=0;
    //
    SpinLock stratum_lock;
    //
    Prokrustean& prokrustean;

    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();

    int block_count;

    vector<RawStratifiedRegionBlock*> raw_data_blocks;

    vector<SpinLock> block_locks;

    uint64_t seq_cnt;

    uint64_t seq_total_length;

    ProkrusteanExtension* prokrustean_optional;

    StratumProjectionWorkspace(Prokrustean &prokrustean, uint64_t seq_cnt, uint64_t seq_total_length, ProkrusteanExtension* prokrustean_optional)
    :prokrustean(prokrustean),seq_cnt(seq_cnt),seq_total_length(seq_total_length){
        this->block_count=seq_total_length/this->block_size+1;
        this->block_locks=vector<SpinLock>(this->block_count);
        this->raw_data_blocks=vector<RawStratifiedRegionBlock*>(this->block_count);
        this->prokrustean_optional=prokrustean_optional;
    }

    StratumId make_stratum(StratumSize size){
        this->stratum_lock.lock();
        /* critical region */
        StratumId new_id=this->new_stratum_id;
        this->new_stratum_id++;
        this->prokrustean.set_stratum_size(size);
        if(this->prokrustean_optional!=nullptr && this->prokrustean_optional->collect_left_right_extensions){
            this->prokrustean_optional->stratum_left_ext_count.resize(this->new_stratum_id, 0);
            this->prokrustean_optional->stratum_right_ext_count.resize(this->new_stratum_id, 0);
        }
        /* critical region */
        stratum_lock.unlock();
        return new_id;
    }


    void add_stratum_optional(StratumId stratum_id, uint8_t left_ext_cnt, uint8_t right_ext_cnt){
        this->prokrustean_optional->stratum_left_ext_count[stratum_id]=left_ext_cnt;
        this->prokrustean_optional->stratum_right_ext_count[stratum_id]=right_ext_cnt;
    }

    void add_projected_region(SuffixArrayIdx sa_idx, StratumId stratum_id, bool is_primary){
        int block_idx= sa_idx/this->block_size;
        this->block_locks[block_idx].lock();
        if(this->raw_data_blocks[block_idx]==nullptr){
            this->raw_data_blocks[block_idx]=new RawStratifiedRegionBlock();
        }
        this->raw_data_blocks[block_idx]->add_projected_region(sa_idx%block_size, stratum_id, is_primary);
        this->block_locks[block_idx].unlock();
    }

    void prepare_prokrustean_spaces(){
        this->prokrustean.set_seq_count(this->seq_cnt);
        this->prokrustean.set_stratum_count(this->new_stratum_id);
    }

};

void decide_representative(TreeWorkspace &ext){
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


void report_representative_locations(FmIndex &index, TreeWorkspace &workspace, StratumProjectionWorkspace &output){

    decide_representative(workspace);
    
    workspace.stratum_id=output.make_stratum(workspace.node.depth);
    // workspace.repr_work.locations.clear();
    workspace.repr_work.sa_indices.clear();
    for(int c=0; c<workspace.characters_cnt; c++){
        // the first suffix index of Wa form is explored by default.
        if(workspace.repr_work.right_repr[c]){
            // workspace.repr_work.locations.push_back(index.locator->get_location(workspace.first_r(c)));
            workspace.repr_work.sa_indices.push_back(workspace.first_r(c));
        }
        // the first suffix index of cW should be found and checked for duplication
        if(workspace.repr_work.left_repr[c]){
            //revert the rank process
            SuffixArrayIdx sa_idx=0;
            CharId a;
            for(int i=0; i<workspace.characters_cnt; i++){
                if(workspace.c_nodes[c].firsts[i]<workspace.c_nodes[c].firsts[i+1]) {
                    // first suffix array index of cW 
                    sa_idx=workspace.c_nodes[c].firsts[i];
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
            if(workspace.repr_work.right_repr[a] && prev_sa_idx==workspace.first_r(a)){
            } else {
                // workspace.repr_work.locations.push_back(index.locator->get_location(prev_sa_idx));
                workspace.repr_work.sa_indices.push_back(prev_sa_idx);
            }
        }
    }
    
    /*important. If termination is placed in both ends, all are representative. 
    decide_repr_sa_extensions marks it by simply including 0 in left character. 
    Then in here we should collect all suffixes where the form is #W#.
    */
    for(auto sa_idx: workspace.both_ext_terms){
        // workspace.repr_work.locations.push_back(index.locator->get_location(sa_idx));
        workspace.repr_work.sa_indices.push_back(sa_idx);
    }

    // find primary
    // auto loc_cnt=workspace.repr_work.locations.size();
    auto idx_cnt=workspace.repr_work.sa_indices.size();
    auto primary_idx = 0;
    for(int i=1; i<idx_cnt; i++){
        // if(workspace.repr_work.locations[primary_idx] < workspace.repr_work.locations[i]){
        //     primary_idx=i;
        // }
        if(workspace.repr_work.sa_indices[i]<workspace.repr_work.sa_indices[primary_idx]){
            primary_idx=i;
        }
    }
    for(int i=0; i<idx_cnt; i++){
        // output.add_projected_regions(workspace.repr_work.sa_indices[i], workspace.stratum_id, primary_idx==i);
        // output.raw_data_by_threads[workspace.thread_idx].add_projected_region(workspace.repr_work.sa_indices[i], workspace.stratum_id, primary_idx==i);
        // output.add_projected_regions_single(workspace.thread_idx, workspace.repr_work.sa_indices[i], workspace.stratum_id, primary_idx==i);
        output.add_projected_region(workspace.repr_work.sa_indices[i], workspace.stratum_id, primary_idx==i);
    }

    if(output.prokrustean_optional!=nullptr && output.prokrustean_optional->collect_left_right_extensions){
        uint8_t left_ext_cnt=0; 
        uint8_t right_ext_cnt=0;
        // for each letter
        for(int c=0; c < workspace.characters_cnt; c++){
            if(c==0) continue;

            if(workspace.node.firsts[c]<workspace.node.firsts[c+1]){
                right_ext_cnt++;
            }
            if(workspace.c_nodes_open[c]) {
                left_ext_cnt++;;
            }
        }
        output.add_stratum_optional(workspace.stratum_id, left_ext_cnt, right_ext_cnt);
    }
}


#endif