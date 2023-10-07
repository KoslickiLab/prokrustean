#ifndef CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#include "../fm_index/tree_new.hpp"
#include "../fm_index/ssa.hpp"
#include <algorithm>
#include <stack>
#include <tuple>

using namespace std;

struct ProjectedStratifiedRegion{
    StratumId stratum_id;
    Pos pos;
    bool is_primary;

    ProjectedStratifiedRegion(StratumId stratum_id, Pos pos, bool is_primary): stratum_id(stratum_id), pos(pos), is_primary(is_primary)
    {}
};


struct StratumProjectionWorkspace{
    /* 
    * 1. collect stratums by each thread
    * 2. collect projected stratified regions
    */
    // manage block to occupy less space for stratum size
    // uint8_t stratum_block_unit=numeric_limits<uint8_t>::max();
    //
    // atomic<StratumId>* stratum_id_generator;
    StratumId new_stratum_id=0;
    //
    SpinLock stratum_lock;
    //
    Prokrustean& prokrustean;
    //
    uint64_t stratum_reserved=0;
    //
    vector<vector<ProjectedStratifiedRegion>> sequence_regions;

    int block_size;

    vector<unordered_map<uint16_t, vector<ProjectedStratifiedRegion>>> raw_region_blocks;
    // vector<vector<tuple<uint16_t, ProjectedStratifiedRegion>>> raw_region_blocks;

    // vector<unordered_map<uint16_t, bool>> raw_region_block_numbers;

    //
    vector<SpinLock> region_block_locks;
    //
    vector<SpinLock> sequence_locks;

    uint64_t seq_cnt;

    uint64_t seq_total_length;

    ProkrusteanOptional* prokrustean_optional;

    StratumProjectionWorkspace(Prokrustean &prokrustean, FmIndex &fm_index, ProkrusteanOptional* prokrustean_optional)
    :prokrustean(prokrustean),seq_cnt(fm_index.seq_cnt()),seq_total_length(fm_index.size()){
        this->sequence_regions=vector<vector<ProjectedStratifiedRegion>>(seq_cnt);
        this->sequence_locks=vector<SpinLock>(seq_cnt);
        this->block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
        this->raw_region_blocks=vector<unordered_map<uint16_t, vector<ProjectedStratifiedRegion>>>(seq_total_length/this->block_size+1);
        this->region_block_locks=vector<SpinLock>(seq_total_length/this->block_size+1);
        this->prokrustean_optional=prokrustean_optional;
        // this->raw_region_blocks=vector<vector<tuple<uint16_t, ProjectedStratifiedRegion>>>(seq_total_length/this->block_size+1);
        // this->raw_region_block_numbers=vector<unordered_map<uint16_t,bool>>(seq_total_length/this->block_size+1);
    }

    StratumId make_stratum(StratumSize size){
        this->stratum_lock.lock();
        /* critical region */
        StratumId new_id=this->new_stratum_id;
        this->new_stratum_id++;
        this->prokrustean.stratums__size.push_back(size);
        if(this->new_stratum_id==this->stratum_reserved){
            this->update_reserve_amount();
            this->prokrustean.stratums__size.reserve(this->stratum_reserved);
        }
        if(this->prokrustean_optional->collect_ext){
            this->prokrustean_optional->stratum_left_ext_count.resize(this->prokrustean.stratums__size.size(), 0);
            this->prokrustean_optional->stratum_right_ext_count.resize(this->prokrustean.stratums__size.size(), 0);
        }
        /* critical region */
        stratum_lock.unlock();
        return new_id;
    }

    void update_reserve_amount(){
        // todo: clever strategy
        this->stratum_reserved+=pow(10, 7);
    }

    // void add_stratified_regions(tuple<SeqId, Pos> loc, StratumId stratum_id, bool is_primary){
    //     this->sequence_locks[get<0>(loc)].lock();
    //     this->sequence_regions[get<0>(loc)].push_back(ProjectedStratifiedRegion(stratum_id, get<1>(loc), is_primary));
    //     this->sequence_locks[get<0>(loc)].unlock();
    // }

    void add_projected_regions(SuffixArrayIdx sa_idx, StratumId stratum_id, bool is_primary){        
        auto block_idx=sa_idx/block_size;
        auto local_sa_idx=sa_idx%block_size;
        
        // lock
        this->region_block_locks[block_idx].lock();
        // add 
        if(this->raw_region_blocks[block_idx].count(local_sa_idx)==0){
            this->raw_region_blocks[block_idx][local_sa_idx]=vector<ProjectedStratifiedRegion>();
        }
        this->raw_region_blocks[block_idx][local_sa_idx].push_back(ProjectedStratifiedRegion(stratum_id, sa_idx, is_primary));
        // this->raw_region_block_numbers[block_idx][local_sa_idx]=true;
        // unlock
        this->region_block_locks[block_idx].unlock();
    }

    void add_stratum_optional(StratumId stratum_id, uint8_t left_ext_cnt, uint8_t right_ext_cnt){
        this->prokrustean_optional->stratum_left_ext_count[stratum_id]=left_ext_cnt;
        this->prokrustean_optional->stratum_right_ext_count[stratum_id]=right_ext_cnt;
    }

    optional<vector<ProjectedStratifiedRegion>*> fetch(SuffixArrayIdx sa_idx){
        auto block_idx=sa_idx/block_size;
        auto local_sa_idx=sa_idx%block_size;
        if(this->raw_region_blocks[block_idx].count(local_sa_idx)>0){
            return &this->raw_region_blocks[block_idx][local_sa_idx];
        } else {
            return std::nullopt;
        }
        // if(this->raw_region_block_numbers[block_idx].count(local_sa_idx)>0){
        //     vector<ProjectedStratifiedRegion*> regions;
        //     for(auto &r: this->raw_region_blocks[block_idx]){
        //         if(get<0>(r)==local_sa_idx){
        //             regions.push_back(&get<1>(r));
        //         }
        //     }
        //     return regions;
        // } else {
        //     return nullopt;
        // }
    }

    void setup_prokrustean(){
        // sequence size is inferred at step2
        this->prokrustean.sequences__size.resize(this->seq_cnt);
        this->prokrustean.sequences__region.resize(this->seq_cnt);
        this->prokrustean.sequences__region_cnt.resize(this->seq_cnt, 0);
        
        // stratum size is already collected
        this->prokrustean.stratums__region.resize(this->prokrustean.stratums__size.size());
        this->prokrustean.stratums__region_cnt.resize(this->prokrustean.stratums__size.size(), 0);
    }

    void dispose(SuffixArrayIdx sa_idx){
        auto block_idx=sa_idx/block_size;
        auto local_sa_idx=sa_idx%block_size;
        this->raw_region_blocks[block_idx][local_sa_idx].clear();
        this->raw_region_blocks[block_idx][local_sa_idx].shrink_to_fit();
    }

    uint64_t get_cardinality() {
        uint64_t cnt=0;
        for (auto& m: this->raw_region_blocks){
            cnt+=m.size();
        }
        return cnt;
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
        // output.add_stratified_regions(workspace.repr_work.locations[i], workspace.stratum_id, primary_idx==i);
        output.add_projected_regions(workspace.repr_work.sa_indices[i], workspace.stratum_id, primary_idx==i);
    }

    if(output.prokrustean_optional->collect_ext){
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