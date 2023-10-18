#ifndef CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#define CONSTRUCTION_ALGO_PROCEDURES_NEW_HPP_
#include <algorithm>
#include <stack>
#include <tuple>
#include "algorithms.stage1_tree.hpp"
#include "../prokrustean.enhance.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"
using namespace std;
using namespace sdsl;

struct ProjectedStratifiedRegion{
    StratumId stratum_id;
    Pos pos;
    bool is_primary;

    ProjectedStratifiedRegion(StratumId stratum_id, Pos pos, bool is_primary): stratum_id(stratum_id), pos(pos), is_primary(is_primary)
    {}
};

struct SuccinctStratifiedData {
    StratumId* data;

    void set_data(uint8_t annot_cnt, vector<StratumId> &stratum_ids, vector<bool> &is_primaries){
        assert(annot_cnt>0);
        uint8_t bits_per_type=(8*sizeof(StratumId));
        uint8_t extra=annot_cnt%bits_per_type>0? annot_cnt/bits_per_type+1 : annot_cnt/bits_per_type;
        data=new StratumId[annot_cnt+extra]; // 8* size means bits available so that is_primaries are added
        for(int i=0; i<annot_cnt; i++){
            data[i]=stratum_ids[i];
        }
        StratumId result=0;
        if(extra==1){
            for (uint8_t i = 0; i < annot_cnt; ++i) {
                result |= (is_primaries[i] ? 1 : 0) << (i % bits_per_type);
            }
            data[annot_cnt]=result;
        } else {
            vector<StratumId> results;
            for (uint8_t i = 0; i < annot_cnt; ++i) {
                result |= (is_primaries[i] ? 1 : 0) << (i % bits_per_type);
                if ((i + 1) % bits_per_type == 0) {
                    results.push_back(result);
                    result = 0;
                }
            }
            results.push_back(result);
            assert(results.size()==extra);
            for(int i=0; i<extra; i++){
                data[annot_cnt+i]=results[i];
            }
        }
    }
    void get_data(int annot_cnt, vector<StratumId> &stratum_ids, vector<bool> &is_primaries){
        stratum_ids.resize(annot_cnt);
        is_primaries.resize(annot_cnt);
        int bits_per_type=(8*sizeof(StratumId));
        uint8_t extra=annot_cnt%bits_per_type>0? annot_cnt/bits_per_type+1 : annot_cnt/bits_per_type;
        for (int8_t i = 0; i < annot_cnt; ++i) {
            stratum_ids[i]=data[i];
        }
        if(extra==1){
            for (int8_t i = 0; i < annot_cnt; ++i) {
                is_primaries[i] = ((data[annot_cnt] >> (i % bits_per_type)) & 1) != 0;
            }
        } else {
            for(int idx=0; idx<extra; idx++){
                for (int8_t i = 0; i < bits_per_type; ++i) {
                    is_primaries[bits_per_type*idx+i] = ((data[annot_cnt+idx] >> (i % bits_per_type)) & 1) != 0;
                }
            }
        }
    }
    void dispose(){
        delete data;
    }
};

struct RawStratifiedRegionBlock{
    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
    // raw
    vector<SuffixArrayIdx_InBlock> raw_sa_indices;
    vector<StratumId> raw_stratum_ids;
    vector<bool> raw_is_primarys;

    void add_projected_region(SuffixArrayIdx_InBlock local_sa_idx, StratumId stratum_id, bool is_primary){
        this->raw_stratum_ids.push_back(stratum_id);
        this->raw_sa_indices.push_back(local_sa_idx);
        this->raw_is_primarys.push_back(is_primary);
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


struct SuffixArrayAnnotationBlock{
    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
    // query
    bit_vector sa_bv;
    rank_support_v<> sa_rb;
    vector<uint8_t> sa_idx_abundances;
    vector<SuccinctStratifiedData> sa_idx_data;
    int sa_idx_queriable_left;

    void set_contents(RawStratifiedRegionBlock &raw_block){
        if(raw_block.raw_sa_indices.size()==0){
            return;
        }
        auto &raw_sa_indices = raw_block.raw_sa_indices;
        auto &raw_is_primarys = raw_block.raw_is_primarys;
        auto &raw_stratum_ids = raw_block.raw_stratum_ids;
        // sa indices
        // unordered_map<SuffixArrayIdx_InBlock, int> counts_of_idx;
        this->sa_bv.resize(block_size);
        for(auto idx: raw_sa_indices){
            sa_bv[idx]=true;
        }
        this->sa_rb=rank_support_v<>(&this->sa_bv);
        
        // content indices
        int sa_idx_cnt=this->sa_rb.rank(this->sa_bv.size());
        this->sa_idx_abundances.resize(sa_idx_cnt, 0);
        this->sa_idx_data.resize(sa_idx_cnt);

        vector<vector<int>> raw_indices_by_content_indices(sa_idx_cnt);
        vector<StratumId> stratum_ids;
        vector<bool> is_primaries;
        for(int i=0; i< raw_sa_indices.size(); i++){
            int content_idx = this->sa_rb.rank(raw_sa_indices[i]);
            this->sa_idx_abundances[content_idx]++;
            raw_indices_by_content_indices[content_idx].push_back(i);
        }
        for(int content_idx=0; content_idx<sa_idx_cnt; content_idx++){
            stratum_ids.clear();
            is_primaries.clear();
            for(auto i: raw_indices_by_content_indices[content_idx]){
                stratum_ids.push_back(raw_stratum_ids[i]);
                is_primaries.push_back(raw_is_primarys[i]);
            }
            this->sa_idx_data[content_idx].set_data(raw_indices_by_content_indices[content_idx].size(), stratum_ids, is_primaries);
        }

        sa_idx_queriable_left=sa_idx_cnt;
    }

    bool fetch(SuffixArrayIdx_InBlock local_sa_idx, int &annot_cnt, vector<StratumId> &stratum_ids, vector<bool> &is_primaries){
        if(!this->sa_bv[local_sa_idx]){
            return false;
        } 
        int sa_idx_rank=this->sa_rb.rank(local_sa_idx);
        annot_cnt=sa_idx_abundances[sa_idx_rank];
        assert(sa_idx_rank< this->sa_idx_data.size());
        this->sa_idx_data[sa_idx_rank].get_data(annot_cnt, stratum_ids, is_primaries);
        
        // will be fetched only once
        delete this->sa_idx_data[sa_idx_rank].data;
        return true;
    }
};


// struct RawDataByThread{
//     vector<SpinLock> locks;
//     int block_size;
//     // vector<vector<ProjectedStratifiedRegion>> raw_data_blocks;
//     vector<vector<StratumId>> stratum_id_blocks;
//     vector<vector<SuffixArrayIdx_InBlock>> sa_idx_blocks;
//     vector<vector<bool>> is_primary_blocks;

//     RawDataByThread(int block_size, int block_count): block_size(block_size){
//         this->locks=vector<SpinLock>(block_count);
//         // this->raw_data_blocks.resize(block_count);
//         this->stratum_id_blocks.resize(block_count);
//         this->sa_idx_blocks.resize(block_count);
//         this->is_primary_blocks.resize(block_count);
//     }

//     void add_projected_region(SuffixArrayIdx sa_idx, StratumId stratum_id, bool is_primary){
//         int block_idx=sa_idx/block_size;
//         SuffixArrayIdx_InBlock local_sa_idx=sa_idx%block_size;
//         // this->raw_data_blocks[block_idx].push_back(ProjectedStratifiedRegion(stratum_id, local_sa_idx, is_primary));
//         this->locks[block_idx].lock();
//         this->stratum_id_blocks[block_idx].push_back(stratum_id);
//         this->sa_idx_blocks[block_idx].push_back(local_sa_idx);
//         this->is_primary_blocks[block_idx].push_back(is_primary);
//         this->locks[block_idx].unlock();
//     }

//     void dispose_block(int block_idx){
//         this->stratum_id_blocks[block_idx].clear();
//         this->stratum_id_blocks[block_idx].shrink_to_fit();
//         this->sa_idx_blocks[block_idx].clear();
//         this->sa_idx_blocks[block_idx].shrink_to_fit();
//         this->is_primary_blocks[block_idx].clear();
//         this->is_primary_blocks[block_idx].shrink_to_fit();
//     }
// };

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

    vector<SuffixArrayAnnotationBlock*> suffix_array_annot_blocks;

    vector<SpinLock> block_locks;

    uint64_t seq_cnt;

    uint64_t seq_total_length;

    ProkrusteanEnhancement* prokrustean_optional;

    StratumProjectionWorkspace(Prokrustean &prokrustean, FmIndex &fm_index, ProkrusteanEnhancement* prokrustean_optional, int thread_cnt=1)
    :prokrustean(prokrustean),seq_cnt(fm_index.seq_cnt()),seq_total_length(fm_index.size()){
        this->block_count=seq_total_length/this->block_size+1;
        this->block_locks=vector<SpinLock>(this->block_count);
        this->raw_data_blocks=vector<RawStratifiedRegionBlock*>(this->block_count);
        // for(int i=0; i<this->block_count; i++){
        //     this->raw_data_blocks[i]=new RawStratifiedRegionBlock();
        // }
        this->suffix_array_annot_blocks=vector<SuffixArrayAnnotationBlock*>(this->block_count);
        // for(auto &block: this->suffix_array_annot_blocks){
        //     block.block_size=this->block_size;
        // }
        this->prokrustean_optional=prokrustean_optional;
    }

    StratumId make_stratum(StratumSize size){
        this->stratum_lock.lock();
        /* critical region */
        StratumId new_id=this->new_stratum_id;
        this->new_stratum_id++;
        this->prokrustean.stratums__size.push_back(size);
        if(this->prokrustean_optional!=nullptr && this->prokrustean_optional->collect_left_right_extensions){
            this->prokrustean_optional->stratum_left_ext_count.resize(this->prokrustean.stratums__size.size(), 0);
            this->prokrustean_optional->stratum_right_ext_count.resize(this->prokrustean.stratums__size.size(), 0);
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


    bool fetch(SuffixArrayIdx sa_idx, int &cnt, vector<StratumId> &stratum_ids, vector<bool> &is_primaries){
        auto block_idx=sa_idx/block_size;
        if(this->suffix_array_annot_blocks[block_idx]==nullptr){
            return false;
        }
        bool fetched=this->suffix_array_annot_blocks[block_idx]->fetch(sa_idx%block_size, cnt, stratum_ids, is_primaries);
        return fetched;
    }

    void set_block(int block_idx){
        if(this->raw_data_blocks[block_idx]==nullptr){
            return;
        }
        
        this->suffix_array_annot_blocks[block_idx]=new SuffixArrayAnnotationBlock();
        this->suffix_array_annot_blocks[block_idx]->set_contents(*this->raw_data_blocks[block_idx]);
        this->raw_data_blocks[block_idx]->dispose_block();
        delete this->raw_data_blocks[block_idx];
    }

    void prepare_prokrustean_spaces(){
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
    }

    uint64_t get_cardinality() {
        uint64_t cnt=0;
        // for (auto& m: this->raw_region_blocks){
        //     cnt+=m.size();
        // }
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