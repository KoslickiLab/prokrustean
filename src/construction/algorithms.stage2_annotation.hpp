#ifndef CONSTRUCTION_ALGO_SUFFIX_ANNOT_HPP_
#define CONSTRUCTION_ALGO_SUFFIX_ANNOT_HPP_
#include <algorithm>
#include <stack>
#include <tuple>
#include "algorithms.stage1_projection.hpp"
#include "../prokrustean.support.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;


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
struct SuffixArrayAnnotationBlock{
    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
    // query
    bit_vector sa_bv;
    rank_support_v<> sa_rb;
    vector<uint8_t> sa_idx_abundances;
    vector<SuccinctStratifiedData> sa_idx_data;
    int sa_idx_queriable_left;

    void set_contents(RawStratifiedRegionBlock &raw_block){
        if(raw_block.content_count==0){
            return;
        }
        // sa indices
        this->sa_bv.resize(block_size);
        for(auto idx: raw_block.raw_sa_indices){
            sa_bv[idx]=true;
        }
        this->sa_rb=rank_support_v<>(&this->sa_bv);
        
        // content indices
        int sa_idx_cnt=this->sa_rb.rank(this->sa_bv.size());
        this->sa_idx_abundances.resize(sa_idx_cnt, 0);
        this->sa_idx_data.resize(sa_idx_cnt);

        vector<vector<int>> raw_indices_by_content_indices(sa_idx_cnt);
        for(int i=0; i< raw_block.content_count; i++){
            int content_idx = this->sa_rb.rank(raw_block.raw_sa_indices[i]);
            this->sa_idx_abundances[content_idx]++;
            raw_indices_by_content_indices[content_idx].push_back(i);
        }

        vector<StratumId> stratum_ids;
        vector<bool> is_primaries;
        for(int content_idx=0; content_idx<sa_idx_cnt; content_idx++){
            stratum_ids.clear();
            is_primaries.clear();
            for(auto i: raw_indices_by_content_indices[content_idx]){
                stratum_ids.push_back(raw_block.raw_stratum_ids[i]);
                is_primaries.push_back(raw_block.raw_is_primarys[i]);
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
        
        this->sa_idx_data[sa_idx_rank].get_data(annot_cnt, stratum_ids, is_primaries);
        // will be fetched only once
        this->sa_idx_data[sa_idx_rank].dispose();
        return true;
    }
};


struct SuffixAnnotationWorkspace{
    /* 
    * 1. collect stratums by each thread
    * 2. collect projected stratified regions
    */
    // manage block to occupy less space for stratum size
    // uint8_t stratum_block_unit=numeric_limits<uint8_t>::max();
    //
    int block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();

    int block_count;

    vector<SuffixArrayAnnotationBlock*> suffix_array_annot_blocks;

    uint64_t seq_total_length;

    SuffixAnnotationWorkspace(uint64_t seq_total_length) :seq_total_length(seq_total_length){
        this->block_count=seq_total_length/this->block_size+1;
        this->suffix_array_annot_blocks=vector<SuffixArrayAnnotationBlock*>(this->block_count);
    }
    bool fetch(SuffixArrayIdx sa_idx, int &cnt, vector<StratumId> &stratum_ids, vector<bool> &is_primaries){
        auto block_idx=sa_idx/block_size;
        if(this->suffix_array_annot_blocks[block_idx]==nullptr){
            return false;
        }
        bool fetched=this->suffix_array_annot_blocks[block_idx]->fetch(sa_idx%block_size, cnt, stratum_ids, is_primaries);
        return fetched;
    }

    void set_block(int block_idx, StratumProjectionWorkspace &stage1_workspace){
        if(stage1_workspace.raw_data_blocks[block_idx]==nullptr){
            return;
        }
        this->suffix_array_annot_blocks[block_idx]=new SuffixArrayAnnotationBlock();
        this->suffix_array_annot_blocks[block_idx]->set_contents(*stage1_workspace.raw_data_blocks[block_idx]);
        stage1_workspace.raw_data_blocks[block_idx]->dispose_block();
        delete stage1_workspace.raw_data_blocks[block_idx];
    }
};

#endif