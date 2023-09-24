#ifndef CONSTRUCTION_MODEL_HPP_
#define CONSTRUCTION_MODEL_HPP_

#include <vector>
#include <list>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

typedef uint32_t StratumIdxUnionIsPrimary;
typedef uint16_t SuffixArrayIdx_InBlock; //2 byte size block
typedef uint8_t StratumIdx_PerReprSA; //how many stratum interval starts at each position. 255 has been enough
const SuffixArrayIdx_InBlock CONSTRUCTION_BLOCK_UNIT = numeric_limits<SuffixArrayIdx_InBlock>::max();  
// to save space at storing stratum ids at repr locations, the isPrimary flag uses one bit.
// need to modify to uint64_t StratumIdxUnionIsPrimary if not enough.
const StratumIdxUnionIsPrimary STRATUM_ID_LIMIT=numeric_limits<StratumIdxUnionIsPrimary>::max()/2;

struct MaximalRepeatAnnotation {
    // 
    uint64_t size;

    vector<SuffixArrayIdx> repr_indexes;

    SuffixArrayIdx first_repr_idx;
};

struct PositionAnnotation {
    Pos pos; 
    
    SuffixArrayIdx sa_idx; 

    // order desc for rep sizes
    vector<MaximalRepeatAnnotation> reps;

    // corresponding ids
    vector<StratumId> rep_ids;
};

struct SequenceAnnotation {
    //
    SeqId id;

    // 
    uint64_t size;

    //
    vector<PositionAnnotation> positions;

    //
    optional<string> sequence; 
};

struct ReprSuffixRank {
    //
    sdsl::bit_vector bv;
    //
    sdsl::rank_support_v<> rb;

    ReprSuffixRank(){}

    void initialize(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        bv.resize(seq_length);
        for(auto rep: repeats){
            for(auto sa_idx: rep.repr_indexes){
                bv[sa_idx]=true;
            }
        }
        
        rank_support_v<> _rb(&bv);
        this->rb = _rb;
    }
    
    uint64_t get_repr_size(){
        return rb.rank(bv.size());
    }

    // 
    bool exists(uint64_t sa_idx){
        return bv[sa_idx];
    }

    uint64_t rank(uint64_t sa_idx){
        return rb.rank(sa_idx);
    }
};

struct ReprSuffixAnnotation {
private:
    ReprSuffixRank sa_rank;
    vector<vector<StratumId>> repr_suffixes;

public:
    void initialize_rank(uint64_t seq_length, vector<MaximalRepeatAnnotation> &repeats){
        sa_rank.initialize(seq_length, repeats);
    }

    void initialize_repr_sa(vector<MaximalRepeatAnnotation> &repeats){
        cout << "--- in suffix annotation --- " << endl;
        repr_suffixes.resize(sa_rank.get_repr_size());
        for(uint64_t id=0; id<repeats.size(); id++){
            for(auto sa_idx: repeats[id].repr_indexes){
                uint64_t r = sa_rank.rank(sa_idx);
                repr_suffixes[r].push_back(id);
            }
        }
    }

    bool exists(uint64_t sa_idx){
        return sa_rank.exists(sa_idx);
    }

    optional<vector<StratumId>> get_repeats(uint64_t sa_idx){
        if(!sa_rank.exists(sa_idx)){
            return nullopt;
        } else {
            auto rank = sa_rank.rank(sa_idx);
            return repr_suffixes[rank];
        }
    }
};



// class ReprSuffixAnnotationParallel {
//     /* thread safe imple*/
//     vector<vector<StratumId>> repr_suffixes;

//     ReprSuffixAnnotationParallel(uint64_t repr_size){
//         repr_suffixes.resize(repr_size);
//     }

//     // 
//     void set_repeats(vector<MaximalRepeatAnnotation> &repeats, uint64_t from, uint64_t to){

//     }

//     // 
//     vector<StratumId> get_repeats(uint64_t sa_idx){

//     }
// };

void print_repeats(vector<MaximalRepeatAnnotation> repeats){
    for(int i=0; i<repeats.size(); i++){
        cout << "R" << i << "("<< repeats[i].size << ")" <<": ";
        for(auto s: repeats[i].repr_indexes){
            cout << s << ", ";
        }
        cout << endl;
    }
}

void print_positions(vector<PositionAnnotation> positions){
    cout << "---- print_positions ---- " << endl;
    for(auto pos: positions){
        cout << "sa: " << pos.sa_idx<< ", ";
        for(int i=0; i< pos.reps.size(); i++){
            cout << "R" << pos.rep_ids[i] << "("<< pos.reps[i].size << ")" << ", ";
        }
        cout << endl;
    }
}



/* ************************************************************ */
/* At step1, threads use a list of blocks to fill in raw data   */
/* Thread each has their own set of blocks                       */
/* ************************************************************ */
struct StratifiedRaw{
    SuffixArrayIdx_InBlock repr_sa_in_block;
    bool is_primary_of_stratum;
    StratumId stratum_id;

    StratifiedRaw(SuffixArrayIdx_InBlock repr_id, StratumId stratum_id, bool is_primary_of_stratum)
    : repr_sa_in_block(repr_id), stratum_id(stratum_id), is_primary_of_stratum(is_primary_of_stratum) {}
};

struct StratifiedRawBlock{
    uint32_t cnt=0;
     // make private for better memory usage
    vector<StratifiedRaw> raw_regions;

    
    void add(SuffixArrayIdx repr_idx, StratumId stratum_id, bool is_primary_of_stratum){
        raw_regions.push_back(StratifiedRaw(repr_idx, stratum_id, is_primary_of_stratum));
        // for(int i=0; i<3; i++){
        //     auto *ptr2 = (StratifiedRaw*)malloc(sizeof(StratifiedRaw));
        //     temps.push_back(ptr2);
        // }
        cnt++;
    }

    void clear(){
        raw_regions.clear();
        raw_regions.shrink_to_fit();
        raw_regions=vector<StratifiedRaw>();
        // for (auto* ptr : temps) {
        //     free(ptr);
        // }
        // stack<StratifiedRaw, vector<StratifiedRaw>>().swap(temps);
    }
    // private:
    //     vector<StratifiedRaw*> temps;
};

struct StratifiedRawWorkspace{
    //
    uint8_t workspace_id;
    //
    vector<StratifiedRawBlock> raw_blocks;

    vector<tuple<StratumId, StratumSize>> stratum_raws;
    
    uint32_t stratum_cnt=0;

    uint64_t stratified_rgn_cnt=0;

    atomic<StratumId>* id_generator;

    StratifiedRawWorkspace(uint8_t workspace_id, uint64_t seq_length, atomic<StratumId> &id_generator){
        this->workspace_id=workspace_id;
        this->raw_blocks=vector<StratifiedRawBlock>(seq_length/CONSTRUCTION_BLOCK_UNIT+1);
        this->id_generator=&id_generator;
    }

    //at step1
    StratumId get_new_stratum(uint32_t size){
        StratumId id = id_generator->fetch_add(1);
        stratum_cnt++;
        stratum_raws.push_back(make_tuple(id, size));
        return id;
    }

    // repr_idx is splitted so that it is inferred.
    void add_repr_raw(SuffixArrayIdx repr_idx, StratumId stratum_id, bool is_primary_of_stratum){
        raw_blocks[repr_idx/CONSTRUCTION_BLOCK_UNIT].add(repr_idx%CONSTRUCTION_BLOCK_UNIT, stratum_id, is_primary_of_stratum);
        stratified_rgn_cnt++;
    }

    //at the middle of step1&2 - distribute stratum ids between workspaces
    void set_prokrutean_stratum(Prokrustean &prokrustean){
        for(auto &pair: stratum_raws){
            prokrustean.stratums[get<0>(pair)].size=get<1>(pair);
        }
        stratum_raws.clear();
        stratum_raws.shrink_to_fit();
        stratum_raws=vector<tuple<StratumId, StratumSize>>();
    }

    void clear_rgn_block(int block_no){
        // raw_blocks[block_no].raw_regions.clear();
        // raw_blocks[block_no].raw_regions.shrink_to_fit();
        // vector<StratifiedRaw>().swap(raw_blocks[block_no].raw_regions);
        raw_blocks[block_no].clear();
    }

    // debugging purpose
    vector<uint64_t> restore_reprs(){
        vector<uint64_t> reprs;
        for(int b=0;b<this->raw_blocks.size(); b++){
            for(auto &raw: this->raw_blocks[b].raw_regions){
                auto original = raw.repr_sa_in_block+b*CONSTRUCTION_BLOCK_UNIT;
                reprs.push_back(original);
            }
        }
        return reprs;
    }
};

/* ************************************************************ */
/* At step2, blocks assgined to threads are merged              */
/* threads are distributed to each block                          */
/* ************************************************************ */
struct Stratified_OfReprPos{
    /* repr_sa is inferred from the Block location and bit vector (repr_exists) */
    uint32_t size; 
    // fixed length whether it is primary(first of repr suffix index) of rep (dynamically allocated) to secure space efficiency
    tuple<StratumId, bool> stratum_id_array[];
};

struct StratifiedBlock{
    //
    uint32_t block_no;
    //
    bit_vector repr_exists;
    //
    rank_support_v<> repr_exists_rank;
    // index is inferred from outside. (repr_sa/block_size, repr_sa % block_size)
    vector<Stratified_OfReprPos*> reprs;
    //
    uint32_t stat__repr_cnt=0;
    //
    uint64_t stat__rgn_cnt=0;
    //
    // vector<StratumIdx_PerReprSA> stra_cnt_per_repr;

    StratifiedBlock(uint32_t block_no){
        this->block_no=block_no;
        this->repr_exists=bit_vector(CONSTRUCTION_BLOCK_UNIT, 0);
    }

    optional<Stratified_OfReprPos*> get_repr(uint16_t repr_sa_in_block){
        if(!repr_exists[repr_sa_in_block]) return nullopt;
        
        return reprs[repr_exists_rank.rank(repr_sa_in_block)];
    }

    void dispose_repr(uint16_t repr_sa_in_block){
            auto i = this->repr_exists_rank.rank(repr_sa_in_block);
            delete this->reprs[i];
            this->reprs[i]=nullptr;
    }

    void index_reprs(vector<StratifiedRawWorkspace> &workspaces){
        // set bit vector
        for(auto &work: workspaces){
            for(auto raw :work.raw_blocks[block_no].raw_regions){
                stat__rgn_cnt++;
                if(!repr_exists[raw.repr_sa_in_block]){
                    repr_exists[raw.repr_sa_in_block]=true;
                }
            }
        }

        // get stats
        this->repr_exists_rank=rank_support_v<>(&this->repr_exists);
        this->stat__repr_cnt = this->repr_exists_rank.rank(this->repr_exists.size());
        // this->stra_cnt_per_repr=vector<StratumIdx_PerReprSA>(repr_cnt, 0);
        vector<StratumIdx_PerReprSA>stat__stra_cnt_per_repr(this->stat__repr_cnt);
        for(auto &work: workspaces){
            for(auto raw: work.raw_blocks[block_no].raw_regions){
                stat__stra_cnt_per_repr[this->repr_exists_rank.rank(raw.repr_sa_in_block)]++;
            }
        }
        // prepare spaces
        this->reprs=vector<Stratified_OfReprPos*>(this->stat__repr_cnt);
        for(int i=0; i<this->stat__repr_cnt; i++){
            auto stra_cnt=stat__stra_cnt_per_repr[i];
            Stratified_OfReprPos* repr_stra = static_cast<Stratified_OfReprPos*>(std::malloc(sizeof(Stratified_OfReprPos)+stra_cnt*sizeof(tuple<StratumId, bool>)));
            repr_stra->size=stra_cnt;
            this->reprs[i]=repr_stra;
        }

        // set stratums. re-utilize the stat information
        for(auto &work: workspaces){
            for(auto raw: work.raw_blocks[block_no].raw_regions){
                auto i = this->repr_exists_rank.rank(raw.repr_sa_in_block);
                auto stratum_idx_in_repr = this->reprs[i]->size - stat__stra_cnt_per_repr[i];
                this->reprs[i]->stratum_id_array[stratum_idx_in_repr]=make_tuple(raw.stratum_id, raw.is_primary_of_stratum);

                stat__stra_cnt_per_repr[i]--;
            }
        }
    }

    // debugging purpose
    void validate(){
        uint64_t total_cnt=0;
        for(int i=0; i<stat__repr_cnt; i++){
            total_cnt+=reprs[i]->size;
        }
        assert(total_cnt==stat__rgn_cnt);
    }
};


struct StratifiedSA_ParallelModel {
    // each thread has own blocks of raw data that will be converted to real block
    vector<StratifiedRawWorkspace> parallel_workspaces;
    // parallel converged - boolean but to make it threadsafe, use 1 byte symbol
    vector<uint8_t> block_converged;
    // blocks that are each merged from correspondings in parallel workspaces
    vector<StratifiedBlock> blocks;
    //
    int parallel_scale;
    // 
    atomic<uint32_t> stratum_id_generator;

    StratifiedSA_ParallelModel(int thread_cnt, uint64_t seq_length){
        parallel_scale=thread_cnt;
        for(int i=0; i<parallel_scale; i++){
            parallel_workspaces.push_back(StratifiedRawWorkspace(i, seq_length, stratum_id_generator));    
        }
        auto block_cnt = seq_length/CONSTRUCTION_BLOCK_UNIT+1;
        blocks.reserve(block_cnt);
        for(int i=0; i<block_cnt; i++){
            blocks.push_back(StratifiedBlock(i));
        }
        block_converged=vector<uint8_t>(block_cnt, 0);
    }

    // uint64_t is too expensive. block-based id leads to uint16_t repr idx
    optional<Stratified_OfReprPos*> query(uint64_t repr_sa){
        return blocks[repr_sa/CONSTRUCTION_BLOCK_UNIT].get_repr(repr_sa%CONSTRUCTION_BLOCK_UNIT);
    }

    void prepare_prokrustean(Prokrustean &prokrustean){
        uint64_t total_stratums=0;
        uint32_t stratum_max_size;
        for(int i=0; i< parallel_workspaces.size(); i++){
            total_stratums+=parallel_workspaces[i].stratum_cnt;
        }
        prokrustean.stratums = vector<Stratum>(total_stratums);
    }

    void converge_block(uint64_t block_no){
        assert(block_converged[block_no]==0);
        block_converged[block_no]=1;
        blocks[block_no].index_reprs(parallel_workspaces);
        for(auto &work: parallel_workspaces){
            work.clear_rgn_block(block_no);
        }
    }
};

#endif