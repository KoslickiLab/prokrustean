#ifndef CONSTRUCTION_MODEL_HPP_
#define CONSTRUCTION_MODEL_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

typedef uint32_t StratumId_ThreadSpecific;
typedef uint16_t ReprIdx_BlockSpecific;

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
    ReprIdx_BlockSpecific repr_sa_in_block;
    bool is_primary_of_rep;
    // starting from thread specific -> later assigned global
    uint32_t stratum_id;

    StratifiedRaw(ReprIdx_BlockSpecific repr_id, StratumId_ThreadSpecific stratum_id, bool is_primary_of_rep)
    : repr_sa_in_block(repr_id), stratum_id(stratum_id), is_primary_of_rep(is_primary_of_rep) {}
};

struct StratifiedRawWorkspace{
    //
    uint8_t workspace_id;
    //
    int block_unit;
    //
    vector<vector<StratifiedRaw>> blocks_of_raws;
    // stratum id -> index, only store size 
    vector<StratumSize> stratum_sizes;
    // each size -> index, store the aubndance of size
    vector<uint32_t> stratum_size_stats;
    // 
    uint32_t stratum_cnt=0;

    StratifiedRawWorkspace(uint8_t workspace_id, uint64_t seq_length, int block_unit){
        this->workspace_id=workspace_id;
        this->block_unit=block_unit;
        this->blocks_of_raws=vector<vector<StratifiedRaw>>(seq_length/block_unit+1);

    }
    //at step1
    StratumId_ThreadSpecific get_new_stratum(uint32_t size){
        if(size>=stratum_size_stats.size()){
            stratum_size_stats.resize(size+1, 0);
        }
        stratum_size_stats[size]++;
        stratum_sizes.push_back(size);
        stratum_cnt++;
        return stratum_cnt;
    }

    // repr_idx is splitted so that it is inferred.
    void add_repr_raw(SuffixArrayIdx repr_idx, StratumId_ThreadSpecific stratum_id, bool is_primary_of_rep){
        blocks_of_raws[repr_idx/block_unit].push_back(StratifiedRaw(repr_idx%block_unit, stratum_id, is_primary_of_rep));
    }

    //at the middle of step1&2 - distribute stratum ids between workspaces
    void set_real_stratum_id(vector<StratifiedRawWorkspace> &whole_workspaces, Prokrustean &prokrustean){
        /* */
        vector<StratumId> real_stratum_ids_per_local(stratum_cnt);
        uint32_t stratum_id=0;
        uint32_t idx=0;
        uint32_t max_stratum_size = stratum_size_stats.size();
        uint8_t workspace_cnt = whole_workspaces.size();
        vector<uint32_t> max_stratum_sizes_of_workspaces;
        vector<uint32_t> current_stratum_ids_per_size(max_stratum_size);
        for(auto &w: whole_workspaces){
            max_stratum_sizes_of_workspaces.push_back(w.stratum_size_stats.size());
        }
        // first, setup starting indices of stratums per size of this workspace.
        // Stratum indices are disjoint among workspaces by stat calculation
        for(uint32_t size=0; size<max_stratum_size; size++){
            for(int i=0; i<workspace_cnt; i++){
                if(size>=max_stratum_sizes_of_workspaces[i]){
                    continue;
                }
                if(i==this->workspace_id){
                    current_stratum_ids_per_size[size]=stratum_id;
                }
                stratum_id+=whole_workspaces[i].stratum_size_stats[size];
            }
        }

        // second, setup local__to__global stratum id map
        for(uint32_t i=0; i<stratum_cnt; i++){
            auto size_of_the_stratum = stratum_sizes[i];
            //for each stratum of size -> map the corresponding real stratum id
            real_stratum_ids_per_local[i]=current_stratum_ids_per_size[size_of_the_stratum];
            current_stratum_ids_per_size[size_of_the_stratum]++;
        }

        //third, set stratum sizes to prokrustean structure
        assert(stratum_cnt <= prokrustean.stratums.size());
        for(uint64_t i=0; i<stratum_cnt; i++){
            prokrustean.stratums.at(real_stratum_ids_per_local[i]).size=stratum_sizes[i];
        }

        //lastly, switch the stratum ids in the raw data
        for(auto &raws: blocks_of_raws){
            for(auto &raw: raws){
                raw.stratum_id=real_stratum_ids_per_local[raw.stratum_id];
            }
        }
    }

    void clear_local_stratum(){
        stratum_sizes.clear();
        stratum_sizes.shrink_to_fit();
        stratum_size_stats.clear();
        stratum_size_stats.shrink_to_fit();
    }

    void clear_block(int block_no){
        blocks_of_raws[block_no].clear();
        blocks_of_raws[block_no].shrink_to_fit();
    }

    // debugging purpose
    uint64_t get_cardinality(){
        uint64_t cardinality=0;
        for(int i=0;i<this->blocks_of_raws.size(); i++){
            cardinality+=this->blocks_of_raws[i].size();
        }
        return cardinality;
    }

    // debugging purpose
    vector<uint64_t> restore_reprs(){
        vector<uint64_t> reprs;
        for(int b=0;b<this->blocks_of_raws.size(); b++){
            for(auto raw: this->blocks_of_raws[b]){
                auto original = raw.repr_sa_in_block+b*block_unit;
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
struct StratifiedReprBased{
    /* repr_sa is inferred from the Block location and bit vector (repr_exists) */
    int count=0;
    // fixed length rep_id (dynamically allocated) to secure space efficiency
    StratumId *stratum_id_array;
    // fixed length whether it is primary(first of repr suffix index) of rep (dynamically allocated) to secure space efficiency
    bool *stratum_region_is_primary_array;
};

struct StratifiedBlock{
    uint32_t block_no;
    bit_vector repr_exists;
    rank_support_v<> repr_exists_rank;
    // index is inferred from outside. (repr_sa/block_size, repr_sa % block_size)
    vector<StratifiedReprBased> reprs;
    optional<StratifiedReprBased*> check_and_get_repr(uint16_t repr_sa_in_block){
        if(!repr_exists[repr_sa_in_block]) return nullopt;
        
        return &reprs[repr_exists_rank.rank(repr_sa_in_block)];
    }
    void setup(int block_unit, uint32_t block_no){
        this->block_no=block_no;
        this->repr_exists=bit_vector(block_unit, 0);
    }

    void index_reprs(vector<StratifiedRawWorkspace> &workspaces){
        // set bit vector
        for(int i=0; i<workspaces.size(); i++){
            auto raw_cnt = workspaces[i].blocks_of_raws[block_no].size();
            for(auto &raw :workspaces[i].blocks_of_raws[block_no]){
                if(!repr_exists[raw.repr_sa_in_block]){
                    repr_exists[raw.repr_sa_in_block]=true;
                }
            }
        }
        // // set ranks
        this->repr_exists_rank=rank_support_v<>(&this->repr_exists);
        auto repr_cnt = this->repr_exists_rank.rank(this->repr_exists.size());
        this->reprs = vector<StratifiedReprBased>(repr_cnt);
        // 
        vector<int> stra_cnt_per_repr(repr_cnt, 0);
        for(int i=0; i<workspaces.size(); i++){
            for(auto &raw :workspaces[i].blocks_of_raws[block_no]){
                stra_cnt_per_repr[this->repr_exists_rank.rank(raw.repr_sa_in_block)]++;
            }
        }

        for(int i=0; i<repr_cnt; i++){
            auto stra_cnt=stra_cnt_per_repr[i];
            // auto repr_rank = this->repr_exists_rank.rank(i);
            this->reprs[i].count=stra_cnt;
            this->reprs[i].stratum_id_array= (StratumId*) malloc(stra_cnt*sizeof(StratumId));
            this->reprs[i].stratum_region_is_primary_array= (bool*) malloc(stra_cnt*sizeof(bool));
        }

        vector<int> repr_stat(repr_cnt, 0);
        for(int i=0; i<workspaces.size(); i++){
            for(auto &raw :workspaces[i].blocks_of_raws[block_no]){
                auto repr_rank = this->repr_exists_rank.rank(raw.repr_sa_in_block);
                auto stra_idx = this->reprs[repr_rank].count-stra_cnt_per_repr[repr_rank];
                this->reprs[repr_rank].stratum_id_array[stra_idx]=raw.stratum_id;
                this->reprs[repr_rank].stratum_region_is_primary_array[stra_idx]=raw.is_primary_of_rep;

                stra_cnt_per_repr[repr_rank]--;
            }
        }
    }
};


struct StratifiedSA_ParallelModel {
    // each thread has own blocks of raw data that will be converted to real block
    vector<StratifiedRawWorkspace> parallel_workspaces;
    // parallel converged - boolean but to make it threadsafe, use 1 byte symbol
    vector<uint8_t> block_already_converged;
    // blocks that are each merged from correspondings in parallel workspaces
    vector<StratifiedBlock> blocks;
    // block size is 65535
    int block_unit = numeric_limits<uint16_t>::max();
    int parallel_scale;

    StratifiedSA_ParallelModel(int thread_cnt, uint64_t seq_length){
        parallel_scale=thread_cnt;
        for(int i=0; i<parallel_scale; i++){
            parallel_workspaces.push_back(StratifiedRawWorkspace(i, seq_length, block_unit));    
        }
        blocks=vector<StratifiedBlock>(seq_length/block_unit+1);
        for(int i=0; i<blocks.size(); i++){
            blocks[i].setup(block_unit, i);
        }
        block_already_converged=vector<uint8_t>(seq_length/block_unit+1, 0);
    }

    optional<StratifiedReprBased*> query(uint64_t repr_sa){
        return blocks[repr_sa/block_unit].check_and_get_repr(repr_sa%block_unit);
    }

    void converge_stratums(Prokrustean &prokrustean){
        uint64_t total_stratums=0;
        uint32_t stratum_max_size;
        for(int i=0; i< parallel_workspaces.size(); i++){
            total_stratums+=parallel_workspaces[i].stratum_cnt;
        }
        prokrustean.stratums = vector<Stratum>(total_stratums);

        for(int i=0; i< parallel_scale; i++){
            // will parallelize if too slow.
            parallel_workspaces[i].set_real_stratum_id(parallel_workspaces, prokrustean);
            parallel_workspaces[i].clear_local_stratum();
        }
    }

    void converge_block(uint64_t block_no){
        assert(block_already_converged[block_no]==0);
        block_already_converged[block_no]=1;

        blocks[block_no].index_reprs(parallel_workspaces);

        for(int i=0;i<parallel_scale; i++){
            parallel_workspaces[i].clear_block(block_no);
        }
    }
};

#endif