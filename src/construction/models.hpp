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

typedef uint16_t StratumId_ThreadSpecific;
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
    StratumId_ThreadSpecific stratum_id_in_thread;
    bool is_primary_of_rep;

    StratifiedRaw(ReprIdx_BlockSpecific repr_id, StratumId_ThreadSpecific stratum_id, bool is_primary_of_rep)
    : repr_sa_in_block(repr_id), stratum_id_in_thread(stratum_id), is_primary_of_rep(is_primary_of_rep) {}
};

struct StratifiedRawWorkspace{
    uint8_t workspace_id;
    int block_unit;
    vector<vector<StratifiedRaw>> blocks_of_raws;
    // stratum id -> index, only store size 
    vector<StratumSize> stratum_sizes;
    // each size -> index, store the aubndance of size
    vector<uint32_t> stratum_size_stats;
    // 
    uint32_t stratum_cnt=0;

    vector<StratumId> real_stratum_ids;

    StratifiedRawWorkspace(uint8_t workspace_id, uint64_t seq_length, int block_unit){
        this->workspace_id=workspace_id;
        this->block_unit=block_unit;
        this->blocks_of_raws=vector<vector<StratifiedRaw>>(seq_length/block_unit+1);

    }

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

    void set_real_stratum_id(vector<StratifiedRawWorkspace> &whole_workspaces){
        /* */
        this->real_stratum_ids=vector<StratumId>(stratum_cnt);
        uint32_t stratum_id=0;
        uint32_t idx=0;
        uint32_t max_stratum_size = stratum_size_stats.size();
        uint8_t workspace_cnt = whole_workspaces.size();
        vector<uint32_t> max_stratum_sizes_of_workspaces;
        vector<uint32_t> this_workspace__curr_stratum_ids_per_size(max_stratum_size);
        for(auto &w: whole_workspaces){
            max_stratum_sizes_of_workspaces.push_back(w.stratum_size_stats.size());
        }
        
        for(uint32_t size=0; size<max_stratum_size; size++){
            for(int i=0; i<workspace_cnt; i++){
                if(size>=max_stratum_sizes_of_workspaces[i]){
                    continue;
                }
                if(i==this->workspace_id){
                    this_workspace__curr_stratum_ids_per_size[size]=stratum_id;
                }
                stratum_id+=whole_workspaces[i].stratum_size_stats[size];
            }
        }
        for(uint32_t i=0; i<stratum_cnt; i++){
            auto size_of_the_stratum = stratum_sizes[i];
            //for each stratum of size -> map the corresponding real stratum id
            real_stratum_ids[i]=this_workspace__curr_stratum_ids_per_size[size_of_the_stratum];
            this_workspace__curr_stratum_ids_per_size[size_of_the_stratum]++;
        }
    }

    void set_stratum_sizes(Prokrustean &prokrustean){
        uint64_t max_stratum_size = stratum_sizes.size();
        assert(max_stratum_size <= prokrustean.stratums.size());
        for(uint64_t i=0; i<max_stratum_size; i++){
            prokrustean.stratums.at(real_stratum_ids[i]).size=stratum_sizes[i];
        }
    }
    // void set_real_stratum_id(vector<uint32_t> &stratum_starting_id_by_size){
    //     real_stratum_ids=vector<StratumId>(temporary_stratum_sizes.size());
    //     for(int size=0; size< temporary_stratum_sizes.size()){
    //         for(int i=0; i<size; i++){
    //             stratum_starting_id_by_size[size]
    //         }
    //     }
    // }

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
    StratumId *rep_id_array;
    // fixed length whether it is primary(first of repr suffix index) of rep (dynamically allocated) to secure space efficiency
    bool *rep_id_is_primary_array;
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
        // temporary map to prepare for more structured storage.
        std::unordered_map<int, vector<tuple<StratumId, bool>>> repr_to__stratum_is_primary_pair;
        // set bit vector
        for(int i=0; i<workspaces.size(); i++){
            auto raw_cnt = workspaces[i].blocks_of_raws[block_no].size();
            for(int r=0; r<raw_cnt; r++){
                StratifiedRaw* raw = &workspaces[i].blocks_of_raws[block_no][r];
                if(!repr_exists[raw->repr_sa_in_block]){
                    repr_exists[raw->repr_sa_in_block]=true;
                    repr_to__stratum_is_primary_pair[raw->repr_sa_in_block]=vector<tuple<StratumId, bool>>();
                }
                StratumId real_id = workspaces[i].real_stratum_ids.at(raw->stratum_id_in_thread);
                repr_to__stratum_is_primary_pair[raw->repr_sa_in_block].push_back(make_tuple(real_id, raw->is_primary_of_rep));
            }
        }
        assert(false);
        // set ranks
        this->repr_exists_rank=rank_support_v<>(&this->repr_exists);
        auto repr_cnt = this->repr_exists_rank.rank(this->repr_exists.size());
        this->reprs = vector<StratifiedReprBased>(repr_cnt);
        // set reprs by rep ids
        for (auto it = repr_to__stratum_is_primary_pair.begin(); it != repr_to__stratum_is_primary_pair.end(); ++it) {
            auto repr_rank=this->repr_exists_rank.rank(it->first);
            auto rep_cnt = it->second.size();
            // to optimize space, manually assign memory
            this->reprs[repr_rank].count=rep_cnt;
            this->reprs[repr_rank].rep_id_array= (StratumId*) malloc(rep_cnt*sizeof(StratumId));
            this->reprs[repr_rank].rep_id_is_primary_array= (bool*) malloc(rep_cnt*sizeof(bool));
            for(int i=0; i< rep_cnt; i++){
                // StratumId real_id = get<0>(it->second[i]); 
                // bool is_primary = get<1>(it->second[i]);
                // this->reprs[repr_rank].rep_id_array[i]=real_id;
                // this->reprs[repr_rank].rep_id_is_primary_array[i]=is_primary;
            }
        }
    }  
};


struct StratifiedSA_ParallelModel {
    Prokrustean* prokrustean;
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

        for(int i=0; i< parallel_workspaces.size(); i++){
            // will parallelize if too slow.
            parallel_workspaces[i].set_real_stratum_id(parallel_workspaces);
            parallel_workspaces[i].set_stratum_sizes(prokrustean);
        }
    }

    void converge_block(uint64_t block_no){
        assert(block_already_converged[block_no]==0);
        block_already_converged[block_no]=1;

        StratifiedBlock* output_block=&blocks[block_no];
        vector<vector<StratifiedRaw>*> raws_from_each_workspace(parallel_scale);
        // assuming each thread is assigned blocks exclusively
        for(int i=0;i<parallel_scale; i++){
            raws_from_each_workspace[i]=&parallel_workspaces[i].blocks_of_raws[block_no];
        }
        output_block->index_reprs(parallel_workspaces);
        for(int i=0;i<parallel_scale; i++){
            parallel_workspaces[i].clear_block(block_no);
        }
    }

    void converge_completed(){
        parallel_workspaces.clear();
        parallel_workspaces.shrink_to_fit();
    }
};

#endif