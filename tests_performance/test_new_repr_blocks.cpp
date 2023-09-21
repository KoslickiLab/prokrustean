#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.procedures_new.hpp"
#include "../src/fm_index/tree_new.hpp"
#include "../src/application/kmers.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

#include <iostream>
#include <atomic>
#include <thread>
#include <vector>

/* ************************************************************ */
/* At step1, threads use a list of blocks to fill in raw data   */
/* Thread each has their own set of blocks                       */
/* ************************************************************ */
struct StratifiedRaw{
    uint16_t repr_sa_in_block;
    RepId rep_id;
    bool is_primary_of_rep;

    StratifiedRaw(uint16_t repr_sa_in_block, uint64_t rep_id, bool is_primary_of_rep)
    : repr_sa_in_block(repr_sa_in_block), rep_id(rep_id), is_primary_of_rep(is_primary_of_rep) {}
};

struct StratifiedRawWorkspace{
    vector<vector<StratifiedRaw>> blocks_of_raws;
    int block_unit;
    StratifiedRawWorkspace(uint64_t seq_length, int block_unit){
        this->block_unit=block_unit;
        this->blocks_of_raws=vector<vector<StratifiedRaw>>(seq_length/block_unit+1);

    }
    // repr_idx is splitted so that it is inferred.
    void add_repr_raw(uint64_t repr_idx, uint64_t rep_id, bool is_primary_of_rep){
        blocks_of_raws[repr_idx/block_unit].push_back(StratifiedRaw(repr_idx%block_unit, rep_id, is_primary_of_rep));
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
    RepId *rep_id_array;
    // fixed length whether it is primary(first of repr suffix index) of rep (dynamically allocated) to secure space efficiency
    bool *rep_id_is_primary_array;
};

struct StratifiedBlock{
    rank_support_v<> repr_exists_rank;
    bit_vector repr_exists;
    // index is inferred from outside. (repr_sa/block_size, repr_sa % block_size)
    vector<StratifiedReprBased> reprs;

    StratifiedBlock(int block_unit){
        this->repr_exists=bit_vector(block_unit, 0);
    }
    optional<StratifiedReprBased*> check_and_get_repr(uint16_t repr_sa_in_block){
        if(!repr_exists[repr_sa_in_block]) return nullopt;
        
        return &reprs[repr_exists_rank.rank(repr_sa_in_block)];
    }

    void index_reprs(vector<vector<StratifiedRaw>*>  &raws_from_each_workspace){
        // use temporary map only locally and preserve more structured faster way.
        std::unordered_map<int, vector<StratifiedRaw*>> repr_rep;
        // set bit vector
        for(int i=0; i<raws_from_each_workspace.size(); i++){
            auto raw_cnt = (*raws_from_each_workspace[i]).size();
            for(int r=0; r<raw_cnt; r++){
                auto repr_in_block = (*raws_from_each_workspace[i])[r].repr_sa_in_block;
                if(!repr_exists[repr_in_block]){
                    repr_exists[repr_in_block]=true;
                    repr_rep[repr_in_block]=vector<StratifiedRaw*>();
                }
                repr_rep[repr_in_block].push_back(&(*raws_from_each_workspace[i])[r]);
            }
        }
        // set ranks
        this->repr_exists_rank=rank_support_v<>(&this->repr_exists);
        auto repr_cnt = this->repr_exists_rank.rank(this->repr_exists.size());
        reprs = vector<StratifiedReprBased>(repr_cnt);
        // set reprs by rep ids
        for (auto it = repr_rep.begin(); it != repr_rep.end(); ++it) {
            auto repr_rank=this->repr_exists_rank.rank(it->first);
            auto rep_cnt = it->second.size();
            reprs[repr_rank].count=rep_cnt;
            reprs[repr_rank].rep_id_array= (RepId*) malloc(rep_cnt*sizeof(RepId));
            reprs[repr_rank].rep_id_is_primary_array= (bool*) malloc(rep_cnt*sizeof(bool));
            for(int i=0; i< rep_cnt; i++){
                reprs[repr_rank].rep_id_array[i]=it->second[i]->rep_id;
                reprs[repr_rank].rep_id_is_primary_array[i]=it->second[i]->is_primary_of_rep;
            }
        }
    }
    
};

//sorting 
// #include <iostream>
// #include <algorithm>  // for std::sort

// int main() {
//     // Assume ptr points to a dynamically allocated array of N integers.
//     int N = 5;
//     int* ptr = new int[N] {5, 2, 9, 1, 4};

//     // Sort the array.
//     std::sort(ptr, ptr + N);

//     // Print the sorted array.
//     for (int i = 0; i < N; ++i) {
//         std::cout << ptr[i] << " ";
//     }
//     std::cout << std::endl;

//     // Don't forget to delete the dynamically allocated memory!
//     delete[] ptr;

//     return 0;
// }


// seq_length/numeric_limits<uint16_t>::max()
struct StratifiedSA_ParallelModel {
    // each thread has own blocks of raw data that will be converted to real block
    vector<StratifiedRawWorkspace> parallel_workspaces;
    // blocks that are each merged from correspondings in parallel workspaces
    vector<StratifiedBlock> blocks;
    // block size is 65535
    int block_unit = numeric_limits<uint16_t>::max();
    int parallel_scale;

    StratifiedSA_ParallelModel(int thread_cnt, uint64_t seq_length){
        parallel_scale=thread_cnt;
        parallel_workspaces=vector<StratifiedRawWorkspace>(parallel_scale, StratifiedRawWorkspace(seq_length, block_unit));
        blocks=vector<StratifiedBlock>(seq_length/block_unit+1, StratifiedBlock(block_unit));
    }

    optional<StratifiedReprBased*> query(uint64_t repr_sa){
        return blocks[repr_sa/block_unit].check_and_get_repr(repr_sa%block_unit);
    }

    void converge_block(uint64_t block_no){
        StratifiedBlock* output_block=&blocks[block_no];
        vector<vector<StratifiedRaw>*> raws_from_each_workspace(parallel_scale);
        // assuming each thread is assigned blocks exclusively
        for(int i=0;i<parallel_scale; i++){
            raws_from_each_workspace[i]=&parallel_workspaces[i].blocks_of_raws[block_no];
        }
        output_block->index_reprs(raws_from_each_workspace);
        for(int i=0;i<parallel_scale; i++){
            parallel_workspaces[i].clear_block(block_no);
        }
    }

    void converge_completed(){
        parallel_workspaces.clear();
        parallel_workspaces.shrink_to_fit();
    }
};


void test_stratified_raw_block_operation(){
    uint64_t seq_length=1000001; 
    int block_unit=65535;
    auto workspace=StratifiedRawWorkspace(seq_length, block_unit);
    workspace.add_repr_raw(1234, 9999, true);
    assert(workspace.get_cardinality()==1);
    assert(workspace.restore_reprs()[0]==1234);

    workspace.add_repr_raw(1000002, 9998, false);
    assert(workspace.get_cardinality()==2);
    assert(workspace.restore_reprs()[0]==1234);
    assert(workspace.restore_reprs()[1]==1000002); // all length covered.

    workspace.add_repr_raw(0, 9999, false);
    assert(workspace.get_cardinality()==3);
    assert(workspace.restore_reprs()[0]==1234);
    assert(workspace.restore_reprs()[1]==0); // same block -> 0 comes next to 1234
    assert(workspace.restore_reprs()[2]==1000002);
}

void test_stratified_block_operation(){
    uint64_t seq_length=1000001; 
    int block_unit=65535;
    auto workspace=StratifiedRawWorkspace(seq_length, block_unit);
    auto block=StratifiedBlock(block_unit);
    vector<vector<StratifiedRaw>*>  raws_from_each_workspace(1);
    raws_from_each_workspace[0]=&workspace.blocks_of_raws[0];

    workspace.add_repr_raw(1234, 9999, true);
    workspace.add_repr_raw(65534, 3333, false);
    workspace.add_repr_raw(8322, 4444, true);
    workspace.add_repr_raw(8322, 5555, false);
    block.index_reprs(raws_from_each_workspace);
    auto result = block.check_and_get_repr(1234);
    assert(result.has_value());
    assert(result.value()->count==1);
    assert(result.value()->rep_id_array[0]==9999);
    assert(result.value()->rep_id_is_primary_array[0]);

    result = block.check_and_get_repr(65534);
    assert(result.value()->rep_id_array[0]==3333);
    assert(result.value()->rep_id_is_primary_array[0]==false);

    result = block.check_and_get_repr(8322);
    assert(result.value()->rep_id_array[0]==4444);
    assert(result.value()->rep_id_is_primary_array[0]==true);
    result = block.check_and_get_repr(8322);
    assert(result.value()->rep_id_array[1]==5555);
    assert(result.value()->rep_id_is_primary_array[1]==false);
}

void test_block_operation_with_tree(){
    int Lmin = 30;
    // auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');

    auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto start = std::chrono::steady_clock::now();
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << " wv tree took " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats;
    SuffixArrayNode_NEW root_new = get_root_new(fm_idx);
    vector<MaximalRepeatAnnotation> repeats_new;
    
    start = std::chrono::steady_clock::now();
    // navigate_tree_new<MaximalRepeatAnnotation, get_repeat_annotations_new>(root_new, Lmin, fm_idx, repeats_new);
    navigate_tree_new<MaximalRepeatAnnotation, report_repr_sa>(root_new, Lmin, fm_idx, repeats_new);
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    int repr_cnt=0;
    int repr_distinct_cnt=0;
    auto repr_distinct = vector<bool>(str.size());
    for(int i=0; i<repeats_new.size(); i++){
        for(auto idx: repeats_new[i].repr_indexes){
            repr_cnt++;
            if(!repr_distinct[idx]){
                repr_distinct_cnt++;
                repr_distinct[idx]=true;
            }
        }
    }
    cout << "repeats reprs cnt: "<< repr_cnt <<", distinct count: " << repr_distinct_cnt  << endl;

    start = std::chrono::steady_clock::now();
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    cout << "original algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    assert(repeats.size()==repeats_new.size());
    cout << "repeats are same: " << repeats.size() << endl;
    
    for(int i=0; i<repeats.size(); i++){
        assert(repeats[i].repr_indexes.size()==repeats_new[i].repr_indexes.size());
    }
    cout << "repeats reprs are same" << endl;
}


void main_performance_new_repr_block() {
    test_stratified_raw_block_operation();
    test_stratified_block_operation();
    // test_memory_reserved_stack();
    // test_variable_length_array();
}
