#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
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


void test_stratified_raw_block_operation(){
    uint64_t seq_length=1000001; 
    int block_unit=65535;
    auto workspace=StratifiedRawWorkspace(0, seq_length, block_unit);
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

void test_stratum_convergence(){
    uint64_t seq_length=1000001; 
    int block_unit=65535;
    vector<StratifiedRawWorkspace> workspaces;
    workspaces.push_back(StratifiedRawWorkspace(0, seq_length, block_unit));
    auto workspace=&workspaces[0];
    (*workspace).add_repr_raw(1234, 9999, true);
    (*workspace).add_repr_raw(65534, 3333, false);
    (*workspace).add_repr_raw(8322, 4444, true);
    (*workspace).add_repr_raw(8322, 5555, false);
    workspace->set_real_stratum_id(workspaces);

    workspaces.clear();
    Prokrustean prokrustean;
    workspaces.push_back(StratifiedRawWorkspace(0, seq_length, block_unit));
    workspaces.push_back(StratifiedRawWorkspace(1, seq_length, block_unit));
    workspaces.push_back(StratifiedRawWorkspace(2, seq_length, block_unit));
    auto sid = workspaces[0].get_new_stratum(500);
    workspaces[0].add_repr_raw(1234, sid, true);
    sid = workspaces[0].get_new_stratum(400);
    workspaces[0].add_repr_raw(65534, sid, false);
    sid = workspaces[0].get_new_stratum(300);
    workspaces[0].add_repr_raw(8322, sid, true);
    workspaces[0].add_repr_raw(9999, sid, false);

    sid = workspaces[1].get_new_stratum(300);
    workspaces[1].add_repr_raw(1234, sid, true);
    sid = workspaces[1].get_new_stratum(300);
    workspaces[1].add_repr_raw(2222, sid, false);
    sid = workspaces[1].get_new_stratum(200);
    workspaces[1].add_repr_raw(4444, sid, true);
    sid = workspaces[1].get_new_stratum(100);
    workspaces[1].add_repr_raw(6666, sid, false);

    sid = workspaces[2].get_new_stratum(500);
    workspaces[2].add_repr_raw(1234, sid, true);
    sid = workspaces[2].get_new_stratum(500);
    workspaces[2].add_repr_raw(8322, sid, false);
    sid = workspaces[2].get_new_stratum(700);
    workspaces[2].add_repr_raw(8322, sid, true);
    int total_cnt = 0;
    for(auto &w: workspaces){
        total_cnt += w.stratum_cnt;
    }
    prokrustean.stratums.resize(total_cnt);

    for(auto &w: workspaces){
        w.set_real_stratum_id(workspaces);
        w.set_stratum_sizes(prokrustean);
    }
    for(auto &stratum: prokrustean.stratums){
        assert(stratum.size>0);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> repr_dist(0, seq_length);
    std::uniform_int_distribution<uint32_t> size_dist(1, seq_length/1000);
    workspaces.clear();
    prokrustean=Prokrustean();
    for(int i=0; i<10; i++){
        workspaces.push_back(StratifiedRawWorkspace(i, seq_length, block_unit));
    }
    for(auto &w: workspaces){
        for(int i=0; i<seq_length/workspaces.size()/10; i++){
            auto stratum_id=w.get_new_stratum(size_dist(gen));
            for(int r=0; r<10; r++){
                w.add_repr_raw(repr_dist(gen), stratum_id, true);
            }
        }
    }
    total_cnt = 0;
    for(auto &w: workspaces){
        total_cnt += w.stratum_cnt;
    }

    prokrustean.stratums.resize(total_cnt);
    for(auto &w: workspaces){
        w.set_real_stratum_id(workspaces);
        w.set_stratum_sizes(prokrustean);
    }

    for(auto &stratum: prokrustean.stratums){
        assert(stratum.size>0);
    }
}

void test_stratified_block_operation(){
    uint64_t seq_length=1000001; 
    int block_unit=65535;
    vector<StratifiedRawWorkspace> workspaces;
    workspaces.push_back(StratifiedRawWorkspace(0, seq_length, block_unit));
    auto workspace=&workspaces[0];
    (*workspace).add_repr_raw(1234, 9999, true);
    (*workspace).add_repr_raw(65534, 3333, false);
    (*workspace).add_repr_raw(8322, 4444, true);
    (*workspace).add_repr_raw(8322, 5555, false);
    workspace->set_real_stratum_id(workspaces);

    auto block=StratifiedBlock();
    block.setup(block_unit, 0);
    block.index_reprs(workspaces);
    // auto result = block.check_and_get_repr(1234);
    // assert(result.has_value());
    // assert(result.value()->count==1);
    // assert(result.value()->rep_id_array[0]==9999);
    // assert(result.value()->rep_id_is_primary_array[0]);

    // result = block.check_and_get_repr(65534);
    // assert(result.value()->rep_id_array[0]==3333);
    // assert(result.value()->rep_id_is_primary_array[0]==false);

    // result = block.check_and_get_repr(8322);
    // assert(result.value()->rep_id_array[0]==4444);
    // assert(result.value()->rep_id_is_primary_array[0]==true);
    // result = block.check_and_get_repr(8322);
    // assert(result.value()->rep_id_array[1]==5555);
    // assert(result.value()->rep_id_is_primary_array[1]==false);
}

void test_stratified_parallel_model_coverage(){
    int thread = 3;
    uint64_t seq_length=1000001; 
    auto model = StratifiedSA_ParallelModel(thread, seq_length);
    auto workspace1=&model.parallel_workspaces[0];
    (*workspace1).add_repr_raw(1234, 9999, true);
    (*workspace1).add_repr_raw(65534, 3333, false);
    (*workspace1).add_repr_raw(8322, 4444, true);
    (*workspace1).add_repr_raw(8322, 5555, false);
    model.converge_block(1234/model.block_unit);
    auto result = model.query(1234);
    assert(result.has_value());
    assert(result.value()->count==1);
    assert(result.value()->rep_id_array[0]==9999);
    assert(result.value()->rep_id_is_primary_array[0]);
    result = model.query(8322);
    assert(result.has_value());
    assert(result.value()->count==2);
    assert(result.value()->rep_id_array[0]==4444);
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
    test_stratum_convergence();
    // test_stratified_block_operation();
    // test_stratified_parallel_model_coverage();
    
    // test_memory_reserved_stack();
    // test_variable_length_array();
}
