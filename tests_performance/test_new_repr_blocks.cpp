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
        w.set_real_stratum_id(workspaces, prokrustean);
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
        w.set_real_stratum_id(workspaces, prokrustean);
    }

    for(auto &stratum: prokrustean.stratums){
        assert(stratum.size>0);
    }
}

void test_stratified_block_convergence(){
    uint64_t seq_length=324152; 
    uint64_t total_cnt=0;
    std::unordered_map<SuffixArrayIdx, vector<StratumId>> data;
    auto model=StratifiedSA_ParallelModel(10, seq_length);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> repr_dist(0, seq_length-1);
    std::uniform_int_distribution<uint32_t> size_dist(1, seq_length/1000);
    for(auto &w: model.parallel_workspaces){
        for(int i=0; i<seq_length/model.parallel_workspaces.size()/10; i++){
            auto stratum_id=w.get_new_stratum(size_dist(gen));
            for(int r=0; r<10; r++){
                auto repr = repr_dist(gen);
                w.add_repr_raw(repr, stratum_id, true);
                if(data.find(repr) == data.end()){
                    data[repr]=vector<StratumId>();
                }
                data[repr].push_back(stratum_id);
                total_cnt++;
            }
        }
    }

    Prokrustean prokrustean;
    model.converge_stratums(prokrustean);
    for(int i=0; i<model.blocks.size(); i++){
        model.converge_block(i);
    }

    cout << "total:" << total_cnt << endl;

    uint64_t total_cnt2=0;
    for(uint64_t i=0; i<seq_length; i++){
        auto ret=model.query(i);
        if(ret.has_value()){
            total_cnt2+=ret.value()->count;
            assert(data[i].size()==ret.value()->count);
        }
    }
    cout << "total2:" << total_cnt2 << endl;
}

void test_parallelized_push(){
    uint64_t seq_length=32415200; 
    uint64_t total_cnt=0;
    auto num_threads = 10;
    vector<future<void>> futures;
    auto model=StratifiedSA_ParallelModel(num_threads, seq_length);
    
    // push dummy stratifieds
    auto func__ = [](StratifiedSA_ParallelModel &model, int workspace_id, uint64_t seq_length) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<uint32_t> repr_dist(0, seq_length-1);
            std::uniform_int_distribution<uint32_t> size_dist(1, seq_length/1000);
            auto w = &model.parallel_workspaces[workspace_id];
            for(int i=0; i<seq_length/model.parallel_workspaces.size()/10; i++){
                auto stratum_id=w->get_new_stratum(size_dist(gen));
                for(int r=0; r<10; r++){
                    auto repr = repr_dist(gen);
                    w->add_repr_raw(repr, stratum_id, true);
                }
            }
        };
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__, ref(model), i, seq_length));
    }
    for (auto& future : futures) {
        future.wait();
    }
    futures.clear();

    auto start = std::chrono::steady_clock::now();
    Prokrustean prokrustean;
    uint32_t stratum_max_size;
    for(int i=0; i< model.parallel_workspaces.size(); i++){
        for(auto &raws: model.parallel_workspaces[i].blocks_of_raws){
            total_cnt+=raws.size();
        }
    }
    prokrustean.stratums = vector<Stratum>(total_cnt);

    auto func__converge_stra = [](StratifiedSA_ParallelModel &model, Prokrustean &prokrustean, int workspace_id) {
        model.parallel_workspaces[workspace_id].set_real_stratum_id(model.parallel_workspaces, prokrustean);
    };
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__converge_stra, ref(model), ref(prokrustean), i));
    }
    for (auto& future : futures) {
        future.wait();
    }
    futures.clear();
    for(int i=0; i<num_threads; i++){
        model.parallel_workspaces[i].clear_local_stratum();
    }
    cout << "step2 stratum convergence: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();



    auto func__converge_block = [](StratifiedSA_ParallelModel &model, std::atomic<int> &block_idx_counter) {
        while(true){
            auto idx = block_idx_counter.fetch_add(1);
            if(idx >= model.blocks.size()){
                break;
            }
            model.converge_block(idx);
        }
    };
    std::atomic<int> block_idx_counter;
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__converge_block, ref(model), ref(block_idx_counter)));
    }
    for (auto& future : futures) {
        future.wait();
    }

    uint64_t total_cnt2=0;
    for(uint64_t i=0; i<seq_length; i++){
        auto ret=model.query(i);
        if(ret.has_value()){
            total_cnt2+=ret.value()->count;
        }
    }
    assert(total_cnt==total_cnt2);
    cout << "total:" << total_cnt2 << endl;
    cout << "step2 block convergence: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}


void main_performance_new_repr_block() {
    // test_stratified_raw_block_operation();
    // test_stratum_convergence();
    // test_stratified_block_convergence();
    test_parallelized_push();
}
