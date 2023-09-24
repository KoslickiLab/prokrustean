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
    atomic<StratumId> id_generator;
    uint64_t seq_length=1000001; 
    auto workspace=StratifiedRawWorkspace(0, seq_length, id_generator);
    workspace.add_repr_raw(1234, 9999, true);
    assert(workspace.stratified_rgn_cnt==1);
    assert(workspace.restore_reprs()[0]==1234);

    workspace.add_repr_raw(1000002, 9998, false);
    assert(workspace.stratified_rgn_cnt==2);
    assert(workspace.restore_reprs()[0]==1234);
    assert(workspace.restore_reprs()[1]==1000002); // all length covered.

    workspace.add_repr_raw(0, 9999, false);
    assert(workspace.stratified_rgn_cnt==3);
    assert(workspace.restore_reprs()[0]==1234);
    assert(workspace.restore_reprs()[1]==0); // same block -> 0 comes next to 1234
    assert(workspace.restore_reprs()[2]==1000002);
}

void test_stratum_convergence(){
    uint64_t seq_length=1000001; 
    vector<StratifiedRawWorkspace> workspaces;
    atomic<StratumId> id_generator;
    Prokrustean prokrustean;
    workspaces.push_back(StratifiedRawWorkspace(0, seq_length, id_generator));
    workspaces.push_back(StratifiedRawWorkspace(1, seq_length, id_generator));
    workspaces.push_back(StratifiedRawWorkspace(2, seq_length, id_generator));
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
        w.set_prokrutean_stratum(prokrustean);
    }
    for(auto &stratum: prokrustean.stratums){
        assert(stratum.size>0);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> repr_dist(0, seq_length);
    std::uniform_int_distribution<uint32_t> size_dist(1, seq_length/1000);
    workspaces.clear();
    id_generator.exchange(0);
    prokrustean=Prokrustean();
    for(int i=0; i<10; i++){
        workspaces.push_back(StratifiedRawWorkspace(i, seq_length, id_generator));
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
        w.set_prokrutean_stratum(prokrustean);
    }
    
    for(auto &stratum: prokrustean.stratums){
        assert(stratum.size>0);
    }
}

void test_stratified_block_convergence(){
    uint64_t seq_length=324152; 
    uint64_t total_cnt=0;
    atomic<StratumId> id_generator;
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
    model.prepare_prokrustean(prokrustean);
    uint64_t total_cnt2=0;
    for(int i=0; i<model.blocks.size(); i++){
        model.converge_block(i);
        model.blocks[i].validate();
        total_cnt2+=model.blocks[i].stat__rgn_cnt;
    }

    assert(total_cnt==total_cnt2);
}

void test_parallelized_push(){
    /* Example result
        seq N:10billion
        stratum length:5000
        stratums per thread:9259259 --> this is the key info
        repr per stratum:3
        step2 sample created: 33837ms
        step2 stratum convergence: 7644ms
        step2 block convergence start 
        step2 block convergence: 16438ms
    */
    auto num_threads = 12;
    // performance- influencing factors
    uint64_t how_many_stratums_sampled_per_thread=pow(10, 8)/num_threads; //1 billion total
    uint64_t how_many_reprs_sampled_per_stratum=3;

    uint64_t seq_length=how_many_stratums_sampled_per_thread*20;
    uint64_t max_stratum_length=5000; 
    vector<future<void>> futures;
    auto model=StratifiedSA_ParallelModel(num_threads, seq_length);
    auto start = std::chrono::steady_clock::now();

    vector<uint64_t> params={seq_length, max_stratum_length, how_many_stratums_sampled_per_thread, how_many_reprs_sampled_per_stratum};
    
    cout << "performance push: "<<endl; 
    cout << "stratums per thread:" << params[2]<<endl; 
    cout << "repr per stratum:" << params[3]<<endl; 

    // push dummy stratifieds
    auto func__ = [](StratifiedSA_ParallelModel &model, int workspace_id, vector<uint64_t> params) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<uint32_t> repr_dist(0, params[0]-1);
            std::uniform_int_distribution<uint32_t> size_dist(1, params[1]);
            auto w = &model.parallel_workspaces[workspace_id];
            auto no = params[2];
            for(int i=0; i<no; i++){
                auto stratum_id=w->get_new_stratum(size_dist(gen));
                for(int r=0; r<params[3]; r++){
                    auto repr = repr_dist(gen);
                    w->add_repr_raw(repr, stratum_id, true);
                }
            }
        };
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__, ref(model), i, ref(params)));
    }
    for (auto& future : futures) {
        future.wait();
    }
    futures.clear();
    cout << "step2 sample created: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
    start = std::chrono::steady_clock::now();

    /* parallel version! - space efficiency issue */
    Prokrustean prokrustean;
    int stratum_cnt=0;
    for(auto &work: model.parallel_workspaces){
        stratum_cnt+=work.stratum_cnt;
    }
    prokrustean.stratums = vector<Stratum>(stratum_cnt);

    auto func__converge_stra = [](StratifiedRawWorkspace &workspace, Prokrustean &prokrustean) {
        workspace.set_prokrutean_stratum(prokrustean);
    };
    for(int i=0; i<num_threads; i++){
        // auto f = std::async(std::launch::async, collect_repeat_thread, ref(subtree_roots), Lmin, ref(fm_idx), ref(m));
        futures.push_back(std::async(std::launch::async, func__converge_stra, ref(model.parallel_workspaces[i]), ref(prokrustean)));
    }
    for (auto& future : futures) {
        future.wait();
    }
    futures.clear();
    
    cout << "step2 stratum convergence: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    cout << "sleeping... " << endl;
    std::this_thread::sleep_for(std::chrono::seconds(20));

    cout << "step2 block convergence start " << endl;
    start = std::chrono::steady_clock::now();
    auto func__converge_block = [](StratifiedSA_ParallelModel &model, std::atomic<int> &block_idx_counter) {
        while(true){
            auto idx = block_idx_counter.fetch_add(1);
            if(idx >= model.blocks.size()){
                break;
            }
            model.converge_block(idx);
            // std::this_thread::sleep_for(std::chrono::milliseconds(100));
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
    cout << "step2 block convergence: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // std::this_thread::sleep_for(std::chrono::seconds(20));
    
    uint64_t total_cnt=0;
    uint64_t total_cnt2=0;
    for(auto &work: model.parallel_workspaces){
        total_cnt+=work.stratified_rgn_cnt;
    }
    for(auto &block: model.blocks){
        block.validate();
        total_cnt2+=block.stat__rgn_cnt;
    }
    cout << "total " << total_cnt << " total2 " << total_cnt2 << endl;
    assert(total_cnt==total_cnt2);
}


void test_memory(){

    const size_t num_blocks = 100000;
    const size_t block_size = 1024 * 1024;  // 1 MiB

    // Allocate memory
    std::cout << "Allocating memory..." << std::endl;
    char* blocks[num_blocks];
    for (size_t i = 0; i < num_blocks; ++i) {
        blocks[i] = (char*)malloc(block_size);
        if (!blocks[i]) {
            std::cerr << "Memory allocation failed." << std::endl;
        }
    }

    // Sleep to allow for memory usage observation
    std::cout << "Sleeping for 10 seconds..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));

    // Deallocate memory
    std::cout << "Deallocating memory..." << std::endl;
    for (size_t i = 0; i < num_blocks; ++i) {
        free(blocks[i]);
    }

    // Sleep again to allow for memory usage observation
    std::cout << "Sleeping for 10 seconds..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));

    std::cout << "Done." << std::endl;
}
struct Block{
    uint32_t cnt=0;
     // make private for better memory usage
    vector<StratifiedRaw> raw_regions;
};

void test_memory2(){

    const size_t num_blocks = 100000000;

    // Allocate memory
    std::cout << "Allocating memory..." << std::endl;
    vector<Block> blocks(num_blocks);
    for(auto &block: blocks){
        block.raw_regions.push_back(StratifiedRaw(1, 1, true));
    }
    // Block* blocks[num_blocks];
    // for (size_t i = 0; i < num_blocks; ++i) {
    //     blocks[i] = (Block*)malloc(sizeof(Block));
    //     if (!blocks[i]) {
    //         std::cerr << "Memory allocation failed." << std::endl;
    //     }
    // }

    // Sleep to allow for memory usage observation
    std::cout << "Sleeping for 10 seconds..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));

    // Deallocate memory
    std::cout << "Deallocating memory..." << std::endl;
    blocks.clear();
    blocks.shrink_to_fit();
    blocks=vector<Block>();
    // for (size_t i = 0; i < num_blocks; ++i) {
    //     free(blocks[i]);
    // }

    // Sleep again to allow for memory usage observation
    std::cout << "Sleeping for 10 seconds..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(20));

    std::cout << "Done." << std::endl;
}

void test_memory3(){
    std::cout << sizeof(tuple<uint64_t, bool>) << std::endl;
}

void main_performance_new_repr_block() {
    // test_stratified_raw_block_operation();
    // test_stratum_convergence();
    // test_stratified_block_convergence();
    test_parallelized_push();
    // test_memory();
    // test_memory2();
    // test_memory3();
}
