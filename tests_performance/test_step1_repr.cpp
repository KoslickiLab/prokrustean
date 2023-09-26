#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "util.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/construction/models.hpp"
#include "../src/construction/algorithms.procedures_new.hpp"


using namespace std;
using namespace sdsl;


void test_stratification_lock_push(){
    int stratified_region_factor=6;
    int num_threads=12;
    int cnt=pow(10,6)/num_threads;
    vector<StratificationOutput> blocks(1000000); 
    auto func__protect = [](int cnt, vector<ProjectedStratification> &block) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> block_dist(0, block.size()-1);
        
        for(int i=0; i<cnt; i++){
            auto idx = block_dist(gen);
            block[idx].add_stratified_regions(1,1,true);
        }
    };
    auto start = std::chrono::steady_clock::now();
    vector<future<void>> futures;
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__protect, cnt, ref(blocks)));
    }
    for(int i=0; i<num_threads; i++){
        futures[i].get();
    }
    
    cout << "protected call created: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
}

void main_performance_full() {
    test_stratification_lock_push();
}
