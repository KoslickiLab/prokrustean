#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;


void test_bitvector_ranks(){
    bit_vector bv(10, 0);
    bv[5]=true;

    rank_support_v<> rb(&bv);
    assert(rb.rank(4)==false);
    assert(rb.rank(5)==false);
    assert(rb.rank(6)==true);
}


void test_many_sets(){
    // N=8b, test=1b => sampling 21s, set bits 83s, get bits 66s, rank 85s
    // N=8b, test=1m => sampling 0.021s, set bits 0.071s, get bits 0.067s, rank 0.085s 
    uint64_t N = 8*pow(10,6);;
    uint64_t test_count = pow(10,6);;
    bit_vector bv(N, 0);
    vector<int> numbers;
    bool print = false;
    
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < test_count; i++){
        int no = rand() % N;
        if(no == 0){
            continue;
        }
        numbers.push_back(no);
    }

    if(print) std::cout << "sampling(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
    start = std::chrono::steady_clock::now();

    for(auto no: numbers){
        bv[no]=true;
    }

    if(print) std::cout << "set bits(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;

    start = std::chrono::steady_clock::now();

    int cnt = 0;
    for(auto no: numbers){
        if(bv[no]) cnt++;
    }
    
    if(print) std::cout << "get bits(ms)= (" << cnt << ") " << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;

    start = std::chrono::steady_clock::now();

    rank_support_v<> rb(&bv);

    if(print) std::cout << "indexing(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
    start = std::chrono::steady_clock::now();

    for (auto no: numbers){
        assert(rb.rank(no)<rb.rank(no+1));
    }

    if(print) std::cout << "rank(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
}


void main_bitvector() {
    test_bitvector_ranks();
    test_many_sets();
}
