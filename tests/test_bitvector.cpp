#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/BitMagic/bm.h"

using namespace std;

void test_bitvector_ranks(){
    bm::bvector<>   bv { 1, 20, 30, 31 }; // init a test bit-vector
        
    // construct rank-select index
    std::unique_ptr<bm::bvector<>::rs_index_type>
                                rs_idx(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs_idx.get());

    // lets find a few ranks here
    //
    auto r1 = bv.rank(20, *rs_idx);
    assert(r1==2);
    r1 = bv.rank(21, *rs_idx);
    assert(r1==2);
    r1 = bv.rank(30, *rs_idx);
    assert(r1==3);
    // position value corrected rank
    // one special case of rank function returns rank-1
    // if position bit is set or just rank, otherwise
    //
    // this is an equivalent of
    // bv.count_range(0, n) - bv.text(n)
    // (just faster, because of the fused rank-test)
    //
    auto r1c = bv.rank_corrected(31, *rs_idx); // 3
    assert(r1c==3);
    r1c = bv.rank_corrected(32, *rs_idx); // 4
    assert(r1c==4);
    r1c = bv.rank_corrected(33, *rs_idx); // 4
    assert(r1c==4);

    // now perform a search for a position for a rank
    //
    bm::bvector<>::size_type pos;
    bool found = bv.select(2, pos, *rs_idx);
    assert(found && pos==20);
    found = bv.select(3, pos, *rs_idx);
    assert(found && pos==30);
    found = bv.select(5, pos, *rs_idx);
    assert(~found);
}

void test_fixed_size_bitvector(){
    bm::bvector<> bv(100);
    bv.set_bit(10, true);
    std::unique_ptr<bm::bvector<>::rs_index_type>
                                rs_idx(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs_idx.get());

    assert(bv.rank(9, *rs_idx)==false);
    assert(bv.rank(10, *rs_idx)==true);
    assert(bv.rank(11, *rs_idx)==true);
}

void test_large_size_bitvector(){
    uint64_t N = 80000000000;
    uint64_t test_count = 10000;
    bm::bvector<> bv(N);
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
        bv.set_bit(no);
    }

    if(print) std::cout << "set bits(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
    start = std::chrono::steady_clock::now();

    std::unique_ptr<bm::bvector<>::rs_index_type>
                                rs_idx(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs_idx.get());

    if(print) std::cout << "indexing(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
    start = std::chrono::steady_clock::now();

    for (auto no: numbers){
        assert(bv.rank(no-1, *rs_idx)<bv.rank(no, *rs_idx));
    }

    if(print) std::cout << "check(ms)=" << (std::chrono::steady_clock::now()-start).count()/1000000 << std::endl;
    
    // assert(bv.rank(size-3, *rs_idx)==false);
    // assert(bv.rank(size-2, *rs_idx)==true);
    // assert(bv.rank(size-1, *rs_idx)==true);
}

void main_bitvector() {
    cout << "test_ranks finished" << endl; 
    test_bitvector_ranks();
    test_fixed_size_bitvector();
    test_large_size_bitvector();
}
