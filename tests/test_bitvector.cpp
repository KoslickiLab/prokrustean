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

void main_bitvector() {
    cout << "test_ranks finished" << endl; 
    test_bitvector_ranks();
}
