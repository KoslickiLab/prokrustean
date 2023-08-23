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

void main_bitvector() {
    test_bitvector_ranks();
}
