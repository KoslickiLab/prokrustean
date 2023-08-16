#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"

using namespace std;
using namespace sdsl;

void test_repr_rank(){
    uint64_t idx1 = 5;
    uint64_t idx2 = 98;
    uint64_t idx3 = 99;
    vector<MaximalRepeatAnnotation> annots;
    annots.push_back({123, {idx1}});
    annots.push_back({456, {idx3, idx2}});
    ReprSuffixRank rank;
    rank.initialize(100, annots);
    assert(rank.exists(idx1) == 1);
    assert(rank.exists(idx2) == 1);
    assert(rank.exists(idx3) == 1);
    
    assert(rank.get_repr_size() == 3);

    assert(rank.rank(idx1-1) == 0);
    assert(rank.rank(idx1) == 0);
    assert(rank.rank(idx1+1) == 1);

}


void main_construction_reprrank() {
    test_repr_rank();
}
