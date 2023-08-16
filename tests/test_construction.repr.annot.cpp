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
    uint64_t N = 100;
    uint64_t idx1 = 5;
    uint64_t idx2 = 98;
    uint64_t idx3 = 99;
    vector<MaximalRepeatAnnotation> repeats;
    repeats.push_back({123, {idx1}});
    repeats.push_back({456, {idx3, idx2}});

    ReprSuffixRank rank;
    rank.initialize(N, repeats);
    assert(rank.exists(idx1) == 1);
    assert(rank.exists(idx2) == 1);
    assert(rank.exists(idx3) == 1);
    assert(rank.get_repr_size() == 3);
    assert(rank.rank(idx1-1) == 0);
    assert(rank.rank(idx1) == 0);
    assert(rank.rank(idx1+1) == 1);
}

void test_repr_annotation(){
    uint64_t N = 100;
    uint64_t idx1 = 5;
    uint64_t idx2 = 10;
    uint64_t idx3 = 50;
    uint64_t idx_empty = 66;
    RepId rep1 = 2;
    RepId rep2 = 4;
    
    vector<MaximalRepeatAnnotation> repeats;
    repeats.push_back({rep1, {idx2}});
    repeats.push_back({rep2, {idx1, idx2, idx3}});

    ReprSuffixAnnotation annot;
    annot.initialize_rank(N, repeats);
    annot.initialize_repr_sa(repeats);
    
    assert(annot.get_repeats(idx_empty).has_value() == 0);
    assert(annot.get_repeats(idx1).has_value() == 1);
    assert(annot.get_repeats(idx1).value().size() == 1);
    assert(annot.get_repeats(idx2).value().size() == 2);
}


void main_construction_repr_annot() {
    test_repr_rank();
    test_repr_annotation();
}
