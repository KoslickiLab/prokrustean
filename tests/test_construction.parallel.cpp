#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

void test_basic_parallel_construction(){
    int Lmin = 2;
    int thread = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    Prokrustean pk = build_prokrustean_parallel(fm_idx, thread, Lmin, true);
    assert(pk.stratums.size()>0);
    assert(pk.seqs.size()>0);
    // for(int i=0; i< pk.stratums.size(); i++){
    //     assert(pk.stratums[i].id == i);
    // }
    // for(int i=0; i< pk.seqs.size(); i++){
    //     assert(pk.seqs[i].id == i);
    // }
}

void main_construction_parallel() {
    test_basic_parallel_construction();
}
