#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <set>
#include "const.cpp"	
#include "naive_impl.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.support.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/dbg.count.hpp"

using namespace std;
using namespace sdsl;

void test_unitig_counting_range(){
    int Lmin = 15;
    auto num_threads=12;
    WaveletString str(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanExtension ext(prokrustean);
    ext.collect_left_right_extensions=true;
    
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, Lmin, &ext);
    
    auto start = std::chrono::steady_clock::now();
    
    vector<uint64_t> output;
    count_maximal_unitigs_range_of_k(30, 200, ext, output);
    
    cout << "parallel: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}


void main_application_unitig_count() {
    test_unitig_counting_range();
}
