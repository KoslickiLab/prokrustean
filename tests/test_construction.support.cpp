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

// void test_left_right_extension_counting(){
//     int Lmin = 1;
//     WaveletString str(PATH6_CDBG_SAMPLE2, '$');
//     auto fm_idx = FmIndex(str);
    
//     Prokrustean prokrustean;
//     ProkrusteanExtension ext(prokrustean);
//     ext.collect_left_right_extensions=true;
//     construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin, &ext);
//     vector<CharCount> stratum_left_ext_count_original=ext.stratum_left_ext_count;
//     vector<CharCount> stratum_right_ext_count_original=ext.stratum_right_ext_count;

//     count_left_right_character_extensions(ext);
    
//     // assert(stratum_left_ext_count_original==ext.stratum_left_ext_count);
//     // assert(stratum_right_ext_count_original==ext.stratum_right_ext_count);
//     for(int i=0; prokrustean.stratum_count; i++){
//         cout << "stratum_left_ext_count_original[i] " << (int)stratum_left_ext_count_original[i] << " ext.stratum_left_ext_count " << (int)ext.prokrustean.get_left_cnt(i)<< endl;
//         assert(stratum_left_ext_count_original[i]==ext.prokrustean.get_left_cnt(i));
//     // assert(stratum_right_ext_count_original==ext.stratum_right_ext_count);
//     }
// }

void test_left_right_storage(){
    Prokrustean prokrustean;
    prokrustean.stratums__character_extensions_cnt.resize(1);
    CharCount left=1;
    CharCount right=5;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));

    left=0;
    right=0;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));
    
    left=3;
    right=1;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));

    left=4;
    right=4;
    prokrustean.set_stratum_left_right_extensions(0, left, right);
    assert(left==prokrustean.get_left_cnt(0));
    assert(right==prokrustean.get_right_cnt(0));
}
void main_prokrustean_support(){
    // test_left_right_extension_counting();
    test_left_right_storage();
}
