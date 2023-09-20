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
#include "../src/construction/algorithms.procedures_new.hpp"
#include "../src/fm_index/tree_new.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

void test_comparison(){
    int Lmin = 1;
    // auto str = WaveletString(PATH1_PERFORMANCE_SREAD_SEQ, '$');
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;

    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats;
    SuffixArrayNode_NEW root_new = get_root_new(fm_idx);
    vector<MaximalRepeatAnnotation> repeats_new;
    
    start = std::chrono::steady_clock::now();
    navigate_tree_new<MaximalRepeatAnnotation, get_repeat_annotations_new>(root_new, Lmin, fm_idx, repeats_new);
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    cout << "original algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    assert(repeats.size()==repeats_new.size());
    cout << "repeats are same: " << repeats.size() << endl;
}

void main_performance_new_tree() {
    test_comparison();
}
