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
    // auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto start = std::chrono::steady_clock::now();
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << " wv tree took " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats;
    SuffixArrayNode_NEW root_new = get_root_new(fm_idx);
    vector<MaximalRepeatAnnotation> repeats_new;
    
    start = std::chrono::steady_clock::now();
    // navigate_tree_new<MaximalRepeatAnnotation, get_repeat_annotations_new>(root_new, Lmin, fm_idx, repeats_new);
    navigate_tree_new<MaximalRepeatAnnotation, report_repr_sa>(root_new, Lmin, fm_idx, repeats_new);
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    cout << "original algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    assert(repeats.size()==repeats_new.size());
    cout << "repeats are same: " << repeats.size() << endl;
    
    for(int i=0; i<repeats.size(); i++){
        assert(repeats[i].repr_indexes.size()==repeats_new[i].repr_indexes.size());
    }
    cout << "repeats reprs are same"  << endl;
}

void test_memory_reserved_stack(){
    std::vector<int> myVector;
    myVector.reserve(100);
    
    std::stack<int, std::vector<int>> myStack(myVector);
    cout << "stack size: "<< myStack.size() << endl;
    // Now push and pop without incurring dynamic allocations up to 100 elements
    myStack.push(1);
    cout << "stack size: "<< myStack.size() << endl;
    myStack.push(2);
    cout << "stack size: "<< myStack.size() << endl;
    myStack.pop();
    cout << "stack size: "<< myStack.size() << endl;
}

void main_performance_new_tree() {
    test_comparison();
    // test_memory_reserved_stack();
}
