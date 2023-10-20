#ifndef TEST_TEXT_INDEXING_HPP_
#define TEST_TEXT_INDEXING_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>	
#include "const.cpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/util/string.access.hpp"

using namespace std;

void test_indexing_simple(){
    std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.save_strings(strings);
    assert("orange"==sequnce_access.get_string(1));
    assert("rang"==sequnce_access.get_substring(1, 1, 4));
}

void test_bwt_indexing(){
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<tuple<int, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42)
    };
    DiskSequenceAccess sequnce_access("data1.dat");
    sequnce_access.save_strings(seq_texts);
    for(auto [idx, from, to]: substring_locs){
        assert(seq_texts[idx]==sequnce_access.get_string(idx));
        assert(seq_texts[idx].substr(from, to-from+1)==sequnce_access.get_substring(idx, from, to-from+1));
    }
}

void test_verification(){
    std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.save_strings(strings);
    DiskSequenceAccess sequnce_access1("data.dat");
    assert(sequnce_access1.verify());
}

void main_text_indexing() {
    // Call all tests. Using a test framework would simplify this.
    test_indexing_simple();
    test_bwt_indexing();
    test_verification();
}

#endif