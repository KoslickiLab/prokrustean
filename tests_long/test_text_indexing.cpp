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
#include "../src/construction/algorithms.hpp"

using namespace std;

void test_bwt_prokrustean_indexing_update(){
    WaveletString str(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx = FmIndex(str);
    Prokrustean prokrustean;
    construct_prokrustean_parallel(fm_idx, prokrustean, 12);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42)
    };
    auto start = std::chrono::steady_clock::now();
    DiskSequenceAccess sequence_access("test_bwt_prokrustean_indexing.dat");
    sequence_access.write_open();
    sequence_access.write_metadata(prokrustean);
    sequence_access.write_close();
    start = std::chrono::steady_clock::now();
    cout << "save meta: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    sequence_access.update_open();
    

    for(int i=0; i< seq_texts.size(); i++){
        sequence_access.update_single_sequence(i, seq_texts[i]);
    }
    cout << "save sequences " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    sequence_access.update_close();

    DiskSequenceAccess sequence_access2("test_bwt_prokrustean_indexing.dat");
    sequence_access2.load_metadata();
    sequence_access2.read_open();
    string result;
    for(auto [idx, from, to]: substring_locs){
        sequence_access2.read_seq(idx, result);
        assert(seq_texts[idx]==result);
        sequence_access2.read_seq_substr(idx, from, to-from+1, result);
        assert(seq_texts[idx].substr(from, to-from+1)==result);
    }
    sequence_access2.read_close(); 
    //in memory
    // sequence_access.load_all_strings();
    // for(auto [idx, from, to]: substring_locs){
    //     assert(seq_texts[idx]==sequence_access.get_string(idx));
    //     assert(seq_texts[idx].substr(from, to-from+1)==sequence_access.get_substring(idx, from, to-from+1));
    // }
}

void main_text_indexing() {
    // test_metadata_store();
    // test_sequence_save_and_load();
    test_bwt_prokrustean_indexing_update();
    // test_bwt_prokrustean_indexing();
    // test_verification();
}

#endif