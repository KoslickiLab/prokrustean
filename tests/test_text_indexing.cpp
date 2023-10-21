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

void test_metadata_store(){
    Prokrustean prokrustean;
    prokrustean.file_name="temp.something";
    prokrustean.lmin=20;
    prokrustean.sequence_count=123;
    prokrustean.stratum_count=987;
    prokrustean.sequences__size.resize(prokrustean.sequence_count);
    prokrustean.sequences__size[0]=1000;
    prokrustean.sequences__size[1]=2000;
    prokrustean.sequences__size[2]=3000;
    prokrustean.sequences__size[10]=999;

    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(prokrustean);
    sequnce_access.write_close();
    sequnce_access.load_metadata();
    
    // cout << std::string(sequnce_access.metadata.evidence) << endl;
    assert(sequnce_access.metadata.evidence==PROKRUSTEAN_EVIDENCE);
    prokrustean.file_name.resize(256);
    assert(sequnce_access.metadata.prokrustean_file_name==prokrustean.file_name);
    assert(sequnce_access.metadata.lmin==prokrustean.lmin);
    // cout << "sequnce_access.metadata.sequence_count " << sequnce_access.metadata.sequence_count << endl;
    assert(sequnce_access.metadata.sequence_count==prokrustean.sequence_count);
    assert(sequnce_access.metadata.strata_count==prokrustean.stratum_count);
    assert(sequnce_access.sequence_start_positions.size()==prokrustean.sequence_count);
    for(int i=1; i<prokrustean.sequence_count; i++){
        // sequence content = size + letters
        assert(prokrustean.sequences__size[i-1]+sizeof(SequenceSize) == sequnce_access.sequence_start_positions[i]-sequnce_access.sequence_start_positions[i-1]);
    }
}

void test_sequence_save_and_load(){
    std::vector<std::string> sequences = {"apple", "orange", "banana", "grape"};
    Prokrustean prokrustean;
    prokrustean.file_name="temp.something";
    prokrustean.lmin=20;
    prokrustean.sequence_count=123;
    prokrustean.stratum_count=987;
    prokrustean.sequences__size.resize(prokrustean.sequence_count);
    prokrustean.sequences__size[0]=sequences[0].size();
    prokrustean.sequences__size[1]=sequences[1].size();
    prokrustean.sequences__size[2]=sequences[2].size();
    prokrustean.sequences__size[10]=sequences[3].size();

    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(prokrustean);
    sequnce_access.write_strings(sequences);
    sequnce_access.write_close();
    sequnce_access.read_open();
    string str;
    sequnce_access.read_seq(0, str);
    assert(str==sequences[0]);
    sequnce_access.read_seq_substr(2, 1, 3, str);
    assert(str==sequences[2].substr(1,3));
    sequnce_access.read_seq_substr(3, 3, 2, str);
    assert(str==sequences[3].substr(3,2));
}

// void test_indexing_simple(){
//     std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
//     vector<int> seq_sizes;
//     for(auto &str: strings){
//         seq_sizes.push_back(str.size());
//     }
//     DiskSequenceAccess sequnce_access("data.dat");
//     sequnce_access.write_open();
//     sequnce_access.write_meta(seq_sizes);
//     sequnce_access.write_strings(strings);
//     sequnce_access.write_close();
//     sequnce_access.load_all_strings();
//     // sequnce_access.read_open();
//     cout << sequnce_access.get_metadata() << endl;

//     cout << sequnce_access.get_string(1) << endl;
//     assert("orange"==sequnce_access.get_string(1));
//     assert("rang"==sequnce_access.get_substring(1, 1, 4));
//     sequnce_access.read_close();
// }

void test_bwt_prokrustean_indexing(){
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42)
    };

    DiskSequenceAccess sequnce_access("test_bwt_prokrustean_indexing.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(prokrustean);
    sequnce_access.write_strings(seq_texts);
    sequnce_access.write_close();
    
    DiskSequenceAccess sequnce_access2("test_bwt_prokrustean_indexing.dat");
    sequnce_access2.load_metadata();
    sequnce_access2.read_open();
    string result;
    for(auto [idx, from, to]: substring_locs){
        sequnce_access2.read_seq(idx, result);
        assert(seq_texts[idx]==result);
        sequnce_access2.read_seq_substr(idx, from, to-from+1, result);
        assert(seq_texts[idx].substr(from, to-from+1)==result);
    }
    sequnce_access2.read_close(); 
    //in memory
    // sequnce_access.load_all_strings();
    // for(auto [idx, from, to]: substring_locs){
    //     assert(seq_texts[idx]==sequnce_access.get_string(idx));
    //     assert(seq_texts[idx].substr(from, to-from+1)==sequnce_access.get_substring(idx, from, to-from+1));
    // }
}

// void test_verification(){
//     std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
//     DiskSequenceAccess sequnce_access("data.dat");
//     vector<int> seq_sizes;
//     for(auto &str: strings){
//         seq_sizes.push_back(str.size());
//     }
//     sequnce_access.write_open();
//     sequnce_access.write_meta(seq_sizes);
//     sequnce_access.write_strings(strings);
//     sequnce_access.write_close();
//     DiskSequenceAccess sequnce_access1("data.dat");
//     sequnce_access1.verify();
// }


// void test_bwt_indexing2(){
//     WaveletString str(PATH6_CDBG_SAMPLE2, '$');
//     auto fm_idx = FmIndex(str);
    
//     vector<string> seq_texts;
//     fm_idx.recover_all_texts(seq_texts);
//     // seq_texts={ "orange", "apple"};
//     vector<tuple<int, int, int>> substring_locs={
//         make_tuple(1, 3, 20),
//         make_tuple(3, 8, 37),
//         make_tuple(9, 20, 42)
//     };
//     // vector<tuple<int, int, int>> substring_locs={
//     //     make_tuple(1, 1, 3),
//     //     make_tuple(0, 1, 3)
//     // };
//     Access sequnce_access;
//     sequnce_access.filename="test_data.dat";
//     sequnce_access.outfile=std::ofstream(sequnce_access.filename, std::ios::binary);
//     vector<int> sizes;
//     for(auto &str: seq_texts){
//         sizes.push_back(str.size());
//     }
//     sequnce_access.save_meta(sizes);
//     sequnce_access.save_strings(seq_texts);
//     for(auto [idx, from, size]: substring_locs){
//         cout << sequnce_access.get_string(idx) << endl;
//         assert(seq_texts[idx]==sequnce_access.get_string(idx));
//         assert(seq_texts[idx].substr(from, size)==sequnce_access.get_substring(idx, from, size));
//     }
// }


void main_text_indexing() {
    test_metadata_store();
    test_sequence_save_and_load();
    // test_indexing_simple();
    test_bwt_prokrustean_indexing();
    // test_verification();
}

#endif