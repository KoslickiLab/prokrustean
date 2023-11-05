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
    vector<SequenceSize> sequences__size={10,20,50};
    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(sequences__size);
    sequnce_access.write_close();

    DiskSequenceAccess sequnce_access2("data.dat");
    sequnce_access2.load_metadata();
    
    assert(sequnce_access.metadata.evidence==sequnce_access2.metadata.evidence);
    assert(sequnce_access.metadata.bwt_file_name==sequnce_access2.metadata.bwt_file_name);
    assert(sequnce_access.metadata.sequence_count==sequnce_access2.metadata.sequence_count);
    assert(sequnce_access.metadata.sequence_size_type_byte_size==sequnce_access2.metadata.sequence_size_type_byte_size);
}


void test_sequence_save_and_load(){
    std::vector<std::string> sequences = {"apple", "orange", "banana", "grape"};
    vector<SequenceSize> sequences__size;
    for(auto &seq: sequences){
        sequences__size.push_back(seq.size());
    }
    
    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(sequences__size);
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

void test_bwt_prokrustean_indexing(){
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42),
        make_tuple(seq_sizes.size()-1, 10, 30),
    };

    DiskSequenceAccess sequnce_access("test_bwt_prokrustean_indexing.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(seq_sizes);
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
}

void test_bwt_prokrustean_indexing_loaded(){
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42),
        make_tuple(seq_sizes.size()-1, 10, 30),
    };

    DiskSequenceAccess sequnce_access("test_bwt_prokrustean_indexing.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(seq_sizes);
    sequnce_access.write_strings(seq_texts);
    sequnce_access.write_close();
    
    DiskSequenceAccess sequnce_access2("test_bwt_prokrustean_indexing.dat");
    sequnce_access2.load_metadata();
    sequnce_access2.load_all_strings();
    sequnce_access2.read_open();
    string result;
    for(auto [idx, from, to]: substring_locs){
        sequnce_access2.read_seq(idx, result);
        assert(seq_texts[idx]==result);
        sequnce_access2.read_seq_substr(idx, from, to-from+1, result);
        assert(seq_texts[idx].substr(from, to-from+1)==result);
    }
    sequnce_access2.read_close(); 
}

void test_bwt_prokrustean_indexing_loaded_vs_unloaded(){
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42),
        make_tuple(seq_sizes.size()-1, 10, 30),
    };

    DiskSequenceAccess sequnce_access("test_bwt_prokrustean_indexing.dat");
    sequnce_access.write_open();
    sequnce_access.write_metadata(seq_sizes);
    sequnce_access.write_strings(seq_texts);
    sequnce_access.write_close();
    
    DiskSequenceAccess sequnce_access1("test_bwt_prokrustean_indexing.dat");
    sequnce_access1.load_metadata();
    sequnce_access1.read_open();

    DiskSequenceAccess sequnce_access2("test_bwt_prokrustean_indexing.dat");
    sequnce_access2.load_metadata();
    sequnce_access2.load_all_strings();
    sequnce_access2.read_open();
    string result1;
    string result2;
    for(auto [idx, from, to]: substring_locs){
        sequnce_access1.read_seq(idx, result1);
        sequnce_access2.read_seq(idx, result2);
        assert(result1==result2);
        sequnce_access1.read_seq_substr(idx, from, to-from+1, result1);
        sequnce_access2.read_seq_substr(idx, from, to-from+1, result2);
        assert(result1==result2);
    }
    sequnce_access1.read_close(); 
    sequnce_access2.read_close(); 
}


void main_text_indexing() {
    test_metadata_store();
    test_sequence_save_and_load();
    test_bwt_prokrustean_indexing();
    test_bwt_prokrustean_indexing_loaded();
    test_bwt_prokrustean_indexing_loaded_vs_unloaded();
}

#endif