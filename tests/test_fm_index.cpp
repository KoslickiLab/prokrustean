#ifndef TEST_FMINDEX_HPP_
#define TEST_FMINDEX_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <tuple>
#include "const.cpp"	
#include "naive_impl.cpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/fm_index/index.hpp"

using namespace std;


std::ifstream::pos_type filesize(string filename){
    std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}
/*
* check that the string contains exactly the same characters as the file in path
*/
bool check_content(string path){
    auto str = WaveletString(path, '#');
    ifstream ifs(path);
    bool res = true;
    for(uint64_t i=0;i<str.size();++i){
        char c;
        ifs.read((char*)&c, sizeof(char));
        if(str.operator[](i) != c) {
            res = false;
            break;
        }
    }
    if(res){
        cout << "string content is valid" << endl;
    }else{
        cout << "string content is not valid" << endl;
    }
    return res;
}

void test_strings(){
    IS_TRUE(check_content(PATH1_BWT));
}


void test_recovery(){
    // naive
    vector<string> sequences = get_sequences_naive(PATH1_SEQ);
    // fm_index
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    auto fm_idx_naive = NaiveFmIndex(sequences);
    for (int i=0; i<sequences.size(); i++){
        assert(sequences[i]==fm_idx.recover_text(i, true));
    }
    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i, true)){
            sa[pair.first]=make_pair(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

void test_recovery_unsorted(){
    // naive
    vector<string> sequences = get_sequences_naive(PATH2_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences, 3);
    auto str = WaveletString(PATH2_BWT);
    auto fm_idx = FmIndex(str);
    for (int i=0; i<sequences.size(); i++){
        assert(sequences[i]==fm_idx.recover_text(i, true));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i, true)){
            sa[pair.first]=make_pair(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

void test_recovery_unsorted_tied(){
    // naive
    vector<string> sequences = get_sequences_naive(PATH3_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences);
    auto str = WaveletString(PATH3_BWT);
    auto fm_idx = FmIndex(str);
    for (int i=0; i<sequences.size(); i++){
        assert(sequences[i]==fm_idx.recover_text(i, true));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i, true)){
            sa[pair.first]=make_pair(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

void main_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_recovery();
    test_recovery_unsorted();
    test_recovery_unsorted_tied();
}

#endif