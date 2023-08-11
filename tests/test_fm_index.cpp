#ifndef TEST_FMINDEX_HPP_
#define TEST_FMINDEX_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/rank.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/tree.hpp"

using namespace std;

bool check_rank(SuccintString str){
    ParallelRank p = {};
    bool res = true;

    for(int i=0;i<str.size();++i){
        auto r = str.parallel_rank(i);
        if(p != r){
            res = false;
        }
        p.A += (str.operator[](i)=='A');
        p.C += (str.operator[](i)=='C');
        p.G += (str.operator[](i)=='G');
        p.T += (str.operator[](i)=='T');
    }
    auto r = str.parallel_rank(str.size());
    if(p != r){
        res = false;
    }
    if(res){
        cout << "rank is correct" << endl;
    }else{
        cout << "rank is not correct" << endl;
    }
    return res;

}

std::ifstream::pos_type filesize(string filename){
    std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}
/*
* check that the string contains exactly the same characters as the file in path
*/
bool check_content(string path){
    auto str = SuccintString(path);
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

void test_ranks(){
    auto str = SuccintString(PATH1_BWT);
    IS_TRUE(check_rank(str));
}

void test_recovery(){
    // naive
    vector<string> sequences = get_sequences(PATH1_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences, 3);
    // fm_index
    auto str = SuccintString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==recover_text(fm_idx, i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: recover_suffix_array(fm_idx, i)){
            sa[pair.first]=make_tuple(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        // cout << sa[i].first << ", " << fm_idx_naive.sa[i].first << endl;
        // cout << sa[i].second << ", " << fm_idx_naive.sa[i].second << endl;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

void test_recovery_unsorted(){
    // naive
    vector<string> sequences = get_sequences(PATH2_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences, 3);
    // fm_idx_naive.print_ebwt();
    // fm_index
    auto str = SuccintString(PATH2_BWT);
    auto fm_idx = FmIndex(str);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==recover_text(fm_idx, i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: recover_suffix_array(fm_idx, i)){
            sa[pair.first]=make_tuple(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        // cout << sa[i].first << ", " << fm_idx_naive.sa[i].first << endl;
        // cout << sa[i].second << ", " << fm_idx_naive.sa[i].second << endl;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

void test_recovery_unsorted_tied(){
    // naive
    vector<string> sequences = get_sequences(PATH3_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences);
    // fm_idx_naive.print_ebwt();
    // fm_index
    auto str = SuccintString(PATH3_BWT);
    auto fm_idx = FmIndex(str);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==recover_text(fm_idx, i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: recover_suffix_array(fm_idx, i)){
            sa[pair.first]=make_tuple(text_pos, pair.second);
            text_pos++;
        }
    }
    for (int i=0; i<fm_idx.size(); i++){
        //important: skip terminator only.
        if (i < fm_idx.seq_cnt()) continue;
        // cout << sa[i].first << ", " << fm_idx_naive.sa[i].first << endl;
        // cout << sa[i].second << ", " << fm_idx_naive.sa[i].second << endl;
        assert(sa[i].first==fm_idx_naive.sa[i].first);
        assert(sa[i].second==fm_idx_naive.sa[i].second);
    }
}

tuple<uint64_t, uint64_t> get_sa_range_by_weiner_link(FmIndex &fm_index, string W){
    SuffixArrayInterval interval = get_root(fm_index);
    for(int i=W.size()-1; i>=0; i--){
        CharId c = fm_index.convert_char(W[i]);
        SuffixArrayNode left_ext = to_node(fm_index, interval);
        interval = left_ext.c_intervals[c];
    }
    return make_tuple(interval.firsts[0], interval.firsts[interval.firsts.size()-1]);
}

tuple<uint64_t, uint64_t> get_sa_range(FmIndex &fm_index, string W){
    vector<string> suffixes = recover_suffix_array(fm_index);
    if(W.size()==0){
        return make_tuple(0, suffixes.size());
    }
    vector<uint64_t> list;
    uint64_t i = 0;
    for(auto suffix: recover_suffix_array(fm_index)){
        if(suffix.compare(0, W.size(), W) == 0){
            list.push_back(i);
        }
        i++;
    }
    assert(list.size()!=0);
    return make_tuple(list[0], list[list.size()-1]+1);
}

void test_left_extension(){
    auto str = SuccintString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    // vector<string> suffixes = recover_suffix_array(fm_idx);
    // for(auto suffix: recover_suffix_array(fm_idx)){
    //     cout<< suffix << endl;
    // }
    // test root
    string W = "";
    auto sa_range_wl = get_sa_range_by_weiner_link(fm_idx, W);
    auto sa_range = get_sa_range(fm_idx, W);
    assert(sa_range_wl == sa_range);
    // test single
    W = "G";
    sa_range_wl = get_sa_range_by_weiner_link(fm_idx, W);
    sa_range = get_sa_range(fm_idx, W);
    assert(sa_range_wl == sa_range);
    // test long
    for(auto seq: recover_text(fm_idx)){
        //remove terminator
        W = seq.substr(0, seq.size()-1);
        sa_range_wl = get_sa_range_by_weiner_link(fm_idx, W);
        sa_range = get_sa_range(fm_idx, W);
        // cout << get<0>(sa_range) << ", " << get<1>(sa_range) << endl;
        assert(sa_range_wl == sa_range);
    }
}

void test_left_extension_exhaustive(){
    vector<string> paths = {PATH1_BWT, PATH2_BWT, PATH3_BWT};
    for(auto path: paths){
        auto str = SuccintString(path);
        auto fm_idx = FmIndex(str);
        for(auto seq: recover_text(fm_idx)){
            //remove terminator
            for(int i=0; i<seq.size()-1; i++){
                for(int j=i; j< seq.size()-1; j++){
                    auto W = seq.substr(i, j+1);
                    auto sa_range_wl = get_sa_range_by_weiner_link(fm_idx, W);
                    auto sa_range = get_sa_range(fm_idx, W);
                    assert(sa_range_wl == sa_range);
                }
            }    
            // cout << get<0>(sa_range) << ", " << get<1>(sa_range) << endl;
        }
    }
    
}

void main_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_ranks();
    test_recovery();
    test_recovery_unsorted();
    test_recovery_unsorted_tied();
    test_left_extension();
    test_left_extension_exhaustive();
}

#endif