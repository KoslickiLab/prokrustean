#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/rank.hpp"
#include "../src/fm_index/index.hpp"

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

void print_interval(Interval in){
    cout << in.first_TERM << " , " << in.first_A << " , " << in.first_C << " , " << in.first_G << " , " << in.first_T << " , " << in.last << endl;
}
void print_left_ext_intervals(left_ext_intervals left_exts){
    cout << "left ext TERM: " << ends;
    print_interval(left_exts.TERM);
    cout << "left ext A: " << ends;
    print_interval(left_exts.A);
    cout << "left ext C: " << ends;
    print_interval(left_exts.C);
    cout << "left ext G: " << ends;
    print_interval(left_exts.G);
    cout << "left ext T: " << ends;
    print_interval(left_exts.T);
}

string PATH1_SEQ = "../data/1_sequences.txt";
string PATH1_BWT = "../data/1_ebwt.txt";
string PATH2_SEQ = "../data/2_sequences_unsorted.txt";
string PATH2_BWT = "../data/2_ebwt.txt";
string PATH3_SEQ = "../data/3_sequences_unsorted_tied.txt";
string PATH3_BWT = "../data/3_ebwt.txt";

void test_strings(){
    IS_TRUE(check_content(PATH1_BWT));
}

void test_ranks(){
    auto str = SuccintString(PATH1_BWT);
    IS_TRUE(check_rank(str));
}

void test_LF(){
    auto idx = FmIndex(PATH1_BWT);
    Interval root = idx.root();
    left_ext_intervals left_exts = idx.LF(root);
    // print_left_ext_intervals(left_exts);
    // print_left_ext_intervals(idx.LF(left_exts.TERM));
    // print_left_ext_intervals(idx.LF(left_exts.A));
    // print_left_ext_intervals(idx.LF(left_exts.C));
    // print_left_ext_intervals(idx.LF(left_exts.G));
    // print_left_ext_intervals(idx.LF(left_exts.T));
}

void test_recovery(){
    // naive
    vector<string> sequences = get_sequences(PATH1_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences, 3);
    // fm_index
    auto fm_idx = FmIndex(PATH1_BWT);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==fm_idx.recover_text(i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i)){
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
    auto fm_idx = FmIndex(PATH2_BWT);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==fm_idx.recover_text(i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i)){
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
    auto fm_idx = FmIndex(PATH3_BWT);
    for (int i=0; i<sequences.size(); i++){
        // cout << recover_text(fm_idx, i) <<endl;
        // cout << sequences[i] <<endl;
        assert(sequences[i]==fm_idx.recover_text(i));
    }

    vector<pair<uint64_t, string>> sa(fm_idx.size());
    uint64_t text_pos = 0;
    for (int i=0; i<sequences.size(); i++){
        for(auto pair: fm_idx.recover_suffix_array(i)){
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

void main_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_ranks();
    test_LF();
    test_recovery();
    test_recovery_unsorted();
    test_recovery_unsorted_tied();
}
