#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/rank.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"


using namespace std;

bool check_rank(succint_string str){
    p_rank p = {};
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
    auto str = succint_string(path);
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

void print_interval(interval in){
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
string PATH2_SEQ = "../data/2_sequences.txt";
string PATH2_BWT = "../data/2_ebwt.txt";

void test_strings(){
    IS_TRUE(check_content(PATH1_BWT));
}

void test_ranks(){
    auto str = succint_string(PATH1_BWT);
    IS_TRUE(check_rank(str));
}

void test_LF(){
    auto idx = fm_index(PATH1_BWT);
    interval root = idx.root();
    left_ext_intervals left_exts = idx.LF(root);
    // print_left_ext_intervals(left_exts);
    // print_left_ext_intervals(idx.LF(left_exts.TERM));
    // print_left_ext_intervals(idx.LF(left_exts.A));
    // print_left_ext_intervals(idx.LF(left_exts.C));
    // print_left_ext_intervals(idx.LF(left_exts.G));
    // print_left_ext_intervals(idx.LF(left_exts.T));
}

void test_sampled_suffixes(){
    // naive
    vector<string> sequences = get_sequences(PATH1_SEQ);
    auto fm_idx_naive = NaiveFmIndex(sequences, 3);
    fm_idx_naive.print_sa();
    fm_idx_naive.print_ebwt();
    // fm_index
    auto fm_idx = fm_index(PATH1_BWT);
    sample_by_sa_order(fm_idx, 3);
    for (int i=0; i<fm_idx.size(); i++){
        cout << "sa comparison " << fm_idx_naive.get_suffix(i) << ":" << fm_idx.get_suffix(i) << endl;
        assert(fm_idx_naive.get_suffix(i) == fm_idx.get_suffix(i));
    }
}

void main_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_ranks();
    test_LF();
    test_sampled_suffixes();
}
