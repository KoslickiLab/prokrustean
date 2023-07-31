#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/rank.hpp"
#include "../src/fm_index/index.hpp"


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

string DATA_PATH1 = "../data/simple_ebwt.txt";

void test_strings(){
    IS_TRUE(check_content(DATA_PATH1));
}

void test_ranks(){
    auto str = succint_string(DATA_PATH1);
    IS_TRUE(check_rank(str));
}

void test_LF(){
    auto str = succint_string(DATA_PATH1);
    auto idx = fm_index(DATA_PATH1);
    interval root = idx.root();
    left_ext_intervals left_exts = idx.LF(root);
    print_left_ext_intervals(left_exts);
}

void main_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_ranks();
    test_LF();
}
