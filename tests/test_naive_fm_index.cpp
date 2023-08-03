#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	

using namespace std;


void test_naive_sa(){
    vector<string> sequences = {"ACACACGTGTGT", "AGGCCCTT"};
    auto fm_idx = NaiveFmIndex(sequences, 3);
    // fm_idx.print_sa();
    // fm_idx.print_ssa();
    // fm_idx.print_ebwt();
}


void main_naive_fm_index() {
    // Call all tests. Using a test framework would simplify this.
    test_naive_sa();
}
