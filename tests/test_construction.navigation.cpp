#ifndef TEST_CONST_NAVI_HPP_
#define TEST_CONST_NAVI_HPP_
#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/rank.hpp"
#include "../src/fm_index/tree.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/construction/models.hpp"

using namespace std;

optional<tuple<uint64_t, uint64_t>> _get_lcp(SuffixArrayNode &node, FmIndex &fm_idx){
    if(node.interval.right_maximal()){
        return make_tuple(node.interval.depth, node.interval.firsts[0]);
    } else {
        return nullopt;
    }
}

void test_lcp_equivalence() {
    // fm_index
    auto str = SuccintString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    auto sa = recover_suffix_array(fm_idx);
    auto lcp = recover_lcp(fm_idx);
    // for(auto s: sa){
    //     cout << s << endl;
    // }
    for(int i; i< lcp.size(); i++){
        cout << sa[i].substr(0, lcp[i]) << endl;
    }
    
    uint64_t Lmin = 0;
    SuffixArrayInterval root = get_root(fm_idx);
    auto pairs = navigate_tree<tuple<uint64_t, uint64_t>, _get_lcp>(root, Lmin, fm_idx);
    vector<string> lcp_str;
    for(auto p: pairs){
        uint64_t depth = get<0>(p);
        uint64_t sa_idx = get<1>(p);
        cout << "sa: " << sa_idx << " depth: " << depth << endl;
        // cout << sa[sa_idx].substr(0, depth) << endl;
        lcp_str.push_back(sa[sa_idx].substr(0, depth));
    }

    sort(lcp_str.begin(), lcp_str.end());
    for(auto s: lcp_str){
        cout<< s << endl;
    }
    // SuffixArrayInterval interval = get_root(fm_idx);
    // navi<int, _get_int>(interval, fm_idx);
}

// template<typename T>
// using NodeProcess = optional<T> (*)(NodeLeftExtension&, FmIndex&);

void main_construction_navigation() {
    test_lcp_equivalence();
}

#endif