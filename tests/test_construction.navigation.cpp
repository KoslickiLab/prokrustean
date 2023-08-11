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

void _compare_lcp_recovered_vs_tree(FmIndex &fm_idx){
    auto sa = recover_suffix_array(fm_idx);
    auto lcp = recover_lcp(fm_idx);
    
    vector<string> lcp_str;
    for(int i; i< lcp.size(); i++){
        string str_i = sa[i].substr(0, lcp[i]);
        if(str_i.size()==0 || str_i=="#"){
            continue;
        }
        lcp_str.push_back(str_i);
    }
    sort(lcp_str.begin(), lcp_str.end());
    lcp_str.erase(unique(lcp_str.begin(), lcp_str.end() ), lcp_str.end());
    
    SuffixArrayInterval root = get_root(fm_idx);
    auto pairs = navigate_tree<tuple<uint64_t, uint64_t>, _get_lcp>(root, 0, fm_idx);
    vector<string> tree_lcp_str;
    for(auto p: pairs){
        uint64_t depth = get<0>(p);
        uint64_t sa_idx = get<1>(p);
        tree_lcp_str.push_back(sa[sa_idx].substr(0, depth));
    }
    sort(tree_lcp_str.begin(), tree_lcp_str.end());

    for(int i; i<tree_lcp_str.size(); i++){
        // cout << lcp_str[i] << " , " << tree_lcp_str[i] << endl;
        assert(lcp_str[i]==tree_lcp_str[i]);
    }
}

void test_lcp_equivalence() {
    // fm_index
    auto str = SuccintString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    _compare_lcp_recovered_vs_tree(fm_idx);

    str = SuccintString(PATH2_BWT);
    fm_idx = FmIndex(str);
    _compare_lcp_recovered_vs_tree(fm_idx);

    str = SuccintString(PATH3_BWT);
    fm_idx = FmIndex(str);
    _compare_lcp_recovered_vs_tree(fm_idx);
}

// template<typename T>
// using NodeProcess = optional<T> (*)(NodeLeftExtension&, FmIndex&);

void main_construction_navigation() {
    test_lcp_equivalence();
}

#endif