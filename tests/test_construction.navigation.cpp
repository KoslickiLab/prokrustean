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

optional<tuple<uint64_t, uint64_t>> _get_lcp(SuffixArrayNodeExtension &ext){
    if(ext.node.right_maximal()){
        return make_tuple(ext.node.depth, ext.node.firsts[0]);
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
    
    SuffixArrayNode root = get_root(fm_idx);
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

void test_intervals(){
    int Lmin = 2;
    auto str = SuccintString(PATH1_BWT);
    auto fm_idx = FmIndex(str);

    SuffixArrayNode root = get_root(fm_idx);
    SuffixArrayNodeExtension ext = extend_node(fm_idx, root);
    ext = extend_node(fm_idx, ext.c_nodes[2]);
    // ext = extend_node(fm_idx, ext.c_nodes[1]);
    // ext = extend_node(fm_idx, ext.c_nodes[4]);
    cout << "----- interval -----" << endl;
    for(int l=0; l< ext.c_nodes.size(); l++){
        auto node = ext.c_nodes[l];
        cout << "left: ";
        switch (l)
        {
        case 0: cout << "#"; break;
        case 1: cout << "A"; break;
        case 2: cout << "C"; break;
        case 3: cout << "G"; break;
        case 4: cout << "T"; break;
        default: cout << "LAST"; break;
        }
        cout << endl;
        //"AC"
        for(int i=0; i< node.firsts.size(); i++){
            auto sa_idx = node.firsts[i];
            switch (i)
            {
            case 0: cout << "#"; break;
            case 1: cout << "A"; break;
            case 2: cout << "C"; break;
            case 3: cout << "G"; break;
            case 4: cout << "T"; break;
            default: cout << "LAST"; break;
            }
            cout << ": " << sa_idx << ", ";
        }
        cout << endl;   
    }
}

// template<typename T>
// using NodeProcess = optional<T> (*)(NodeLeftExtension&, FmIndex&);

void main_construction_navigation() {
    // test_lcp_equivalence();
    test_intervals();
}

#endif