#include "index.hpp"
#include <algorithm>
#include <limits>
#include <stack>

using namespace std;

#ifndef FM_INDEX_TREE_HPP_
#define FM_INDEX_TREE_HPP_

/*
 * representation of a right-maximal substring (SA node) as a list of STRING intervals
 */
struct SuffixArrayNode {
    // |characters| + 1 sa positions. The last element is not corresponding to any character
    // ex firsts[0] => sa location of starting terminator
    vector<uint64_t> firsts;
    uint64_t depth;

    vector<tuple<CharId, uint64_t>> distinct_extensions(){
        vector<tuple<CharId, uint64_t>> exts = {};
        for(int i=0; i<firsts.size()-1; i++){
            if(firsts[i]<firsts[i+1]) 
                exts.push_back(make_tuple(i, firsts[i]));
        }
        return exts;
    }

    uint64_t interval_size(){
        return firsts[firsts.size()-1]-firsts[0];
    }

    bool right_maximal(){
        return distinct_extensions().size()>1 || /*suffix count*/ firsts[1]-firsts[0]>1;
    }
};

struct SuffixArrayNodeExtension {
    SuffixArrayNode node;
    // left extended intervals. Each means sa intervals of cW for each character c.
    vector<SuffixArrayNode> c_nodes;
    vector<CharId> c_s;

    vector<tuple<CharId, CharId, uint64_t>> distinct_extensions(){
        vector<tuple<CharId, CharId, uint64_t>> exts = {};
        for(int i=0; i<c_nodes.size(); i++){
            for(auto r: c_nodes[i].distinct_extensions()){
                exts.push_back(make_tuple(i, get<0>(r), get<1>(r)));    
            }
        }
        return exts;
    }

    bool left_maximal(){
        int distinct = 0;
        for(auto interval:c_nodes){
            if(interval.interval_size()>0){
                distinct++;
            }
        }
        return distinct>1 || /*prefix count*/ c_nodes[0].interval_size()>1;
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
SuffixArrayNodeExtension extend_node(FmIndex &index, SuffixArrayNode &node){
    vector<RankArray> left_p_ranks;
    RankArray rank_array;
    for(int i=0; i<node.firsts.size(); i++){
        if(i==0 || node.firsts[i-1]!=node.firsts[i]){
            rank_array = index.STRING->ranks(node.firsts[i]);  
        } 
        left_p_ranks.push_back(rank_array);
    }

    vector<SuffixArrayNode> c_nodes;
    vector<CharId> c_s;
    for(int c=0; c < index.STRING->get_characters().size(); c++){
        uint64_t C = index.C[c];
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = index.C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        c_nodes.push_back({firsts, node.depth+1});
        c_s.push_back(c);
    }
    
    return {node, c_nodes, c_s};
};

/*
* functions for suffix tree navigation
*/
SuffixArrayNode get_root(FmIndex &index){
    vector<uint64_t> firsts;
    for(auto rank: index.C){
        firsts.push_back(rank);
    }
    firsts.push_back(index.STRING->size());
    return {firsts, 0};
};

template<class T> using NodeFunc = optional<T>(*)(SuffixArrayNodeExtension&, FmIndex&);

template<class T, NodeFunc<T> process_node>
vector<T> navigate_tree(SuffixArrayNode &root, int Lmin, FmIndex &fm_idx){
    Lmin = Lmin >= 1? Lmin : 1;

    std::stack<SuffixArrayNode> stack;
    vector<T> Ts;
    
    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();

        auto ext = extend_node(fm_idx, node);
        if(ext.node.depth>=Lmin){
            optional<T> t = process_node(ext, fm_idx);
            if(t.has_value()) {
                Ts.push_back(t.value());
            }
        }

        for(int i=0; i<ext.c_nodes.size(); i++){
            // terminal
            if(ext.c_s[i]==0) 
            continue;
            auto child = ext.c_nodes[i];
            if(child.interval_size()>1){
                stack.push(child);
            }
        }
    }
    return Ts;
}

#endif /* FM_INDEX_TREE_HPP_ */