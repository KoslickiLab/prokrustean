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

    vector<CharId> distinct_extensions(){
        vector<CharId> exts = {};
        for(int i=0; i<firsts.size()-1; i++){
            if(firsts[i]<firsts[i+1]) 
                exts.push_back(i);
        }
        return exts;
    }

    uint64_t interval_size(){
        return firsts[firsts.size()-1]-firsts[0];
    }

    bool right_maximal(){
        return distinct_extensions().size()>1 || /*suffix count*/ firsts[1]-firsts[0]>1;
    }

    SuffixArrayIdx get_valid_first(){
        for(int i=0; i<firsts.size()-1; i++){
            if(firsts[i]<firsts[i+1]) 
                return firsts[i];
        }
        assert(false); //interval might be 0 - invalid
    }
};

struct SuffixArrayNodeExtension {
    SuffixArrayNode node;
    // left extended intervals. Each means sa intervals of cW for each character c.
    vector<SuffixArrayNode> c_nodes;
    // Each means the first sa index of W that corresponds to first sa index of cW for each character c.
    vector<optional<SuffixArrayIdx>> prev_c_firsts;

    vector<tuple<CharId, CharId>> distinct_extensions(){
        vector<tuple<CharId, CharId>> exts = {};
        for(int i=0; i<c_nodes.size(); i++){
            for(auto r: c_nodes[i].distinct_extensions()){
                exts.push_back(make_tuple(i, r));    
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

     SuffixArrayIdx first_l(CharId c){
        assert(prev_c_firsts[c].has_value());
        // cout << "first l: " << prev_c_firsts[c].value() << endl;
        return prev_c_firsts[c].value();
    }

    SuffixArrayIdx first_r(CharId c){
        return node.firsts[c];
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
    for(int c=0; c < index.STRING->get_characters().size(); c++){
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = index.C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        c_nodes.push_back({firsts, node.depth+1});
    }

    vector<optional<SuffixArrayIdx>> prev_c_firsts;
    for(int c=0; c < c_nodes.size(); c++){
        SuffixArrayNode c_node = c_nodes[c];

        if(c_node.interval_size()==0){
            prev_c_firsts.push_back(nullopt);
            continue;
        }

        //revert the rank process
        SuffixArrayIdx sa_idx = c_node.get_valid_first();
        uint64_t rank = sa_idx - index.C[c] + 1;
        SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
        // cout << "obtained " << prev_sa_idx << " for " << sa_idx << endl;
        prev_c_firsts.push_back(prev_sa_idx);
    }
    
    return {node, c_nodes, prev_c_firsts};
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

template<class T> using NodeFunc = optional<T>(*)(SuffixArrayNodeExtension&);

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
            optional<T> t = process_node(ext);
            if(t.has_value()) {
                Ts.push_back(t.value());
            }
        }

        for(int i=0; i<ext.c_nodes.size(); i++){
            // terminal
            if(i==0) continue;

            auto child = ext.c_nodes[i];
            if(child.interval_size()>1){
                stack.push(child);
            }
        }
    }
    return Ts;
}

#endif /* FM_INDEX_TREE_HPP_ */