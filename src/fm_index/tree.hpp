#ifndef FM_INDEX_TREE_HPP_
#define FM_INDEX_TREE_HPP_

#include "index.hpp"
#include <algorithm>
#include <limits>
#include <stack>
#include <tuple>

using namespace std;
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
        int cnt=0;
        for(int i=0; i<firsts.size()-1; i++){
            if(firsts[i]<firsts[i+1]) cnt++;
        }
        return cnt>1 || /*suffix count*/ firsts[1]-firsts[0]>1;
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
    // suffixes such that both left and right extensions are terminations
    vector<SuffixArrayIdx> both_ext_terms;

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


struct WorkSpace{
    vector<RankArray> left_p_ranks;
};
/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
void extend_node(FmIndex &index, SuffixArrayNode &node, SuffixArrayNodeExtension &ext){
    int p_rank_size = 1+index.characters.size();
    vector<RankArray> left_p_ranks(p_rank_size);
    RankArray rank_array;
    for(int i=0; i<node.firsts.size(); i++){
        if(i==0 || node.firsts[i-1]!=node.firsts[i]){
            rank_array = index.STRING->ranks(node.firsts[i]);  
        } 
        left_p_ranks[i]=rank_array;
    }

    // vector<SuffixArrayNode> c_nodes(index.characters.size());
    for(int c=0; c < index.character_cnt; c++){
        vector<uint64_t> firsts(p_rank_size);
        int i=0;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = index.C[c]+p_rank[c];
            // firsts[i]=sa_idx;
            ext.c_nodes[c].firsts[i]=sa_idx;
            i++;
        }
        // ext.c_nodes[c]={firsts, node.depth+1};
        ext.c_nodes[c].depth=node.depth+1;
    }

    // vector<optional<SuffixArrayIdx>> prev_c_firsts(ext.c_nodes.size());
    for(int c=0; c < ext.c_nodes.size(); c++){
        SuffixArrayNode c_node = ext.c_nodes[c];

        if(c_node.interval_size()==0){
            ext.prev_c_firsts[c]=nullopt;
            continue;
        }

        //revert the rank process
        SuffixArrayIdx sa_idx = c_node.get_valid_first();
        uint64_t rank = sa_idx - index.C[c] + 1;
        SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
        // cout << "obtained " << prev_sa_idx << " for " << sa_idx << endl;
        ext.prev_c_firsts[c]=prev_sa_idx;
    }

    // vector<SuffixArrayIdx> both_ext_terms;
    ext.both_ext_terms.clear();
    if(ext.c_nodes[0].firsts[0]!=ext.c_nodes[0].firsts[1]){
        // both ends are termination. collect all
        for(SuffixArrayIdx sa_idx=ext.c_nodes[0].firsts[0]; sa_idx< ext.c_nodes[0].firsts[1]; sa_idx++){
            //revert the rank process
            uint64_t rank = sa_idx - index.C[0] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, 0);
            ext.both_ext_terms.push_back(prev_sa_idx);
        }
    }

    ext.node = node;
    // return {node, c_nodes, prev_c_firsts, both_ext_terms};
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

template<class T> using NodeFunc = T(*)(SuffixArrayNodeExtension&);

template<class T, NodeFunc<T> process_node>
void navigate_tree(SuffixArrayNode &root, int Lmin, FmIndex &fm_idx, vector<T> &Ts){
    Lmin = Lmin >= 1? Lmin : 1;

    std::stack<SuffixArrayNode> stack;
    
    int node_cnt=0;
    int process_cnt=0;

    // vector<SuffixArrayNode> c_nodes(fm_idx.characters.size());
    SuffixArrayNodeExtension ext;
    ext.c_nodes.resize(fm_idx.character_cnt);
    for(int i=0; i< fm_idx.characters.size(); i++){
        ext.c_nodes[i].firsts.resize(fm_idx.character_cnt+1);
    }
    ext.prev_c_firsts.resize(fm_idx.character_cnt);
    // int T_idx=0;
    // int min_margin_size=1000;
    // Ts.resize(min_margin_size);
    
    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();

        extend_node(fm_idx, node, ext);
        if(ext.node.depth>=Lmin && ext.left_maximal()){
            T t = process_node(ext);
            Ts.push_back(t);
            // Ts[T_idx]=t;
            // T_idx++;
            // for performance
            // if(T_idx==Ts.size()){
            //     int margin_size = 0.01*Ts.size() > min_margin_size? 0.01*Ts.size() : min_margin_size;
            //     Ts.resize(Ts.size()+margin_size);
            // }
            process_cnt++;
        }

        for(int i=0; i<ext.c_nodes.size(); i++){
            // terminal
            if(i==0) continue;

            auto child = ext.c_nodes[i];
            if(child.right_maximal()){
                stack.push(child);
            }
        }
        node_cnt++;
        // if(node_cnt%10000==0){
        //     cout << "navigated " << node_cnt << endl;
        // }
        // if(process_cnt>100000 && process_cnt%100000==0){
        //     cout << "processed " << process_cnt << endl;
        // }
    }
    // Ts.resize(T_idx+1);
}


vector<SuffixArrayNode> collect_nodes(SuffixArrayNode root, FmIndex &fm_idx, int depth_max){
    std::stack<SuffixArrayNode> stack;
    vector<SuffixArrayNode> nodes;
    // vector<SuffixArrayNode> c_nodes(fm_idx.characters.size());
    SuffixArrayNodeExtension ext;
    ext.c_nodes.resize(fm_idx.characters.size());
    for(int i=0; i< fm_idx.characters.size(); i++){
        ext.c_nodes[i].firsts.resize(fm_idx.character_cnt+1);
    }
    ext.prev_c_firsts.resize(fm_idx.characters.size());

    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();
        if(node.depth>=depth_max){
            nodes.push_back(node);
        } else {
            extend_node(fm_idx, node, ext);
            for(int i=0; i<ext.c_nodes.size(); i++){
                // terminal
                if(i==0) continue;

                auto child = ext.c_nodes[i];
                if(child.right_maximal()){
                    stack.push(child);
                }
            }
        }
    }
    return nodes;
}


#endif /* FM_INDEX_TREE_HPP_ */