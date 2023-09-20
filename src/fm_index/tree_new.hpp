#ifndef FM_INDEX_TREE_NEW_HPP_
#define FM_INDEX_TREE_NEW_HPP_

#include "index.hpp"
#include <algorithm>
#include <limits>
#include <stack>
#include <tuple>

using namespace std;
/*
 * representation of a right-maximal substring (SA node) as a list of STRING intervals
 */
struct SuffixArrayNode_NEW {
    // |characters| + 1 sa positions. The last element is not corresponding to any character
    // ex firsts[0] => sa location of starting terminator
    vector<SuffixArrayIdx> firsts;
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

struct SuffixArrayNodeExtension_NEW {
    SuffixArrayNode_NEW node;
    // left extended intervals. Each means sa intervals of cW for each character c.
    vector<SuffixArrayNode_NEW> c_nodes;
    // vector<bool> c_nodes_exist;
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

void push_node_to_left_maximal(FmIndex &index, SuffixArrayNode_NEW &node, int &cnt){
    vector<uint64_t> ext_ranks(index.characters_cnt+1);
    vector<SuffixArrayNode_NEW> c_nodes(index.characters_cnt);
    vector<bool> c_nodes_exist(index.characters_cnt);
    bool left_maximal=false;
    while(!left_maximal){
        left_maximal=true;
        for(int c=0; c<index.characters_cnt; c++){
            index.STRING->ranks(c, node.firsts, ext_ranks);
            // interval empty
            if(ext_ranks[0]==ext_ranks[index.characters_cnt]){
                c_nodes_exist[c]=false;
                continue;
            }
            c_nodes_exist[c]=true;
            c_nodes[c]=SuffixArrayNode_NEW();
            c_nodes[c].firsts=vector<SuffixArrayIdx>(index.characters_cnt+1);
            for(int i=0; i<index.characters_cnt+1; i++){
                c_nodes[c].firsts[i]=index.C[c]+ext_ranks[i];
                c_nodes[c].depth=node.depth+1;
            }
            assert(node.interval_size() >= c_nodes[c].interval_size());
            if(c!=0 && node.interval_size()==c_nodes[c].interval_size()){
                // left non-maximal
                node = c_nodes[c];
                left_maximal=false;
                cnt++;
                break;
            }
        }
    }
}
/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
SuffixArrayNodeExtension_NEW extend_node_new(FmIndex &index, SuffixArrayNode_NEW &node){
    vector<RankArray> left_p_ranks(index.characters_cnt+1);
    RankArray rank_array(index.characters_cnt);
    for(int i=0; i<index.characters_cnt+1; i++){
        if(i==0 || node.firsts[i-1]!=node.firsts[i]){
            index.STRING->ranks_new(node.firsts[i], rank_array);
        } 
        left_p_ranks[i]=rank_array;
    }

    vector<SuffixArrayNode_NEW> c_nodes(index.characters_cnt);
    for(int c=0; c < index.characters_cnt; c++){
        vector<uint64_t> firsts(index.characters_cnt+1);
        for(int i=0; i<index.characters_cnt+1; i++){
            uint64_t sa_idx = index.C[c]+left_p_ranks[i][c];
            firsts[i]=sa_idx;
        }
        c_nodes[c]={firsts, node.depth+1};
    }

    vector<optional<SuffixArrayIdx>> prev_c_firsts;
    for(int c=0; c < c_nodes.size(); c++){
        SuffixArrayNode_NEW c_node = c_nodes[c];

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

    vector<SuffixArrayIdx> both_ext_terms;
    if(c_nodes[0].firsts[0]!=c_nodes[0].firsts[1]){
        // both ends are termination. collect all
        for(SuffixArrayIdx sa_idx=c_nodes[0].firsts[0]; sa_idx< c_nodes[0].firsts[1]; sa_idx++){
            //revert the rank process
            uint64_t rank = sa_idx - index.C[0] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, 0);
            both_ext_terms.push_back(prev_sa_idx);
        }
    }

    return {node, c_nodes, prev_c_firsts, both_ext_terms};
};

/*
* functions for suffix tree navigation
*/
SuffixArrayNode_NEW get_root_new(FmIndex &index){
    vector<uint64_t> firsts;
    for(auto rank: index.C){
        firsts.push_back(rank);
    }
    firsts.push_back(index.STRING->size());
    return {firsts, 0};
};

template<class T> using NodeFunc_NEW = optional<T>(*)(SuffixArrayNodeExtension_NEW&);

template<class T, NodeFunc_NEW<T> process_node>
void navigate_tree_new(SuffixArrayNode_NEW &root, int Lmin, FmIndex &fm_idx, vector<T> &Ts){
    Lmin = Lmin >= 1? Lmin : 1;

    std::stack<SuffixArrayNode_NEW> stack;
    
    int node_cnt=0;
    int process_cnt=0;
    int not_left_maximal_cnt=0;
    int left_push_cnt=0;

    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();
        
        push_node_to_left_maximal(fm_idx, node, left_push_cnt);
        auto ext = extend_node_new(fm_idx, node);
        if(ext.node.depth>=Lmin){
            optional<T> t = process_node(ext);
            if(t.has_value()) {
                Ts.push_back(t.value());
            }
            process_cnt++;
        }
        if(!ext.left_maximal()){
            not_left_maximal_cnt++;
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
    cout << "left pushed " << left_push_cnt << ", non left maximal: " << not_left_maximal_cnt << endl;
}


vector<SuffixArrayNode_NEW> collect_nodes(SuffixArrayNode_NEW root, FmIndex &fm_idx, int depth_max){
    std::stack<SuffixArrayNode_NEW> stack;
    vector<SuffixArrayNode_NEW> nodes;

    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();
        if(node.depth>=depth_max){
            nodes.push_back(node);
        } else {
            auto ext = extend_node_new(fm_idx, node);
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


#endif /* FM_INDEX_TREE_NEW_HPP_ */