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
    SuffixArrayIdx valid_first;
    uint64_t depth;
    bool right_maximal;

    SuffixArrayIdx get_valid_first(){
        for(int i=0; i<firsts.size()-1; i++){
            if(firsts[i]<firsts[i+1]) 
                return firsts[i];
        }
        assert(false); //interval might be 0 - invalid
    }
};

struct ReprSuffixArrayIndexWorkspace {
    // For each character c, how many distinct cW form exists. e.g. AWT, AWG -> left cnt is 2
    vector<int> left_cnts;
    // For each character c, example of a of cWa form - only meaningful when exclusive
    vector<CharId> left_paired_example;
    // For each character a, how many distinct Wa form exists.
    vector<int> right_cnts;
    // For each character a, example of c of cWa form.
    vector<CharId> right_paired_example;
    // For each c, if representative.
    vector<bool> left_repr;
    // For each a, if representative.
    vector<bool> right_repr;
};

struct SuffixArrayNodeExtension_NEW {
    /* Includes all information required to process, so that memory is reused.*/
    SuffixArrayNodeExtension_NEW(int characters_cnt){
        c_nodes=vector<SuffixArrayNode_NEW>(characters_cnt);
        c_nodes_open=vector<bool>(characters_cnt);
        c_first_ranks= vector<uint64_t>(characters_cnt+1);
        // prev_c_firsts=vector<optional<SuffixArrayIdx>>(characters_cnt);
        for(int i=0; i<characters_cnt; i++){
            c_nodes[i].firsts=vector<SuffixArrayIdx>(characters_cnt+1);
        }
        this->characters_cnt=characters_cnt;
        this->repr_sa_workspace.left_cnts=vector<int>(characters_cnt);
        this->repr_sa_workspace.left_paired_example=vector<CharId>(characters_cnt);
        this->repr_sa_workspace.right_cnts=vector<int>(characters_cnt);
        this->repr_sa_workspace.right_paired_example=vector<CharId>(characters_cnt);
        this->repr_sa_workspace.left_repr=vector<bool>(characters_cnt);
        this->repr_sa_workspace.right_repr=vector<bool>(characters_cnt);
    }

    // The node in interest
    SuffixArrayNode_NEW node;
    // left extensions. Each means sa intervals of cW for each character c.
    vector<SuffixArrayNode_NEW> c_nodes;
    // whether the corresponding c_node is open (the branch exists)
    vector<bool> c_nodes_open;
    
    ReprSuffixArrayIndexWorkspace repr_sa_workspace;

    // alphabet size
    int characters_cnt;
    // temporary first ranks gathered when c_nodes are computed
    vector<uint64_t> c_first_ranks;
    
    // Each means the first sa index of W that corresponds to first sa index of cW for each character c.
    // vector<optional<SuffixArrayIdx>> prev_c_firsts;
    // suffixes such that both left and right extensions are terminations
    vector<SuffixArrayIdx> both_ext_terms;
    
    // for profiling purposes
    vector<int> any_measure;

    // SuffixArrayIdx first_l(CharId c){
    //     assert(prev_c_firsts[c].has_value());
    //     // cout << "first l: " << prev_c_firsts[c].value() << endl;
    //     return prev_c_firsts[c].value();
    // }

    SuffixArrayIdx first_r(CharId c){
        return node.firsts[c];
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
void extend_node_new(FmIndex &index, SuffixArrayNodeExtension_NEW &ext){
    // auto start = std::chrono::steady_clock::now();

    // push node to left maximal. everytime only one left child found, update node.
    bool left_maximal=false;
    bool left_maximal_checked=false;
    bool same_interval_size_with_the_child=false;
    int right_distinct=0;
    while(!left_maximal){
        left_maximal=true;
        left_maximal_checked=false;
        for(int c=0; c<index.characters_cnt; c++){
            // this function updates c_first_ranks that is, for a target left c, the list of ranks of firsts for each right character.
            index.STRING->ranks(c, ext.node.firsts, ext.c_first_ranks);
            // interval is empty -> left node for c is empty
            if(ext.c_first_ranks[0]==ext.c_first_ranks[index.characters_cnt]){
                ext.c_nodes_open[c]=false;
                continue;
            }
            // check if left non-maximal
            if(c!=0 && !left_maximal_checked){
                bool same_interval_size_with_the_child = ext.node.firsts[ext.characters_cnt]-ext.node.firsts[0] == ext.c_first_ranks[ext.characters_cnt]-ext.c_first_ranks[0];
                if(same_interval_size_with_the_child){
                    // move to the child directly
                    ext.node.depth=ext.node.depth+1;
                    for(int i=0; i<index.characters_cnt+1; i++){
                        ext.node.firsts[i]=index.C[c]+ext.c_first_ranks[i];
                    }
                    left_maximal=false;
                    // ext.any_measure[0]++;
                    // jump to next phase
                    break;
                } else {
                    left_maximal_checked=true;
                }
            }
            
            // update c_nodes, especially right maximal
            ext.c_nodes_open[c]=true;
            ext.c_nodes[c].right_maximal=false;
            ext.c_nodes[c].depth=ext.node.depth+1;
            right_distinct=0;
            for(int i=0; i<index.characters_cnt+1; i++){
                if(i>0 && ext.c_first_ranks[i] == ext.c_first_ranks[i-1]){
                    ext.c_nodes[c].firsts[i]=ext.c_nodes[c].firsts[i-1];
                } else {
                    ext.c_nodes[c].firsts[i]=index.C[c]+ext.c_first_ranks[i];
                    if(i>0) right_distinct++;
                }
            }
            // compute right maximal
            if(right_distinct>1 || ext.c_nodes[c].firsts[1] - ext.c_nodes[c].firsts[0]>1){
                ext.c_nodes[c].right_maximal=true;
            }
        }
    }

    // ext.any_measure[5]+=(std::chrono::steady_clock::now()-start).count();

    // for(int c=0; c < index.characters_cnt; c++){
    //     if(!ext.c_nodes_open[c]){
    //         ext.prev_c_firsts[c]=nullopt;
    //         continue;
    //     } 

    //     //revert the rank process
    //     SuffixArrayIdx sa_idx = ext.c_nodes[c].get_valid_first();
    //     uint64_t rank = sa_idx - index.C[c] + 1;
    //     SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, c);
    //     // cout << "obtained " << prev_sa_idx << " for " << sa_idx << endl;
    //     ext.prev_c_firsts[c]=prev_sa_idx;
    // }
    
    ext.both_ext_terms.clear();
    if(ext.c_nodes_open[0] && ext.c_nodes[0].firsts[0]!=ext.c_nodes[0].firsts[1]){
        // both ends are termination. collect all
        for(SuffixArrayIdx sa_idx=ext.c_nodes[0].firsts[0]; sa_idx< ext.c_nodes[0].firsts[1]; sa_idx++){
            //revert the rank process
            uint64_t rank = sa_idx - index.C[0] + 1;
            SuffixArrayIdx prev_sa_idx = index.STRING->select(rank, 0);
            ext.both_ext_terms.push_back(prev_sa_idx);
        }
    }
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

template<class T> using NodeFunc_NEW = void(*)(FmIndex&, SuffixArrayNodeExtension_NEW&, vector<T>&);

template<class T, NodeFunc_NEW<T> process_node>
void navigate_tree_new(SuffixArrayNode_NEW &root, int Lmin, FmIndex &fm_idx, vector<T> &Ts){
    Lmin = Lmin >= 1? Lmin : 1;

    stack<SuffixArrayNode_NEW> stack;
    
    SuffixArrayNodeExtension_NEW ext(fm_idx.characters_cnt);
    ext.any_measure=vector<int>(10,0);
    // auto start = std::chrono::steady_clock::now();

    stack.push(root);
    while(!stack.empty()){
        // start = std::chrono::steady_clock::now();
        ext.node = stack.top();
        stack.pop();
        // ext.any_measure[1]+=(std::chrono::steady_clock::now()-start).count();
        
        
        // start = std::chrono::steady_clock::now();
        extend_node_new(fm_idx, ext);
        // ext.any_measure[2]+=(std::chrono::steady_clock::now()-start).count();
        if(ext.node.depth>=Lmin){
            // start = std::chrono::steady_clock::now();
            process_node(fm_idx, ext, Ts);
            // ext.any_measure[3]+=(std::chrono::steady_clock::now()-start).count();
            
        }
        
        // start = std::chrono::steady_clock::now();
        for(int i=0; i<fm_idx.characters_cnt; i++){
            // terminal
            if(i==0) continue;

            if(!ext.c_nodes_open[i]) continue;

            if(ext.c_nodes[i].right_maximal){
                stack.push(ext.c_nodes[i]);
            }
        }
        // ext.any_measure[1]+=(std::chrono::steady_clock::now()-start).count();
    }
    cout << "left pushed " << ext.any_measure[0] << endl;
    for(int i=0; i<ext.any_measure.size(); i++){
        cout << "measure " << i << ": " << ext.any_measure[i]/1000000 << "ms" << endl;
    }
}


vector<SuffixArrayNode_NEW> collect_nodes(SuffixArrayNode_NEW root, FmIndex &fm_idx, int depth_max){
    std::stack<SuffixArrayNode_NEW> stack;
    vector<SuffixArrayNode_NEW> nodes;
    int cnt=0;

    SuffixArrayNodeExtension_NEW ext(fm_idx.characters_cnt);

    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();
        if(node.depth>=depth_max){
            nodes.push_back(node);
        } else {
            extend_node_new(fm_idx, ext);
            for(int i=0; i<ext.c_nodes.size(); i++){
                // terminal
                if(i==0) continue;

                if(!ext.c_nodes_open[i]) continue;

                auto child = ext.c_nodes[i];
                if(child.right_maximal){
                    stack.push(child);
                }
            }
        }
    }
    return nodes;
}


#endif /* FM_INDEX_TREE_NEW_HPP_ */