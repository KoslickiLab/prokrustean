#ifndef FM_INDEX_TREE_NEW_HPP_
#define FM_INDEX_TREE_NEW_HPP_

#include "../data_types.hpp"
#include "../fm_index/index.hpp"
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
    vector<SuffixArrayIdx> firsts;
    SuffixArrayIdx valid_first;
    uint64_t depth;
    bool right_maximal;

    /*deprecated*/
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
    vector<CharId> left_paired_a_char;
    // For each character a, how many distinct Wa form exists.
    vector<int> right_cnts;
    // For each character a, example of c of cWa form.
    vector<CharId> right_paired_c_char;
    // For each c, if representative.
    vector<bool> left_repr;
    // For each a, if representative.
    vector<bool> right_repr;

    vector<tuple<SeqId, Pos>> locations;

    vector<SuffixArrayIdx> sa_indices;
};

struct TreeWorkspace {
    /* Includes all information required to process, so that memory is reused.*/
    TreeWorkspace(int characters_cnt, int thread_idx=0){
        c_nodes=vector<SuffixArrayNode>(characters_cnt);
        c_nodes_open=vector<bool>(characters_cnt);
        c_first_ranks= vector<uint64_t>(characters_cnt+1);
        // prev_c_firsts=vector<optional<SuffixArrayIdx>>(characters_cnt);
        for(int i=0; i<characters_cnt; i++){
            c_nodes[i].firsts=vector<SuffixArrayIdx>(characters_cnt+1);
        }
        this->thread_idx=thread_idx;
        this->characters_cnt=characters_cnt;
        this->repr_work.left_cnts=vector<int>(characters_cnt);
        this->repr_work.left_paired_a_char=vector<CharId>(characters_cnt);
        this->repr_work.right_cnts=vector<int>(characters_cnt);
        this->repr_work.right_paired_c_char=vector<CharId>(characters_cnt);
        this->repr_work.left_repr=vector<bool>(characters_cnt);
        this->repr_work.right_repr=vector<bool>(characters_cnt);
    }

    // The node in interest
    SuffixArrayNode node;
    // left extensions. Each means sa intervals of cW for each character c.
    vector<SuffixArrayNode> c_nodes;
    // whether the corresponding c_node is open (the branch exists)
    vector<bool> c_nodes_open;
    
    ReprSuffixArrayIndexWorkspace repr_work;

    StratumId stratum_id;

    int thread_idx;
    // alphabet size
    int characters_cnt;
    // temporary first ranks gathered when c_nodes are computed
    vector<uint64_t> c_first_ranks;
    
    // Each means the first sa index of W that corresponds to first sa index of cW for each character c.
    // vector<optional<SuffixArrayIdx>> prev_c_firsts;
    // suffixes such that both left and right extensions are terminations
    vector<SuffixArrayIdx> both_ext_terms;

    // for profiling purposes
    vector<uint64_t> any_measure;

    SuffixArrayIdx first_r(CharId c){
        return node.firsts[c];
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
void extend_node(FmIndex &index, TreeWorkspace &ext){
    // auto start = std::chrono::steady_clock::now();

    // push node to left maximal. everytime only one left child found, update node.
    bool left_maximal=false;
    bool left_maximal_checked=false;
    bool same_interval_size_with_the_child=false;
    int right_distinct=0;
    uint64_t amount_of_interval_left_to_be_captured=0;
    while(!left_maximal){
        left_maximal=true;
        left_maximal_checked=false;
        // start from node interval
        amount_of_interval_left_to_be_captured=ext.node.firsts[index.characters_cnt]-ext.node.firsts[0];
        for(auto c: index.characters_ranked_by_abundance){
            if(amount_of_interval_left_to_be_captured==0){
                ext.c_nodes_open[c]=false;
                continue;
            }

            // auto start2 = std::chrono::steady_clock::now();
            index.STRING->ranks(c, ext.node.firsts, ext.c_first_ranks);
            // ext.any_measure[1]+=(std::chrono::steady_clock::now()-start2).count();
            
            // interval is empty -> left node for c is empty
            if(ext.c_first_ranks[0]==ext.c_first_ranks[index.characters_cnt]){
                ext.c_nodes_open[c]=false;
                continue;
            } else {
                amount_of_interval_left_to_be_captured -= ext.c_first_ranks[ext.characters_cnt]-ext.c_first_ranks[0];
                assert(amount_of_interval_left_to_be_captured>=0);
            }
            // check if left non-maximal. if yes, then casually move to the only branch.
            if(c!=0 && !left_maximal_checked){
                bool same_interval_size_with_the_child = ext.node.firsts[ext.characters_cnt]-ext.node.firsts[0] == ext.c_first_ranks[ext.characters_cnt]-ext.c_first_ranks[0];
                if(same_interval_size_with_the_child){
                    // move to the child directly
                    ext.node.depth=ext.node.depth+1;
                    for(int i=0; i<index.characters_cnt+1; i++){
                        ext.node.firsts[i]=index.C[c]+ext.c_first_ranks[i];
                    }
                    left_maximal=false;
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
SuffixArrayNode get_root(FmIndex &index){
    vector<uint64_t> firsts;
    for(auto rank: index.C){
        firsts.push_back(rank);
    }
    firsts.push_back(index.STRING->size());
    return {firsts, 0};
};

template<class T> using NodeFunc_NEW = void(*)(FmIndex&, TreeWorkspace&, T&);

template<class T, NodeFunc_NEW<T> process_node>
void navigate_strata(SuffixArrayNode &root, int Lmin, FmIndex &fm_idx, T &t, int thread_idx=0, bool verbose=false){
    Lmin = Lmin >= 1? Lmin : 1;
    // assert(fm_idx.locator!=nullptr);
    // cout << "warning: sample x" << endl;

    TreeWorkspace ext(fm_idx.characters_cnt, thread_idx);
    ext.any_measure=vector<uint64_t>(10,0);

    auto start = std::chrono::steady_clock::now();

    stack<SuffixArrayNode> stack;
    stack.push(root);
    while(!stack.empty()){
        if(verbose) start = std::chrono::steady_clock::now();
        
        ext.node = stack.top();
        stack.pop();
        
        if(verbose) ext.any_measure[1]+=(std::chrono::steady_clock::now()-start).count();
        if(verbose) start = std::chrono::steady_clock::now();
        
        extend_node(fm_idx, ext);
        
        if(verbose) ext.any_measure[2]+=(std::chrono::steady_clock::now()-start).count();
        if(verbose) start = std::chrono::steady_clock::now();

        if(ext.node.depth>=Lmin){    
            process_node(fm_idx, ext, t);
        }
        
        if(verbose) ext.any_measure[3]+=(std::chrono::steady_clock::now()-start).count();
        if(verbose) start = std::chrono::steady_clock::now();

        for(int i=0; i<fm_idx.characters_cnt; i++){
            if(i==0 // skip terminal
                || !ext.c_nodes_open[i] // not existing branch
                || !ext.c_nodes[i].right_maximal // not right maximal
                ) continue;

            stack.push(ext.c_nodes[i]);
        }
        
        if(verbose) ext.any_measure[1]+=(std::chrono::steady_clock::now()-start).count();
    }
    if(verbose){
        for(int i=0; i<ext.any_measure.size(); i++){
            cout << "measure " << i << ": " << ext.any_measure[i]/1000000 << "ms" << endl;
        }
    }
}

template<class T, NodeFunc_NEW<T> process_node>
vector<SuffixArrayNode> collect_roots_while_navigate_strata(int Lmin, FmIndex &fm_idx, T &t, int depth_max){
    assert(Lmin>=1 && depth_max >=1);

    TreeWorkspace ext(fm_idx.characters_cnt);
    vector<SuffixArrayNode> roots;
    
    stack<SuffixArrayNode> stack;
    stack.push(get_root(fm_idx));
    while(!stack.empty()){
        ext.node = stack.top();
        stack.pop();
        
        if(ext.node.depth>=depth_max){
            roots.push_back(ext.node);
        } else {
            extend_node(fm_idx, ext);
            if(ext.node.depth>=Lmin){    
                process_node(fm_idx, ext, t);
            }
        
            for(int i=0; i<fm_idx.characters_cnt; i++){
                if(i==0 // skip terminal
                || !ext.c_nodes_open[i] // not existing branch
                || !ext.c_nodes[i].right_maximal // not right maximal
                ) continue;

                stack.push(ext.c_nodes[i]);
            }
        }
    }
    return roots;
}


#endif /* FM_INDEX_TREE_NEW_HPP_ */