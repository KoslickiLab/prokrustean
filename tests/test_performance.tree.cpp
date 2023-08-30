#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

vector<SuffixArrayNode> _extend_node__only_interval(WaveletString &str, CArray C, SuffixArrayNode &node){
    vector<RankArray> left_p_ranks;
    RankArray rank_array;
    for(int i=0; i<node.firsts.size(); i++){
        if(i==0 || (i==node.firsts.size()-1 && node.firsts[0]!=node.firsts[i])){
            rank_array = str.ranks(node.firsts[i]);  
        } 
        left_p_ranks.push_back(rank_array);
    }

    vector<SuffixArrayNode> c_nodes;
    for(int c=0; c < str.get_characters().size(); c++){
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        c_nodes.push_back({firsts, node.depth+1});
    }
    return c_nodes;
};


vector<SuffixArrayNode> _extend_node__wt(WaveletString &str, CArray C, SuffixArrayNode &node){
    vector<RankArray> left_p_ranks;
    RankArray rank_array;
    for(int i=0; i<node.firsts.size(); i++){
        if(i==0 || node.firsts[i-1]!=node.firsts[i]){
            rank_array = str.ranks(node.firsts[i]);  
        } 
        left_p_ranks.push_back(rank_array);
    }

    vector<SuffixArrayNode> c_nodes;
    for(int c=0; c < str.get_characters().size(); c++){
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        c_nodes.push_back({firsts, node.depth+1});
    }
    return c_nodes;
};

/* extend scenarios with left firsts */
vector<SuffixArrayNode> _extend_node__full(FmIndex &index, SuffixArrayNode &node){
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
    if(prev_c_firsts.size()>0 || both_ext_terms.size()>0){
        cout << ""; // dummy to make compiler not skip the functionality.
    }

    return c_nodes;
};

vector<SuffixArrayNode> _extend_node__c_node_only(FmIndex &index, SuffixArrayNode &node){
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
    return c_nodes;
};

vector<SuffixArrayNode> _extend_node__full_but_prev_l(FmIndex &index, SuffixArrayNode &node){
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

    return c_nodes;
};


void _navigate_tree(SuffixArrayNode &root, FmIndex &fm_idx, WaveletString &wt, int scenario_no){
    int node_cnt = 0;

    std::stack<SuffixArrayNode> stack;
    
    stack.push(root);
    while(!stack.empty()){
        auto node = stack.top();
        stack.pop();

        vector<SuffixArrayNode> c_nodes;
        switch (scenario_no)
        {
        case 0: // full extend
            c_nodes = _extend_node__full(fm_idx, node);
            break;
        case 1:
            c_nodes = _extend_node__full_but_prev_l(fm_idx, node);
            break;
        case 2:
            c_nodes = _extend_node__c_node_only(fm_idx, node);
            break;
        case 3:
            c_nodes = _extend_node__wt(wt, fm_idx.C, node);
            break;
        case 4:
            c_nodes = _extend_node__only_interval(wt, fm_idx.C, node);
            break;
        default:
            assert(false);
        }
        
        for(int i=0; i<c_nodes.size(); i++){
            // terminal
            if(i==0) continue;

            auto child = c_nodes[i];
            if(child.right_maximal()){
                stack.push(child);
            }
        }
        node_cnt++;
    }

    cout << "node visiteds: " << node_cnt << endl;
}


void _navigate_tree_components_record(SuffixArrayNode &root, FmIndex &fm_idx, wt_blcd<> &wt){
    int node_cnt = 0;
    int max_node_cnt = 0;
    int maximal_cnt = 0;
    int c_max = fm_idx.characters.size();
    int first_max = root.firsts.size();
    auto start_extend_node = std::chrono::steady_clock::now();
    auto start_extend_node_stage = std::chrono::steady_clock::now();
    auto start_process_node = std::chrono::steady_clock::now();
    auto extend_node_acc = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto extend_node_acc1 = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto extend_node_acc2 = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto extend_node_acc3 = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto extend_node_acc4 = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto process_node_acc = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();
    auto process_node_acc1 = std::chrono::steady_clock::now()-std::chrono::steady_clock::now();

    vector<MaximalRepeatAnnotation> annots;

    vector<optional<SuffixArrayNode>> nodes(1);
    vector<int> available;
    // std::iota(available.begin(), available.end(), 1);
    std::stack<int> stack;
    int curr_idx=0;
    nodes[curr_idx]=root;
    stack.push(curr_idx);
    
    // std::stack<SuffixArrayNode> stack;
    // stack.push(root);
    while(!stack.empty()){
        // auto node = stack.top();
        auto curr_idx = stack.top();
        SuffixArrayNode node = nodes[curr_idx].value();
        stack.pop();

        start_extend_node = std::chrono::steady_clock::now();
        start_extend_node_stage = std::chrono::steady_clock::now();
        vector<RankArray> left_p_ranks(node.firsts.size());
        RankArray rank_array(c_max);
        for(int i=0; i<node.firsts.size(); i++){
            if(i==0 || node.firsts[i-1]!=node.firsts[i]){
                rank_array = fm_idx.STRING->ranks(node.firsts[i]);  
            } 
            left_p_ranks[i]=rank_array;
        }
        extend_node_acc1+= std::chrono::steady_clock::now()-start_extend_node_stage;

        start_extend_node_stage = std::chrono::steady_clock::now();
        vector<SuffixArrayNode> c_nodes(fm_idx.characters.size());
        for(int c=0; c < fm_idx.characters.size(); c++){
            vector<uint64_t> firsts(left_p_ranks.size());
            for(int f=0; f< left_p_ranks.size(); f++){
                uint64_t sa_idx = fm_idx.C[c]+left_p_ranks[f][c];
                firsts[f]=sa_idx;
            }
            c_nodes[c] = {firsts, node.depth+1};
        }
        extend_node_acc2+= std::chrono::steady_clock::now()-start_extend_node_stage;

        start_extend_node_stage = std::chrono::steady_clock::now();
        vector<optional<SuffixArrayIdx>> prev_c_firsts(c_nodes.size());
        for(int c=0; c < c_nodes.size(); c++){
            SuffixArrayNode c_node = c_nodes[c];

            if(c_node.interval_size()==0){
                prev_c_firsts[c]=nullopt;
                continue;
            }
            //revert the rank process
            SuffixArrayIdx sa_idx = c_node.get_valid_first();
            uint64_t rank = sa_idx - fm_idx.C[c] + 1;
            SuffixArrayIdx prev_sa_idx = fm_idx.STRING->select(rank, c);
            // cout << "obtained " << prev_sa_idx << " for " << sa_idx << endl;
            prev_c_firsts[c]=prev_sa_idx;
        }
        extend_node_acc3+= std::chrono::steady_clock::now()-start_extend_node_stage;

        start_extend_node_stage = std::chrono::steady_clock::now();
        vector<SuffixArrayIdx> both_ext_terms;
        if(c_nodes[0].firsts[0]!=c_nodes[0].firsts[1]){
            // both ends are termination. collect all
            for(SuffixArrayIdx sa_idx=c_nodes[0].firsts[0]; sa_idx< c_nodes[0].firsts[1]; sa_idx++){
                //revert the rank process
                uint64_t rank = sa_idx - fm_idx.C[0] + 1;
                SuffixArrayIdx prev_sa_idx = fm_idx.STRING->select(rank, 0);
                both_ext_terms.push_back(prev_sa_idx);
            }
        }
        extend_node_acc4+= std::chrono::steady_clock::now()-start_extend_node_stage;

        SuffixArrayNodeExtension ext = {node, c_nodes, prev_c_firsts, both_ext_terms};
        extend_node_acc+= std::chrono::steady_clock::now()-start_extend_node;

        start_process_node = std::chrono::steady_clock::now();
        auto annot = get_repeat_annotations(ext);
        if(annot.has_value()){
            maximal_cnt++;
            annots.push_back(annot.value());
        }
        process_node_acc+= std::chrono::steady_clock::now()-start_process_node;

        start_process_node = std::chrono::steady_clock::now();
        if(ext.left_maximal() && ext.node.right_maximal()){
            auto repr_extensions = decide_repr_sa_extensions(ext.c_nodes.size(), ext.distinct_extensions());
            maximal_cnt++;
        }
        process_node_acc1+= std::chrono::steady_clock::now()-start_process_node;

        for(int i=0; i<ext.c_nodes.size(); i++){
            // terminal
            if(i==0) continue;

            auto child = ext.c_nodes[i];
            if(child.right_maximal()){
                // stack.push(child);

                if(available.size()>0){
                    auto c_idx = available[available.size()-1];
                    available.pop_back();
                    stack.push(c_idx);
                    nodes[c_idx]=child;
                } else {
                    nodes.push_back(child);
                    stack.push(nodes.size()-1);
                }
            }
        }
        available.push_back(curr_idx);
        nodes[curr_idx] = nullopt;
        max_node_cnt = max_node_cnt<nodes.size()? nodes.size(): max_node_cnt;
        node_cnt++;
    }
    cout << "max node count: " << max_node_cnt;
    cout << "node visiteds: " << node_cnt;
    cout << ", extend_node: " << extend_node_acc.count()/1000000 << "ms";
    cout << ", process_node: " << process_node_acc.count()/1000000 << "ms" << endl;
    cout << "ext1: " << extend_node_acc1.count()/1000000 << "ms";
    cout << ", ext2: " << extend_node_acc2.count()/1000000 << "ms";
    cout << ", ext3: " << extend_node_acc3.count()/1000000 << "ms";
    cout << ", ext4: " << extend_node_acc4.count()/1000000 << "ms";
    cout << ", process decide: " << process_node_acc1.count()/1000000 << "ms";
    
    cout << endl;
    cout << "dummy" << maximal_cnt;
}


void test_compare_various_tree_exploration(){
    int Lmin = 1;
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_SEQ, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;

    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode root = get_root(fm_idx);
    
    start = std::chrono::steady_clock::now();
    _navigate_tree(root, fm_idx, str, 0);
    cout << "navigate tree with full extend: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
    start = std::chrono::steady_clock::now();
    _navigate_tree(root, fm_idx, str, 1);
    cout << "navigate tree with full but prev: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    _navigate_tree(root, fm_idx, str, 2);
    cout << "navigate tree c_node only: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    _navigate_tree(root, fm_idx, str, 3);
    cout << "navigate tree wt : " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
    start = std::chrono::steady_clock::now();
    _navigate_tree(root, fm_idx, str, 4);
    cout << "navigate tree interval : " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void test_tree_exploration_investigation(){
    int Lmin = 1;
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_SEQ, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;

    auto start = std::chrono::steady_clock::now();
    SuffixArrayNode root = get_root(fm_idx);
    _navigate_tree_components_record(root, fm_idx, str.wt);
    cout << "whole navi : " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    vector<MaximalRepeatAnnotation> repeats;
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    cout << "original navi : " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void main_performance_tree() {
    // test_compare_various_tree_exploration();
    test_tree_exploration_investigation();
}
