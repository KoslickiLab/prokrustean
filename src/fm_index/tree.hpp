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
struct SuffixArrayInterval {
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

    uint64_t count(){
        return firsts[firsts.size()-1]-firsts[0];
    }

    uint64_t suffix_count(){
        return firsts[1]-firsts[0];
    }

    bool right_maximal(){
        return distinct_extensions().size()>1 || suffix_count()>1;
    }
};

struct SuffixArrayNode {
    SuffixArrayInterval interval;
    // left extended intervals. Each means sa intervals of cW for each character c.
    vector<SuffixArrayInterval> c_intervals;
    vector<CharId> c_s;

    uint64_t first_l(CharId l){
        // The first sa index always appears at the 0. 
        // Even if '#' and 'A' is not an extension but 'G' is, for example, then firsts[0]=firsts[1]=firsts[2]!=firsts[3].
        return c_intervals[l].firsts[0];
    }

    uint64_t first_r(CharId r){
        uint64_t min = std::numeric_limits<uint64_t>::max();
        for(int i=0; i<c_intervals.size(); i++){
            if(c_intervals[i].firsts[r]<min){
                min = c_intervals[i].firsts[r];
            }
        }
        return min;
    }

    vector<tuple<CharId, CharId, uint64_t>> distinct_extensions(){
        vector<tuple<CharId, CharId, uint64_t>> exts = {};
        for(int i=0; i<c_intervals.size(); i++){
            for(auto r: c_intervals[i].distinct_extensions()){
                exts.push_back(make_tuple(i, get<0>(r), get<1>(r)));    
            }
        }
        return exts;
    }

    uint64_t prefix_count(){
        return c_intervals[0].count();
    }

    bool dup_sequence_exists(){
        return c_intervals[0].suffix_count()>1;
    }

    bool left_maximal(){
        int distinct = 0;
        for(auto interval:c_intervals){
            if(interval.count()>0){
                distinct++;
            }
        }
        return distinct>1 || prefix_count()>1;
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
SuffixArrayNode to_node(FmIndex &index, SuffixArrayInterval &interval){
    vector<RankArray> left_p_ranks;
    RankArray rank_array;
    for(int i=0; i<interval.firsts.size(); i++){
        if(i==0 || interval.firsts[i-1]!=interval.firsts[i]){
            rank_array = index.STRING->ranks(interval.firsts[i]);  
        } 
        left_p_ranks.push_back(rank_array);
    }

    vector<SuffixArrayInterval> left_intervals;
    vector<CharId> c_s;
    for(int c=0; c < index.STRING->get_characters().size(); c++){
        uint64_t C = index.C[c];
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = index.C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        left_intervals.push_back({firsts, interval.depth+1});
        c_s.push_back(c);
    }
    
    return {interval, left_intervals, c_s};
};

// left_extension left_extend(FmIndex &index, Interval interval){
//     SuccintString STRING = index.STRING;
//     CArray C = index.C;

//     ParallelRank rank = STRING.parallel_rank(interval.first_TERM);
//     ParallelRank before_TERM = rank;

//     if(interval.first_A != interval.first_TERM) rank = STRING.parallel_rank(interval.first_A);
//     ParallelRank before_A = rank;

//     if(interval.first_C != interval.first_A) rank = STRING.parallel_rank(interval.first_C);
//     ParallelRank before_C = rank;

//     if(interval.first_G != interval.first_C) rank = STRING.parallel_rank(interval.first_G);
//     ParallelRank before_G = rank;

//     if(interval.first_T != interval.first_G) rank = STRING.parallel_rank(interval.first_T);
//     ParallelRank before_T = rank;

//     if(interval.last != interval.first_T) rank = STRING.parallel_rank(interval.last);
//     ParallelRank before_end = rank;

//     return {
//         {   before_TERM.TERM,    before_A.TERM,    before_C.TERM,    before_G.TERM,    before_T.TERM,    before_end.TERM, interval.depth+1},
//         {C.A + before_TERM.A, C.A + before_A.A, C.A + before_C.A, C.A + before_G.A, C.A + before_T.A, C.A + before_end.A, interval.depth+1},
//         {C.C + before_TERM.C, C.C + before_A.C, C.C + before_C.C, C.C + before_G.C, C.C + before_T.C, C.C + before_end.C, interval.depth+1},
//         {C.G + before_TERM.G, C.G + before_A.G, C.G + before_C.G, C.G + before_G.G, C.G + before_T.G, C.G + before_end.G, interval.depth+1},
//         {C.T + before_TERM.T, C.T + before_A.T, C.T + before_C.T, C.T + before_G.T, C.T + before_T.T, C.T + before_end.T, interval.depth+1}
//     };
// };

/*
* functions for suffix tree navigation
*/
SuffixArrayInterval get_root(FmIndex &index){
    vector<uint64_t> firsts;
    for(auto rank: index.C){
        firsts.push_back(rank);
    }
    firsts.push_back(index.STRING->size());
    return {firsts, 0};
};

template<class T> using NodeFunc = optional<T>(*)(SuffixArrayNode&, FmIndex&);

template<class T, NodeFunc<T> process_node>
vector<T> navigate_tree(SuffixArrayInterval &root, int Lmin, FmIndex &fm_idx){
    Lmin = Lmin >= 1? Lmin : 1;

    std::stack<SuffixArrayInterval> interval_stack;
    vector<T> Ts;
    
    interval_stack.push(root);
    while(!interval_stack.empty()){
        auto interval = interval_stack.top();
        interval_stack.pop();

        auto node = to_node(fm_idx, interval);
        if(node.interval.depth>=Lmin){
            optional<T> t = process_node(node, fm_idx);
            if(t.has_value()) {
                Ts.push_back(t.value());
            }
        }
        
        for(int i=0; i<node.c_intervals.size(); i++){
            // terminal
            if(node.c_s[i]==0) 
            continue;
            auto c_interval = node.c_intervals[i];
            if(c_interval.count()>1){
                interval_stack.push(c_interval);
            }
        }
    }
    return Ts;
}

#endif /* FM_INDEX_TREE_HPP_ */