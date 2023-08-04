#include "index.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_TREE_HPP_
#define FM_INDEX_TREE_HPP_

/*
 * representation of a right-maximal substring (SA node) as a list of STRING intervals
 */
struct Interval {
    // |characters| + 1 sa positions. The last element is not corresponding to any character
    // ex firsts[0] => sa location of starting terminator
    vector<uint64_t> firsts;
	uint64_t depth;

    vector<CharId> collect_extensions(){
        vector<CharId> exts = {};
        for(int i=0; i<firsts.size(); i++){
            if(i<firsts.size()-1){
                if(firsts[i+1]-firsts[i]) 
                exts.push_back(i);
            }
        }
        return exts;
    }
};

struct left_extension {

    // left extended intervals. Each means sa intervals of cW for each character c.
    vector<Interval> intervals;

    uint64_t first(CharId l, CharId r){
        return intervals[l].firsts[r];
    }

    vector<tuple<CharId,CharId>> collect_extentions(){
        vector<tuple<CharId,CharId>> exts = {};
        for(int i=0; i<intervals.size(); i++){
            for(auto r: intervals[i].collect_extensions()){
                exts.push_back(make_tuple(i, r));    
            }
        }
        return exts;
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
left_extension left_extend(FmIndex &index, Interval interval){
    vector<RankArray> left_p_ranks;
    RankArray rank_array;
    for(int i=0; i<interval.firsts.size(); i++){
        if(i==0 || interval.firsts[i-1]!=interval.firsts[i]){
            rank_array = index.STRING->ranks(interval.firsts[i]);  
        } 
        left_p_ranks.push_back(rank_array);
    }

    vector<Interval> left_intervals;
    for(int c=0; c < index.STRING->get_characters().size(); c++){
        uint64_t C = index.C[c];
        vector<uint64_t> firsts;
        for(auto p_rank:left_p_ranks){
            uint64_t sa_idx = index.C[c]+p_rank[c];
            firsts.push_back(sa_idx);
        }
        left_intervals.push_back({firsts, interval.depth+1});
    }
    
    return {left_intervals};
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
Interval get_root(FmIndex &index){
    vector<uint64_t> firsts;
    for(auto rank: index.C){
        firsts.push_back(rank);
    }
    firsts.push_back(index.STRING->size());
    return {firsts, 0};
};

#endif /* FM_INDEX_TREE_HPP_ */