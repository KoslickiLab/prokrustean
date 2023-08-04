
#include "rank.hpp"
#include "index.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_TREE_HPP_
#define FM_INDEX_TREE_HPP_

/*
 * representation of a right-maximal substring (SA node) as a list of STRING intervals
 */
struct Interval {

	uint64_t first_TERM;
	uint64_t first_A;
	uint64_t first_C;
	uint64_t first_G;
	uint64_t first_T;
	uint64_t last;

	//depth = |W|
	uint64_t depth;

    uint64_t first(char c){
        uint64_t idx;
        switch (c)
        {
        case 'A': idx = first_A; break;
        case 'C': idx = first_C; break;
        case 'G': idx = first_G; break;
        case 'T': idx = first_T; break;
        default: idx = first_TERM; break;
        }
        return idx;
    }

	uint64_t key(){
		return first_TERM;
	}

    // uint64_t occurrence_cnt(char c){
    //     uint64_t cnt;
    //     switch (c)
    //     {
    //     case 'A': cnt = first_C-first_A; break;
    //     case 'C': cnt = first_G-first_C; break;
    //     case 'G': cnt = first_T-first_G; break;
    //     case 'T': cnt = last-first_T; break;
    //     // term
    //     default: cnt = first_A-first_TERM; break;
    //     }
    //     return cnt;
    // }
    vector<char> collect_extensions(){
        vector<char> exts = {};
        if(first_A!=first_C) exts.push_back('A');
        if(first_C!=first_G) exts.push_back('C');
        if(first_G!=first_T) exts.push_back('G');
        if(first_T!=last) exts.push_back('T');
        if(first_TERM!=first_A) exts.push_back(0);
        return exts;
    }
};

struct left_extension {

    Interval TERM;
	Interval A;
	Interval C;
	Interval G;
	Interval T;

    uint64_t first(char l, char r){
        uint64_t idx;
        switch (l)
        {
        case 'A': idx = A.first(r); break;
        case 'C': idx = C.first(r); break;
        case 'G': idx = G.first(r); break;
        case 'T': idx = T.first(r); break;
        default: idx = TERM.first(r); break;
        }
        return idx;
    }

    vector<tuple<char,char>> collect_extentions(){
        vector<tuple<char,char>> exts = {};
        for (auto r: TERM.collect_extensions()){
            exts.push_back(make_tuple(0, r));
        }
        for (auto r: A.collect_extensions()){
            exts.push_back(make_tuple('A', r));
        }
        for (auto r: C.collect_extensions()){
            exts.push_back(make_tuple('C', r));
        }
        for (auto r: G.collect_extensions()){
            exts.push_back(make_tuple('G', r));
        }
        for (auto r: T.collect_extensions()){
            exts.push_back(make_tuple('T', r));
        }
        
        return exts;
    }
};

/*
* Input: suffix tree node N.
* Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
*/
left_extension left_extend(FmIndex &index, Interval interval){
    SuccintString STRING = index.STRING;
    CArray C = index.C;

    ParallelRank rank = STRING.parallel_rank(interval.first_TERM);
    ParallelRank before_TERM = rank;

    if(interval.first_A != interval.first_TERM) rank = STRING.parallel_rank(interval.first_A);
    ParallelRank before_A = rank;

    if(interval.first_C != interval.first_A) rank = STRING.parallel_rank(interval.first_C);
    ParallelRank before_C = rank;

    if(interval.first_G != interval.first_C) rank = STRING.parallel_rank(interval.first_G);
    ParallelRank before_G = rank;

    if(interval.first_T != interval.first_G) rank = STRING.parallel_rank(interval.first_T);
    ParallelRank before_T = rank;

    if(interval.last != interval.first_T) rank = STRING.parallel_rank(interval.last);
    ParallelRank before_end = rank;

    return {
        {   before_TERM.TERM,    before_A.TERM,    before_C.TERM,    before_G.TERM,    before_T.TERM,    before_end.TERM, interval.depth+1},
        {C.A + before_TERM.A, C.A + before_A.A, C.A + before_C.A, C.A + before_G.A, C.A + before_T.A, C.A + before_end.A, interval.depth+1},
        {C.C + before_TERM.C, C.C + before_A.C, C.C + before_C.C, C.C + before_G.C, C.C + before_T.C, C.C + before_end.C, interval.depth+1},
        {C.G + before_TERM.G, C.G + before_A.G, C.G + before_C.G, C.G + before_G.G, C.G + before_T.G, C.G + before_end.G, interval.depth+1},
        {C.T + before_TERM.T, C.T + before_A.T, C.T + before_C.T, C.T + before_G.T, C.T + before_T.T, C.T + before_end.T, interval.depth+1}
    };
};

/*
* functions for suffix tree navigation
*/
Interval get_root(FmIndex &index){

    return {
        0,
        index.C.A,
        index.C.C,
        index.C.G,
        index.C.T,
        index.STRING.size(),
        0
    };

};

#endif /* FM_INDEX_TREE_HPP_ */