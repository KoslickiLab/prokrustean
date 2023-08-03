

#include "rank.hpp"
#include "tree.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

/*
 * factory. For potential optimization later.
 */
// FmIndex create_fm_index(string path, char TERM = '#'){
// 	FmIndex idx = FmIndex(path, TERM);
// 	return idx;
// }

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

	uint64_t key(){
		return first_TERM;
	}
    bool has_term(){
        return first_TERM != first_A;
    }
    bool has_A(){
        return first_A != first_C;
    }
    bool has_C(){
        return first_C != first_G;
    }
    bool has_G(){
        return first_G != first_T;
    }
    bool has_T(){
        return first_T != last;
    }
};

struct left_ext_intervals {

    Interval TERM;
	Interval A;
	Interval C;
	Interval G;
	Interval T;

};

struct CArray {
	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;
};

class FmIndex{

public:
	FmIndex(){};

	/*
	 * constructor path of a STRING file containing the STRING in ASCII format
	 */
	FmIndex(string path, char TERM = '#'){

		this->TERM = TERM;

        STRING = SuccintString(path, TERM);
		auto r = STRING.parallel_rank(STRING.size());
        C = {
			r.TERM, 
			r.TERM+r.A, 
			r.TERM+r.A+r.C, 
			r.TERM+r.A+r.C+r.G
		};
    }
    /*
    * Input: suffix tree node N.
    * Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
    */
    left_ext_intervals LF(Interval interval){

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
    }

    uint64_t LF(uint64_t r){
        uint64_t f;
        // cout << STRING[r] << endl;
        switch (STRING[r])
        {
            case 'A': f=C.A+STRING.rank(r, 'A'); break;
            case 'C': f=C.C+STRING.rank(r, 'C'); break;
            case 'G': f=C.G+STRING.rank(r, 'G'); break;
            case 'T': f=C.T+STRING.rank(r, 'T'); break;
            default: //TERM, otherwise invalidated already.
                f=STRING.rank(r, TERM);
                break;
        }
        return f;
    }
    
    char get_character(uint64_t pos){
        return STRING[pos];
    }

    uint64_t size(){
        return STRING.size();
    }

	uint64_t seq_cnt(){
        return C.A;
    }

	// /*
	//  * functions for suffix tree navigation
	//  */

	Interval root(){

		return {
			0,
			C.A,
			C.C,
			C.G,
			C.T,
			STRING.size(),
			0
		};

	}

	/*
	debugging purpose
	*/
	string recover_text(int seq_no){
		uint64_t L = seq_no;
		uint64_t F = LF(L);

		string seq;
		while(F >= seq_cnt()){
			seq = get_character(L) + seq;
			L = F;
			F = LF(L);
		}
		// must be terminator symbol
		seq += get_character(L); 
		return seq;
	}

	/*
	debugging purpose
	*/
	vector<pair<uint64_t, string>> recover_suffix_array(int seq_no){
		uint64_t L = seq_no;
		uint64_t F = LF(L);

		string seq(1, TERM);
		vector<pair<uint64_t, string>> sa;
		while(F >= seq_cnt()){
			seq = get_character(L) + seq;
			sa.push_back(make_tuple(F, seq));
			L = F;
			F = LF(L);
		}
		// important: this can be misleading because F is randomly (in lexicographical order) chosen
		string term(1, TERM);
		sa.insert(sa.begin(), make_tuple(F, term));
		reverse(sa.begin(), sa.end());

		return sa;
	}
	/*
	 * depth = LCP inside the leaf.
	 */
	// sa_leaf first_leaf(){

	// 	return {{0, C.A}, 0};

	// }

	// void next_leaves(sa_leaf & L, vector<sa_leaf> & TMP_LEAVES, int & t, int min_n_children){

	// 	t = 0;

	// 	p_range ext = LF(L.rn);

	// 	sa_leaf A = {ext.A, L.depth+1};
	// 	sa_leaf C = {ext.C, L.depth+1};
	// 	sa_leaf G = {ext.G, L.depth+1};
	// 	sa_leaf T = {ext.T, L.depth+1};

	// 	if(leaf_size(A)>=min_n_children) TMP_LEAVES[t++] = A;
	// 	if(leaf_size(C)>=min_n_children) TMP_LEAVES[t++] = C;
	// 	if(leaf_size(G)>=min_n_children) TMP_LEAVES[t++] = G;
	// 	if(leaf_size(T)>=min_n_children) TMP_LEAVES[t++] = T;

	// 	std::sort( TMP_LEAVES.begin(), TMP_LEAVES.begin()+t, [ ]( const sa_leaf& lhs, const sa_leaf& rhs )
	// 	{
	// 		return leaf_size(lhs) < leaf_size(rhs);
	// 	});

	// }

	// void next_nodes(sa_node & N, vector<sa_node> & TMP_NODES, int & t){

	// 	p_node left_exts = LF(N);

	// 	sa_node A = left_exts.A;
	// 	sa_node C = left_exts.C;
	// 	sa_node G = left_exts.G;
	// 	sa_node T = left_exts.T;

	// 	t = 0;

	// 	if(number_of_children(A) >= 2) TMP_NODES[t++] = A;
	// 	if(number_of_children(C) >= 2) TMP_NODES[t++] = C;
	// 	if(number_of_children(G) >= 2) TMP_NODES[t++] = G;
	// 	if(number_of_children(T) >= 2) TMP_NODES[t++] = T;

	// 	//push right-maximal nodes on stack in decreasing size (i.e. interval length) order

	// 	std::sort( TMP_NODES.begin(), TMP_NODES.begin()+t, [ ]( const sa_node& lhs, const sa_node& rhs )
	// 	{
	// 		return node_size(lhs) < node_size(rhs);
	// 	});

	// }

private:

	char TERM = '#';

    CArray C;

	SuccintString STRING;
};

#endif /* FM_INDEX_INDEX_HPP_ */
