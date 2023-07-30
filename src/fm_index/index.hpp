

#include "rank.hpp"
#include "tree.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

/*
 * representation of a right-maximal substring (SA node) as a list of STRING intervals
 */
struct interval {

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

    interval TERM;
	interval A;
	interval C;
	interval G;
	interval T;

};

class fm_index{

public:
	fm_index(){};

	/*
	 * constructor path of a STRING file containing the STRING in ASCII format
	 */
	fm_index(string path, char TERM = '#'){

		this->TERM = TERM;

        STRING = succint_string(path, TERM);
        
		//build F column
		for(uint64_t i=0;i<STRING.size();++i){

			// assert(STRING[i] < 256);

			F_A += (STRING[i]==TERM);
			F_C += (STRING[i]=='A');
			F_G += (STRING[i]=='C');
			F_T += (STRING[i]=='G');

		}

		F_C += F_A;
		F_G += F_C;
		F_T += F_G;

	}

    /*
    * Input: suffix tree node N.
    * Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
    */
    left_ext_intervals LF(interval interval){

        p_rank rank = STRING.parallel_rank(interval.first_TERM);
        p_rank before_TERM = rank;

        if(interval.first_A != interval.first_TERM) rank = STRING.parallel_rank(interval.first_A);
        p_rank before_A = rank;

        if(interval.first_C != interval.first_A) rank = STRING.parallel_rank(interval.first_C);
        p_rank before_C = rank;

        if(interval.first_G != interval.first_C) rank = STRING.parallel_rank(interval.first_G);
        p_rank before_G = rank;

        if(interval.first_T != interval.first_G) rank = STRING.parallel_rank(interval.first_T);
        p_rank before_T = rank;

        if(interval.last != interval.first_T) rank = STRING.parallel_rank(interval.last);
        p_rank before_end = rank;

        return {
            {   before_TERM.TERM,    before_A.TERM,    before_C.TERM,    before_G.TERM,    before_T.TERM,    before_end.TERM, interval.depth+1},
            {F_A + before_TERM.A, F_A + before_A.A, F_A + before_C.A, F_A + before_G.A, F_A + before_T.A, F_A + before_end.A, interval.depth+1},
            {F_C + before_TERM.C, F_C + before_A.C, F_C + before_C.C, F_C + before_G.C, F_C + before_T.C, F_C + before_end.C, interval.depth+1},
            {F_G + before_TERM.G, F_G + before_A.G, F_G + before_C.G, F_G + before_G.G, F_G + before_T.G, F_G + before_end.G, interval.depth+1},
            {F_T + before_TERM.T, F_T + before_A.T, F_T + before_C.T, F_T + before_G.T, F_T + before_T.T, F_T + before_end.T, interval.depth+1}
        };

    }

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;
        uint64_t n = STRING.size();
		out.write((char*)&n,sizeof(n));
		out.write((char*)&F_A,sizeof(uint64_t));
		out.write((char*)&F_C,sizeof(uint64_t));
		out.write((char*)&F_G,sizeof(uint64_t));
		out.write((char*)&F_T,sizeof(uint64_t));

		w_bytes += sizeof(n) + sizeof(uint64_t)*4;

		w_bytes += STRING.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	// void load(std::istream& in) {

	// 	in.read((char*)&n,sizeof(n));

	// 	in.read((char*)F_A,sizeof(uint64_t));
	// 	in.read((char*)F_C,sizeof(uint64_t));
	// 	in.read((char*)F_G,sizeof(uint64_t));
	// 	in.read((char*)F_T,sizeof(uint64_t));

	// 	STRING.load(in);

	// }

	// void save_to_file(string path){

	// 	std::ofstream out(path);
	// 	serialize(out);
	// 	out.close();

	// }

	/*
	 * path = path of an index file
	 */
	// void load_from_file(string path){

	// 	std::ifstream in(path);
	// 	load(in);
	// 	in.close();

	// }


	// /*
	//  * functions for suffix tree navigation
	//  */

	interval root(){

		return {
			0,
			F_A,
			F_C,
			F_G,
			F_T,
			STRING.size(),
			0
		};

	}

	/*
	 * depth = LCP inside the leaf.
	 */
	// sa_leaf first_leaf(){

	// 	return {{0, F_A}, 0};

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

	uint64_t F_A=0; //F array
	uint64_t F_C=0; //F array
	uint64_t F_G=0; //F array
	uint64_t F_T=0; //F array

	//vector<uint64_t> F;
	succint_string STRING;

};

#endif /* FM_INDEX_INDEX_HPP_ */
