

#include "rank.hpp"
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

struct CArray {
	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;
};

class FmIndex{

public:
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
	char TERM = '#';
    CArray C;
	SuccintString STRING;

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
};

#endif /* FM_INDEX_INDEX_HPP_ */
