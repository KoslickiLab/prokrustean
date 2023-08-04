#include "index.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_LOCATE_HPP_
#define FM_INDEX_LOCATE_HPP_

/*
debugging purpose
*/
string recover_text(FmIndex &fm_idx, int seq_no){
	uint64_t L = seq_no;
	uint64_t F = fm_idx.LF(L);

	string seq;
	while(F >= fm_idx.seq_cnt()){
		seq = fm_idx.get_character(L) + seq;
		L = F;
		F = fm_idx.LF(L);
	}
	// must be terminator symbol
	seq += fm_idx.get_character(L); 
	return seq;
}

/*
debugging purpose
*/
vector<pair<uint64_t, string>> recover_suffix_array(FmIndex &fm_idx, int seq_no){
	uint64_t L = seq_no;
	uint64_t F = fm_idx.LF(L);

	string seq(1, fm_idx.TERM);
	vector<pair<uint64_t, string>> sa;
	while(F >= fm_idx.seq_cnt()){
		seq = fm_idx.get_character(L) + seq;
		sa.push_back(make_tuple(F, seq));
		L = F;
		F = fm_idx.LF(L);
	}
	// important: this can be misleading because F is randomly (in lexicographical order) chosen
	string term(1, fm_idx.TERM);
	sa.insert(sa.begin(), make_tuple(F, term));
	reverse(sa.begin(), sa.end());

	return sa;
}

#endif /* FM_INDEX_LOCATE_HPP_ */
