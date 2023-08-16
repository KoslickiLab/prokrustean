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

vector<string> recover_text(FmIndex &fm_idx){
	vector<string> seqs;
	for(int i=0; i<fm_idx.seq_cnt(); i++){
		seqs.push_back(recover_text(fm_idx, i));
	}
	return seqs;
}

/*
debugging purpose
index is position in the sequence,
first: sa index
second: actual suffix
*/
vector<pair<uint64_t, string>> recover_suffix_array(FmIndex &fm_idx, int seq_no, bool with_term=true){
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
	if(with_term){
		string term(1, fm_idx.TERM);
		sa.insert(sa.begin(), make_tuple(F, term));
	}
	reverse(sa.begin(), sa.end());

	return sa;
}

/*
debugging purpose
index is sa index
*/
vector<string> recover_suffix_array(FmIndex &fm_idx, bool with_term=true){
	vector<pair<uint64_t, string>> sa;
	for(int i=0; i<fm_idx.seq_cnt(); i++){
		for(auto pair: recover_suffix_array(fm_idx, i, with_term)){
			sa.push_back(pair);
		}
	}
	std::sort(sa.begin(), sa.end(), 
        [](tuple<int, string> const &t1, tuple<int, string> const &t2) {
            return get<0>(t1) < get<0>(t2); 
        });
	
	vector<string> suffixes;
	for(auto pair: sa){
		suffixes.push_back(pair.second);
	}
	return suffixes;
}

vector<int> recover_lcp(FmIndex &fm_idx){
	vector<int> lcp;
	auto sa = recover_suffix_array(fm_idx);
	for(int i=1; i<sa.size(); i++){
		int pre = 0;
		for (int j = 0; j < min(sa[i-1].size(), sa[i].size()); j++) {
			if (sa[i-1][j] != sa[i][j] || sa[i][j] == '#')
				break;
			pre++;
		}
		lcp.push_back(pre);
	}
	return lcp;
}

#endif /* FM_INDEX_LOCATE_HPP_ */
