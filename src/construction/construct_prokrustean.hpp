
#include "../fm_index/tree.hpp"
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include "repr_sa.hpp"
#include "repr_occ.hpp"
#include <algorithm>

using namespace std;

#ifndef CONSTRUCTION_PROKRUSTEAN_HPP_
#define CONSTRUCTION_PROKRUSTEAN_HPP_

void construct_repr_sa(left_ext_intervals intervals){
    
}

void construct_repr_occ(left_ext_intervals intervals){
    
}

void construct_prokrustean(FmIndex &fm_idx, int seq_no, ReprSABase &repr_sa, Prokrustean &prokrustean){
    // uint64_t L = seq_no;
    // uint64_t F = LF(L);

    // string seq(1, TERM);
    // vector<pair<uint64_t, string>> sa;
    // while(F >= seq_cnt()){
    //     seq = get_character(L) + seq;
    //     sa.push_back(make_tuple(F, seq));
    //     L = F;
    //     F = LF(L);
    // }
    // // important: this can be misleading because F is randomly (in lexicographical order) chosen
    // string term(1, TERM);
    // sa.insert(sa.begin(), make_tuple(F, term));
    // reverse(sa.begin(), sa.end());

    // return sa;
}

#endif