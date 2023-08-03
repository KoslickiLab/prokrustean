#include "index.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_LOCATE_HPP_
#define FM_INDEX_LOCATE_HPP_


// void sample_by_sa_order(fm_index &index, uint64_t s){
//     //consider between two strategies
//     //https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/csa_sampling_strategy.hpp
//     vector<uint64_t> ssa(index.size()/s);
//     // assumption of lexicographical order
//     ssa[0] = index.size() -1;

//     uint64_t L = 0;
//     for(uint64_t i = (index.size()-1); i > 0; i--){
//         uint64_t F = index.LF(L);
//         // uint64_t sa = F-1;
//         uint64_t pos = i-1;
//         if(F % s == 1){
//             ssa[F/s]=pos;
//         }
//         cout << "c "<< F << ":" << pos << ":" << index.get_character(L) << endl;
//         L = F;
//     }
//     uint64_t F = index.LF(L);
//     cout << "c last"<< F << ":"  << ":" << index.get_character(L) << endl;
//     for(uint64_t i = 0;i< ssa.size(); i++){
//         cout << "sample "<< i*s+1 << ":" << ssa[i] << endl;
//     }
//     index.set_sampled_suffixes(ssa, s);
// }

// void sample_by_sa_order(fm_index &index, uint64_t s){
//     //consider between two strategies
//     //https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/csa_sampling_strategy.hpp
//     vector<uint64_t> ssa(index.size()/s);
//     // assumption of lexicographical order
//     ssa[0] = index.size() -1;

//     uint64_t L = 0;
//     for(uint64_t i = (index.size()-1); i > 0; i--){
//         uint64_t F = index.LF(L);
//         // uint64_t sa = F-1;
//         uint64_t pos = i-1;
//         if(F % s == 1){
//             ssa[F/s]=pos;
//         }
//         cout << "c "<< F << ":" << pos << ":" << index.get_character(L) << endl;
//         L = F;
//     }
//     uint64_t F = index.LF(L);
//     cout << "c last"<< F << ":"  << ":" << index.get_character(L) << endl;
//     for(uint64_t i = 0;i< ssa.size(); i++){
//         cout << "sample "<< i*s+1 << ":" << ssa[i] << endl;
//     }
//     index.set_sampled_suffixes(ssa, s);
// }

/*
debugging purpose
*/
string recover_text(FmIndex &index, int seq_no){
    uint64_t L = seq_no;
    uint64_t F = index.LF(L);

    string seq;
    while(F >= index.seq_cnt()){
        seq = index.get_character(L) + seq;
        L = F;
        F = index.LF(L);
    }
    // must be terminator symbol
    seq += index.get_character(L); 
    return seq;
}

/*
debugging purpose
*/
vector<pair<uint64_t, string>> recover_suffix_array(FmIndex &index, int seq_no){
    uint64_t L = seq_no;
    uint64_t F = index.LF(L);

    string term = "#";
    string seq = term;
    vector<pair<uint64_t, string>> sa;
    while(F >= index.seq_cnt()){
        seq = index.get_character(L) + seq;
        sa.push_back(make_tuple(F, seq));
        L = F;
        F = index.LF(L);
    }
    // important: this can be misleading because F is randomly (in lexicographical order) chosen
    sa.insert(sa.begin(), make_tuple(F, term));
    reverse(sa.begin(), sa.end());

    return sa;
}

#endif