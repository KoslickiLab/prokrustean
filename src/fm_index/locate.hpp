#include "index.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_LOCATE_HPP_
#define FM_INDEX_LOCATE_HPP_

/*
For more optimization:
https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/csa_sampling_strategy.hpp
text-order require another bitvector
*/
void sample_by_sa_order(fm_index &index, uint64_t s){
    //consider between two strategies
    //https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/csa_sampling_strategy.hpp
    int L_term = 0;
    int seq_no = 3;
    for(uint64_t i = 0; i < seq_no; i++){
        vector<pair<uint64_t, uint64_t>> sa_revpos;
        uint64_t L = i;
        uint64_t F = index.LF(L);
        while(F!=i){
            cout << "c "<< F << ":" << index.get_character(L) << endl;
            L = F;
            F = index.LF(L);
        }
    }
    vector<uint64_t> ssa(index.size()/s);
    // // assumption of lexicographical order
    // uint64_t L = 0;
    // for(uint64_t i = 0; i < index.size(); i++){
    //     char c = index.get_character(i);
    //     if (c == '#'){
    //         uint64_t L = i;
    //         uint64_t pos = i-1;
    //         uint64_t F = index.LF(L);
    //         while(F!=i){
    //             cout << "c "<< F << ":" <<  pos << ":" << index.get_character(F) << endl;
    //             pos--;
    //             L = F;
    //             F = index.LF(L);
    //         }
    //     }
    // }
    index.set_sampled_suffixes(ssa, s);
}

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

#endif