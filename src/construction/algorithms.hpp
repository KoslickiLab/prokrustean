
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include "algorithms.procedures.hpp"
#include <algorithm>

using namespace std;

Prokrustean build_prokrustean(FmIndex &fm_idx, uint64_t Lmin=1){
    // step1: collect representative suffix array
    SuffixArrayNode root = get_root(fm_idx);
    // vector<MaximalRepeat> repeats = {};
    vector<MaximalRepeatAnnotation> repeats = navigate_tree<MaximalRepeatAnnotation, get_rep_annot>(root, Lmin, fm_idx);

    // step2 get repr structure
    auto repr_annotation = ReprSuffixAnnotation();
    repr_annotation.initialize_rank(fm_idx.size(), repeats);
    repr_annotation.initialize_repr_sa(repeats);
    
    // step3 build structure
    // parallelize by sequences
    vector<MinCover> covers;
    for(uint64_t i; i < fm_idx.seq_cnt(); i++){
        SequenceAnnotation annot = get_seq_annot(i, fm_idx, repr_annotation);
        get_min_covers(annot);
    }
    return Prokrustean();
}

#endif