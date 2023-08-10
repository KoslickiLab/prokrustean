
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include "algorithms.procedures.hpp"
#include <algorithm>

using namespace std;

Prokrustean build_prokrustean(FmIndex &fm_idx, uint64_t Lmin=1){
    // step1: collect representative suffix array
    SuffixArrayInterval root = get_root(fm_idx);
    // vector<MaximalRepeat> repeats = {};
    vector<MaximalRepeat> repeats = navigate_tree<MaximalRepeat, get_repeat_and_repr_sa>(root, Lmin, fm_idx);

    // step2 get repr structure
    auto repr_annotation = ReprSuffixAnnotation(fm_idx.size());
    repr_annotation.set_repr_sa_rank(repeats);
    repr_annotation.set_repr_sa_annotation(repeats);

    // step3 build structure
    // parallelize by sequences
    vector<MinCover> covers;
    for(uint64_t i; i < fm_idx.seq_cnt(); i++){
        SeqAnnotation annot = get_annotation_structure(i, fm_idx, repr_annotation);
        get_min_covers(annot);
    }
    return Prokrustean();
}

#endif