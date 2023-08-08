
#ifndef CONSTRUCTION_PROKRUSTEAN_HPP_
#define CONSTRUCTION_PROKRUSTEAN_HPP_
#include "../fm_index/tree.hpp"
#include "../fm_index/index.hpp"
#include "../prokrustean.hpp"
#include <algorithm>

using namespace std;

struct MaxRepeat {
public:
    // 
    uint64_t size;
    vector<uint64_t> repr_sa_indexes;
};

struct SeqAnnot {
public:
    // 
    uint64_t size;
    vector<tuple<Pos, vector<RepId>>> annotations;
};

class AbsReprSuffPositions {
public:
    // 
    virtual uint64_t count() = 0;

    // 
    virtual void set(uint64_t sa_idx) = 0;

    // 
    virtual bool exists(uint64_t sa_idx) = 0;

    //
    virtual uint64_t rank(uint64_t sa_idx) = 0;
};

class AbsReprSuff {
public:
    virtual void initialize(uint64_t repr_sa_cnt) = 0;

    //
    virtual void set(uint64_t sa_index, RepId rep_id) = 0;

    // 
    virtual vector<RepId> get(uint64_t sa_index, bool remove=true) = 0;
};


class AbsProkrusteanProcedures {
public:
    //
    virtual vector<MaxRepeat> collect_repeats_from_subtree(Interval &curr_interval, FmIndex &fm_idx, AbsReprSuffPositions &repr_sa_pos) = 0;

    //
    virtual void index_representative_suffix_array(vector<MaxRepeat> repeats, AbsReprSuff &repr_sa) = 0;

    //
    virtual SeqAnnot get_annotation_structure(SeqId id, FmIndex &fm_idx, AbsReprSuffPositions &repr_sa_pos, AbsReprSuff &repr_sa) = 0;

    // 
    virtual vector<MinCover> get_min_covers(SeqAnnot annotation) = 0;
};

Prokrustean build_prokrustean(FmIndex &fm_idx, AbsProkrusteanProcedures &procedures, AbsReprSuffPositions &repr_sa_pos, AbsReprSuff &repr_sa){
    // step1: collect representative suffix array
    Interval root = get_root(fm_idx);
    // parallelize by subtree
    vector<MaxRepeat> repeats = procedures.collect_repeats_from_subtree(root, fm_idx, repr_sa_pos);

    // step2 flip structure
    // parallize by repeats partitions
    procedures.index_representative_suffix_array(repeats, repr_sa);

    // step3 build structure
    // parallelize by sequences
    vector<MinCover> covers;
    for(uint64_t i; i < fm_idx.seq_cnt(); i++){
        SeqAnnot annot = procedures.get_annotation_structure(i, fm_idx, repr_sa_pos, repr_sa);
        procedures.get_min_covers(annot);
    }
    return Prokrustean();
}

#endif