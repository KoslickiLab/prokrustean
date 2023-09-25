#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "util.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.procedures_new.hpp"
#include "../src/fm_index/tree_new.hpp"
#include "../src/application/kmers.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;


struct Step3RegionIdx{
    uint32_t pidx;
    uint16_t sidx;
    Step3RegionIdx(){}
    Step3RegionIdx(uint32_t pidx, uint16_t sidx): pidx(pidx), sidx(sidx)
    {}
};

struct Step3PositionAnnot{
    Pos pos;

    uint32_t stratum_cnt;
    
    StratifiedReprBased* repr_stratification_ref;
};

// struct RegionAnnot{
//     StratumId stratum_id;
//     bool is_primary;
// };

struct Step3SequenceAnnot {
    //
    SeqId id;

    // 
    uint64_t size;

    /* Due to bwt logic, the annotations are grouped by positions, 
    where positions are desc, stratums are asc by size */
    // vector<RegionAnnot> region_annots;

    //
    vector<Step3PositionAnnot> position_annots;

    //
    optional<string> sequence;
};

struct Step3ConsumingRef{
    bool this_pidx_is_exhausted;
    uint32_t pidx;
    uint16_t sidx;
    optional<uint32_t> post_pidx;

    void set_ref(uint32_t pidx, uint16_t sidx, optional<uint32_t> post_pidx=nullopt){
        this->this_pidx_is_exhausted=false;
        this->pidx=pidx;
        this->sidx=sidx;
        this->post_pidx=post_pidx;
    }
    void exhausted(){
        this->this_pidx_is_exhausted=true;
    }
};


struct Step3WorkSpace {
    //
    Step3SequenceAnnot seq_annot;
    //
    Step3ConsumingRef primary_ref;
    // (rgn, post-pid)
    vector<Step3ConsumingRef> consume_refs;

    void reset_by_annot(){
        consume_refs.clear();
        consume_refs=vector<Step3ConsumingRef>(seq_annot.position_annots.size());
        int pidx=0;
        for(auto &ref: consume_refs){
            if(pidx+1<seq_annot.position_annots.size()){
                ref.set_ref(pidx, 0, pidx+1);
            } else {
                ref.set_ref(pidx, 0);
            }
            pidx++;
        }
    }
    
    StratumId get_stratum_id(Step3RegionIdx &rgn){
        auto &pos_annot=seq_annot.position_annots[rgn.pidx];
        return get_region_annot(rgn.pidx, rgn.sidx).stratum_id;
    }

    Pos position(uint32_t pidx){
        return seq_annot.position_annots[pidx].pos;
    }

    StratifiedRegionAnnot get_region_annot(uint32_t pidx, uint16_t sidx){
        return seq_annot.position_annots[pidx].repr_stratification_ref->stratum_id_array[sidx];
    }

    uint16_t get_stratum_cnt(uint32_t pidx){
        return seq_annot.position_annots[pidx].stratum_cnt;
    }

    optional<Step3RegionIdx> child(Step3RegionIdx &rgn){
        if(rgn.sidx==0){
            return nullopt;
        }
        return Step3RegionIdx(rgn.pidx, rgn.sidx-1);
    }
    
    optional<Step3RegionIdx> parent(Step3RegionIdx &rgn){
        if(rgn.sidx+1>=get_stratum_cnt(rgn.pidx)){
            return nullopt;
        }
        return Step3RegionIdx(rgn.pidx, rgn.sidx+1);
    }

    optional<Step3RegionIdx> first_primary(){
        uint32_t pidx = seq_annot.position_annots.size()-1;
        uint16_t sidx = 0;
        bool is_primary=get_region_annot(pidx, sidx).is_primary;
        auto rgn = Step3RegionIdx(pidx, sidx);
        if(is_primary){
            return rgn;
        } else {
            return next_primary(rgn);
        }
    }
    optional<Step3RegionIdx> next_primary(Step3RegionIdx &rgn){
        uint32_t pidx = rgn.pidx;
        uint16_t sidx = rgn.sidx;
        while(true){
            if(sidx+1<seq_annot.position_annots[pidx].stratum_cnt){
                sidx++;
            } else if(pidx>0){
                pidx--;
                sidx=0;
            } else {
                return nullopt;
            }
            bool is_primary=get_region_annot(pidx, sidx).is_primary;
            if(is_primary){
                return Step3RegionIdx(pidx, sidx);
            }
        }
        return nullopt;
    }

    bool check_inclusion(Step3RegionIdx &a, Step3RegionIdx &b,  Prokrustean &prokrustean){
        /* a includes b? */
        auto a_right=position(a.pidx)+prokrustean.stratums[get_region_annot(a.pidx, a.sidx).stratum_id].size; 
        auto b_right=position(b.pidx)+prokrustean.stratums[get_region_annot(b.pidx, b.sidx).stratum_id].size;
        return a_right>=b_right;
    }

    void set_sequence_output(vector<Step3RegionIdx> &regions, Prokrustean &prokrustean){
        Region arr[regions.size()];
        for(int i=0; i<regions.size(); i++){
            arr[i].pos=seq_annot.position_annots[regions[i].pidx].pos;
            arr[i].stratum_id=get_stratum_id(regions[i]);
        }
        prokrustean.seqs[seq_annot.id].size=seq_annot.size;
        prokrustean.seqs[seq_annot.id].regions2=arr;
    }

    void set_stratum_output(Step3RegionIdx &primary, vector<Step3RegionIdx> &regions, Prokrustean &prokrustean){
        auto stratum_id = get_stratum_id(primary);
        
        Region arr[regions.size()];
        for(int i=0; i<regions.size(); i++){
            arr[i].pos=position(regions[i].pidx)-position(primary.pidx);
            arr[i].stratum_id=get_stratum_id(regions[i]);
            assert(arr[i].pos>=0);
        }
        prokrustean.stratums[stratum_id].regions2=arr;
    }
};

void build_prokrustean(Step3SequenceAnnot &seq_anot, Step3WorkSpace &work, Prokrustean &prokrustean){
    /*
    * - Input: Projected regions on a sequence, ordered and grouped by starting position, and ordered by its length,
    *         i.e. a list of position:(stratumId, isPrimary). With stratumId the size is accessed in constant time.
    *         isPrimary means the suffix starting from the region has the first index in the suffix array among those
    *         having the stratum as prefix.
    * - Output: Stratified region sets, where one set is stratified regions of the sequence and others are that of stratums.
    */
    // 1. initialize primaries
    work.reset_by_annot();
    vector<Step3RegionIdx> regions;
    optional<Step3RegionIdx> primary_rgn=work.first_primary();
    optional<uint32_t> primary_post_pidx;
    Step3RegionIdx post_rgn;
    optional<Step3RegionIdx> parent_rgn;
    while(primary_rgn.has_value()){
        regions.clear();
        primary_post_pidx=work.consume_refs[primary_rgn.value().pidx].post_pidx;
        //leftmost case
        auto child = work.child(primary_rgn.value());
        if(child.has_value()){
            regions.push_back(child.value());
        }
        //non-leftmost cases
        while(primary_post_pidx.has_value()){
            post_rgn.pidx=work.consume_refs[primary_post_pidx.value()].pidx;
            post_rgn.sidx=work.consume_refs[primary_post_pidx.value()].sidx;

            // (1) If post.rgn ∉ primary.rgn,
            if(!work.check_inclusion(primary_rgn.value(), post_rgn, prokrustean)){
                break;
            }

            // find the largest but included in the primary
            parent_rgn = work.parent(post_rgn);
            while(parent_rgn.has_value() && work.check_inclusion(primary_rgn.value(), parent_rgn.value(), prokrustean)){
                post_rgn=parent_rgn.value();
                parent_rgn=work.parent(post_rgn);
            }
            // add right region
            regions.push_back(post_rgn);

            // move consuming pointer
            if(parent_rgn.has_value()){
                work.consume_refs[primary_post_pidx.value()].pidx=parent_rgn.value().pidx;
                work.consume_refs[primary_post_pidx.value()].sidx=parent_rgn.value().sidx;
            } else {
                // all exhausted
                work.consume_refs[primary_post_pidx.value()].exhausted();
                work.consume_refs[primary_post_pidx.value()].post_pidx=work.consume_refs[post_rgn.pidx].post_pidx;
            }
            // go further to the right.
            primary_post_pidx=work.consume_refs[primary_post_pidx.value()].post_pidx;
        }

        work.set_stratum_output(primary_rgn.value(), regions, prokrustean);
        // move to next primary
        primary_rgn = work.next_primary(primary_rgn.value());
    }
    regions.clear();
    for(auto &ref: work.consume_refs){
        if(!ref.this_pidx_is_exhausted){
            /*important. there may be two non primary regions that somehow have inclusion relationship*/
            ref.sidx=work.get_stratum_cnt(ref.pidx)-1;
            regions.push_back(Step3RegionIdx(ref.pidx, ref.sidx));
        }
    }
    work.set_sequence_output(regions, prokrustean);
}

StratifiedReprBased* get_random_repr_stratification(int stra_cnt){
    auto seq_length=1000000;
    StratifiedReprBased* repr_stra = static_cast<StratifiedReprBased*>(std::malloc(sizeof(StratifiedReprBased)+stra_cnt*sizeof(tuple<StratumId, bool>)));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> id_dist(0, seq_length);
    for(int i=0; i<stra_cnt; i++){
        repr_stra->stratum_id_array[i].stratum_id=id_dist(gen);
        repr_stra->stratum_id_array[i].stratum_id=id_dist(gen)%2;
    }
    return repr_stra;
}


void test_basic_step3_operation(){
    auto pos_cnt=10;
    auto stra_cnt=3;
    Step3SequenceAnnot seq_annot;
    Step3WorkSpace work;
    work.seq_annot=seq_annot;
    work.seq_annot.position_annots=vector<Step3PositionAnnot>(pos_cnt);
    for(auto &pos_annot: work.seq_annot.position_annots){
        pos_annot.repr_stratification_ref=get_random_repr_stratification(stra_cnt);
        pos_annot.stratum_cnt=stra_cnt;
    }
    //reset function
    work.reset_by_annot();
    assert(work.consume_refs.size()==work.seq_annot.position_annots.size());
    
    //child function
    auto rgn = Step3RegionIdx(5, 2);
    auto child = work.child(rgn);
    assert(child.has_value());
    rgn = Step3RegionIdx(5, 0);
    child = work.child(rgn);
    assert(!child.has_value());

    //parent function
    rgn = Step3RegionIdx(5, stra_cnt-2);
    auto parent = work.parent(rgn);
    assert(parent.has_value());
    rgn = Step3RegionIdx(5, stra_cnt-1);
    parent = work.parent(rgn);
    assert(!parent.has_value());

    //first primary function
    //set all primary false
    for(auto &pos: work.seq_annot.position_annots){
        for(int i=0;i<pos.stratum_cnt; i++){
            pos.repr_stratification_ref->stratum_id_array[i].is_primary=false;
        }
    }

    work.seq_annot.position_annots[5].repr_stratification_ref->stratum_id_array[2].is_primary=true;
    auto primary=work.first_primary();
    assert(primary.has_value());
    assert(primary.value().pidx==5 && primary.value().sidx==2);
    // earlier 
    work.seq_annot.position_annots[6].repr_stratification_ref->stratum_id_array[2].is_primary=true;
    primary=work.first_primary();
    assert(primary.has_value());
    assert(primary.value().pidx==6 && primary.value().sidx==2);
    // earlier 
    work.seq_annot.position_annots[6].repr_stratification_ref->stratum_id_array[1].is_primary=true;
    primary=work.first_primary();
    assert(primary.has_value());
    assert(primary.value().pidx==6 && primary.value().sidx==1);
    
    // next 
    primary=work.first_primary();
    primary=work.next_primary(primary.value());
    assert(primary.value().pidx==6 && primary.value().sidx==2);
    primary=work.next_primary(primary.value());
    assert(primary.value().pidx==5 && primary.value().sidx==2);
    primary=work.next_primary(primary.value());
    assert(!primary.has_value());
    
    //check_inclusion
    Prokrustean prokrustean;
    prokrustean.seqs=vector<Sequence>(3);
    prokrustean.stratums=vector<Stratum>(3);
    prokrustean.stratums[0].size=5;
    prokrustean.stratums[1].size=7;
    prokrustean.stratums[2].size=10;
    auto rgn1=Step3RegionIdx(2, 1);
    auto rgn2=Step3RegionIdx(3, 1);
    auto rgn3=Step3RegionIdx(4, 1);
    work.seq_annot.position_annots[rgn1.pidx].pos=3;
    work.seq_annot.position_annots[rgn2.pidx].pos=4;
    work.seq_annot.position_annots[rgn3.pidx].pos=5;
    work.seq_annot.position_annots[rgn1.pidx].repr_stratification_ref->stratum_id_array[1].stratum_id=1;
    work.seq_annot.position_annots[rgn2.pidx].repr_stratification_ref->stratum_id_array[1].stratum_id=1;
    work.seq_annot.position_annots[rgn3.pidx].repr_stratification_ref->stratum_id_array[1].stratum_id=0;
    assert(!work.check_inclusion(rgn1, rgn2, prokrustean));
    assert(work.check_inclusion(rgn1, rgn3, prokrustean));
    assert(work.check_inclusion(rgn2, rgn3, prokrustean));
    // set output
    work.seq_annot.id=1;
    auto regions=vector<Step3RegionIdx>({rgn1, rgn2});
    work.set_sequence_output(regions, prokrustean);
    assert(prokrustean.seqs[1].regions2[0].stratum_id==work.get_stratum_id(rgn1));
    assert(prokrustean.seqs[1].regions2[1].stratum_id==work.get_stratum_id(rgn2));

    // set output stratum
    regions=vector<Step3RegionIdx>({rgn3});
    work.set_stratum_output(rgn1, regions, prokrustean);
    assert(prokrustean.stratums[work.get_stratum_id(rgn1)].regions2[0].stratum_id==work.get_stratum_id(rgn3));
}

void main_performance_build_prokrustean() {
    test_basic_step3_operation();
}
