#ifndef CONSTRUCTION_ALGO_STEP2_NEW_HPP_
#define CONSTRUCTION_ALGO_STEP2_NEW_HPP_
#include <algorithm>
#include <stack>
#include <tuple>
#include <unordered_set>
#include "algorithms.step1_project_stratums.hpp"

using namespace std;

struct RegionIdx{
    uint32_t pidx;
    uint16_t sidx;
    RegionIdx(){}
    RegionIdx(uint32_t pidx, uint16_t sidx): pidx(pidx), sidx(sidx)
    {}
};

struct RegionAnnotation{
    StratumId stratum_id;
    StratumSize stratum_size;
    bool is_primary; 
    RegionAnnotation()
    {}
    RegionAnnotation(StratumId stratum_id, StratumSize stratum_size, bool is_primary):stratum_id(stratum_id), stratum_size(stratum_size), is_primary(is_primary)
    {}
};

struct PositionAnnotation{
    Pos pos;
    
    vector<RegionAnnotation> regions;
};


struct SequenceAnnotation {
    //
    SeqId seq_id;

    uint32_t seq_size;
    //
    vector<PositionAnnotation> position_annots;
};

struct ConsumingReference{
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


struct StratificationWorkSpace {
    //
    SequenceAnnotation seq_annot;
    //
    ConsumingReference primary_ref;
    // (rgn, post-pid)
    vector<ConsumingReference> consume_refs;
    // //
    // std::unordered_map<Pos, int> pos_stats;
    // //
    // std::unordered_map<Pos, int> pos_inverse;
    // //
    // std::vector<Pos> pos_vector;

    // void update_contexts_by_seq(SeqId seq_id, vector<ProjectedStratifiedRegion> &projected_regions, vector<StratumSize> &stratum_sizes){
    //     this->seq_annot.seq_id=seq_id;

    //     this->pos_stats.clear();
    //     this->pos_inverse.clear();
    //     this->pos_vector.clear();
    //     // collect positions
    //     for(auto &region: projected_regions){
    //         if(pos_stats.count(region.pos)==0){
    //             this->pos_stats[region.pos]=1;
    //             this->pos_vector.push_back(region.pos);
    //         } else {
    //             this->pos_stats[region.pos]++;
    //         }
    //     }
    //     // reuse the space as efficient as possible
    //     this->seq_annot.position_annots.resize(pos_stats.size());
        
    //     // set positions and sort
    //     sort(this->pos_vector.begin(), this->pos_vector.end());
    //     int i=0;
    //     for(auto pos: pos_vector){
    //         this->seq_annot.position_annots[i].pos=pos;
    //         this->seq_annot.position_annots[i].regions.resize(this->pos_stats[pos]);
    //         this->pos_inverse[pos]=i;
    //         i++;
    //     }
    //     // set regions and sort
    //     int region_i;
    //     for(auto &region: projected_regions){
    //         this->pos_stats[region.pos]--;
    //         region_i=this->pos_stats[region.pos];
    //         auto& region_annot=this->seq_annot.position_annots[this->pos_inverse[region.pos]].regions[region_i];
    //         region_annot.stratum_id=region.stratum_id;
    //         region_annot.stratum_size=stratum_sizes[region.stratum_id];
    //         region_annot.is_primary=region.is_primary;
    //     }
    //     for(auto &annot: this->seq_annot.position_annots){
    //         sort(annot.regions.begin(), annot.regions.end(), [](const RegionAnnotation& lhs, const RegionAnnotation& rhs) { 
    //             return lhs.stratum_size < rhs.stratum_size;
    //         });
    //     }
        
    //     // clear raw data
    //     projected_regions.clear();
    //     projected_regions.shrink_to_fit();

    //     // reset references
    //     consume_refs.clear();
    //     consume_refs.resize(seq_annot.position_annots.size());
    //     int pidx=0;
    //     for(auto &ref: consume_refs){
    //         if(pidx+1<seq_annot.position_annots.size()){
    //             ref.set_ref(pidx, 0, pidx+1);
    //         } else {
    //             ref.set_ref(pidx, 0);
    //         }
    //         pidx++;
    //     }
    // }
    
    void update_contexts_for_seq(SeqId seq_id, FmIndex &fm_index, StratumProjectionOutput &output, vector<StratumSize> &stratum_sizes){
        this->seq_annot.seq_id=seq_id;

        // this->pos_stats.clear();
        // this->pos_inverse.clear();
        // this->pos_vector.clear();
        this->seq_annot.position_annots.clear();
        
        uint64_t seq_cnt=fm_index.seq_cnt();
        SuffixArrayIdx L = seq_id;
	    SuffixArrayIdx F = fm_index.LF(L);
        uint64_t idx;
        Pos reverse_pos=0;
        vector<Pos> reverse_positions;
        while(F >= seq_cnt){
            optional<vector<ProjectedStratifiedRegion>*> regions=output.fetch(F);
            if(regions.has_value()){
                auto annot = PositionAnnotation();
                annot.pos=reverse_pos;
                // vector<ProjectedStratifiedRegion> aa =regions.value();
                for(auto &r: *regions.value()){
                    annot.regions.push_back(RegionAnnotation(r.stratum_id, stratum_sizes[r.stratum_id], r.is_primary));
                }
                output.dispose(F);
                this->seq_annot.position_annots.push_back(annot);
            }

            L = F;
            F = fm_index.LF(L);
            reverse_pos++;
        }

        // reuse the space as efficient as possible
        this->seq_annot.seq_size=reverse_pos;
        for(auto &annot: this->seq_annot.position_annots){
            annot.pos=this->seq_annot.seq_size-annot.pos-1;
            sort(annot.regions.begin(), annot.regions.end(), [](const RegionAnnotation& lhs, const RegionAnnotation& rhs) { 
                return lhs.stratum_size < rhs.stratum_size;
            });
        }
        sort(this->seq_annot.position_annots.begin(), this->seq_annot.position_annots.end(), [](const PositionAnnotation& lhs, const PositionAnnotation& rhs) { 
            return lhs.pos < rhs.pos;
        });

        // reset references
        consume_refs.clear();
        consume_refs.resize(seq_annot.position_annots.size());
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
    
    StratumId get_stratum_id(RegionIdx &rgn){
        auto &pos_annot=seq_annot.position_annots[rgn.pidx];
        return get_region_annot(rgn.pidx, rgn.sidx)->stratum_id;
    }

    Pos position(uint32_t pidx){
        return seq_annot.position_annots[pidx].pos;
    }

    RegionAnnotation* get_region_annot(uint32_t pidx, uint16_t sidx){
        return &seq_annot.position_annots[pidx].regions[sidx];
    }

    uint16_t get_stratum_cnt(uint32_t pidx){
        return seq_annot.position_annots[pidx].regions.size();
    }

    optional<RegionIdx> child(RegionIdx &rgn){
        if(rgn.sidx==0){
            return nullopt;
        }
        return RegionIdx(rgn.pidx, rgn.sidx-1);
    }
    
    optional<RegionIdx> parent(RegionIdx &rgn){
        if(rgn.sidx+1>=get_stratum_cnt(rgn.pidx)){
            return nullopt;
        }
        return RegionIdx(rgn.pidx, rgn.sidx+1);
    }

    optional<RegionIdx> first_primary(){
        uint32_t pidx = seq_annot.position_annots.size()-1;
        uint16_t sidx = 0;
        bool is_primary=get_region_annot(pidx, sidx)->is_primary;
        auto rgn = RegionIdx(pidx, sidx);
        if(is_primary){
            return rgn;
        } else {
            return next_primary(rgn);
        }
    }
    optional<RegionIdx> next_primary(RegionIdx &rgn){
        uint32_t pidx = rgn.pidx;
        uint16_t sidx = rgn.sidx;
        while(true){
            if(sidx+1<seq_annot.position_annots[pidx].regions.size()){
                sidx++;
            } else if(pidx>0){
                pidx--;
                sidx=0;
            } else {
                return nullopt;
            }
            if(get_region_annot(pidx, sidx)->is_primary){
                return RegionIdx(pidx, sidx);
            }
        }
        return nullopt;
    }

    bool check_inclusion(RegionIdx &a, RegionIdx &b){
        /* a includes b? */
        auto a_right=position(a.pidx)+get_region_annot(a.pidx, a.sidx)->stratum_size; 
        auto b_right=position(b.pidx)+get_region_annot(b.pidx, b.sidx)->stratum_size;
        return a_right>=b_right;
    }

    void set_sequence_output(Prokrustean &prokrustean, vector<RegionIdx>* regions=nullptr){
        if(regions!=nullptr){
            Region arr[(*regions).size()];
            for(int i=0; i<(*regions).size(); i++){
                arr[i].pos=seq_annot.position_annots[(*regions)[i].pidx].pos;
                arr[i].stratum_id=get_stratum_id((*regions)[i]);
            }
            prokrustean.sequences__region[this->seq_annot.seq_id]=arr;
            prokrustean.sequences__region_cnt[this->seq_annot.seq_id]=(*regions).size();
        }
        
        prokrustean.sequences__size[this->seq_annot.seq_id]=this->seq_annot.seq_size;
    }

    void set_stratum_output(RegionIdx &primary, vector<RegionIdx> &regions, Prokrustean &prokrustean){
        auto stratum_id = get_stratum_id(primary);
        
        Region arr[regions.size()];
        for(int i=0; i<regions.size(); i++){
            arr[i].pos=position(regions[i].pidx)-position(primary.pidx);
            arr[i].stratum_id=get_stratum_id(regions[i]);
        }
        prokrustean.stratums__region[stratum_id]=arr;
        prokrustean.stratums__region_cnt[stratum_id]=regions.size();
    }
};

void build_prokrustean(StratificationWorkSpace &workspace, Prokrustean &prokrustean){
    /*
    * - Input: Projected regions on a sequence, ordered and grouped by starting position, and ordered by its length,
    *         i.e. a list of position:(stratumId, isPrimary). With stratumId the size is accessed in constant time.
    *         isPrimary means the suffix starting from the region has the first index in the suffix array among those
    *         having the stratum as prefix.
    * - Output: Stratified region sets, where one set is stratified regions of the sequence and others are that of stratums.
    */
   
    if(workspace.seq_annot.position_annots.size()==0){
        workspace.set_sequence_output(prokrustean);
        return;
    }
    // 1. initialize primaries
    vector<RegionIdx> regions;
    optional<RegionIdx> primary_rgn=workspace.first_primary();
    optional<uint32_t> primary_post_pidx;
    RegionIdx post_rgn;
    optional<RegionIdx> parent_rgn;
    while(primary_rgn.has_value()){
        regions.clear();
        //leftmost case
        auto child = workspace.child(primary_rgn.value());
        if(child.has_value()){
            regions.push_back(child.value());
        }
        //non-leftmost cases
        primary_post_pidx=workspace.consume_refs[primary_rgn.value().pidx].post_pidx;
        while(primary_post_pidx.has_value()){
            post_rgn.pidx=workspace.consume_refs[primary_post_pidx.value()].pidx;
            post_rgn.sidx=workspace.consume_refs[primary_post_pidx.value()].sidx;

            // (1) If post.rgn ∉ primary.rgn,
            if(!workspace.check_inclusion(primary_rgn.value(), post_rgn)){
                break;
            }

            // find the largest but included in the primary
            parent_rgn = workspace.parent(post_rgn);
            while(parent_rgn.has_value() && workspace.check_inclusion(primary_rgn.value(), parent_rgn.value())){
                post_rgn=parent_rgn.value();
                parent_rgn=workspace.parent(post_rgn);
            }
            // add right region
            regions.push_back(post_rgn);

            // move consuming pointer
            if(parent_rgn.has_value()){
                workspace.consume_refs[primary_post_pidx.value()].pidx=parent_rgn.value().pidx;
                workspace.consume_refs[primary_post_pidx.value()].sidx=parent_rgn.value().sidx;
            } else {
                // all exhausted
                workspace.consume_refs[primary_post_pidx.value()].exhausted();
                workspace.consume_refs[primary_post_pidx.value()].post_pidx=workspace.consume_refs[post_rgn.pidx].post_pidx;
            }
            // go further to the right.
            primary_post_pidx=workspace.consume_refs[primary_post_pidx.value()].post_pidx;
        }
        workspace.set_stratum_output(primary_rgn.value(), regions, prokrustean);
        // move to next primary
        primary_rgn = workspace.next_primary(primary_rgn.value());
    }
    regions.clear();
    for(auto &ref: workspace.consume_refs){
        if(!ref.this_pidx_is_exhausted){
            /*important. there may be two non primary regions that somehow have inclusion relationship*/
            ref.sidx=workspace.get_stratum_cnt(ref.pidx)-1;
            regions.push_back(RegionIdx(ref.pidx, ref.sidx));
        }
    }
    workspace.set_sequence_output(prokrustean, &regions);
}


#endif