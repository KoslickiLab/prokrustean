#ifndef CONSTRUCTION_ALGO_STEP2_NEW_HPP_
#define CONSTRUCTION_ALGO_STEP2_NEW_HPP_
#include <algorithm>
#include <stack>
#include <tuple>
#include <unordered_set>
#include "algorithms.stage1_projection.hpp"

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

    void print(){
        cout << "pos: " << pos << ", ";
        for(auto &r: regions){
            cout << r.stratum_id << "("<< r.stratum_size << ") ";
        }
        cout << endl;
    }
};


struct SequenceAnnotation {
    //
    SeqId seq_id;

    uint32_t seq_size;
    //
    vector<PositionAnnotation> position_annots;

    void print(){
        cout << "seq annot: " << seq_id << endl;
        for(auto &pos_annot: position_annots){
            pos_annot.print();
        }
    }
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

    // shared memory
    vector<bool> work_is_primaries;
    vector<StratumId> work_stratum_ids;
    
    void update_contexts_for_seq(SeqId seq_id, FmIndex &fm_index, StratumProjectionWorkspace &output, Prokrustean &prokrustean){
        this->seq_annot.seq_id=seq_id;
        uint64_t pos_idx=0;
        
        uint64_t seq_cnt=fm_index.seq_cnt();
        SuffixArrayIdx L = seq_id;
	    SuffixArrayIdx F = fm_index.LF(L);
        uint64_t idx;
        Pos reverse_pos=0;
        int cnt;
        while(F >= seq_cnt){
            if(output.fetch(F, cnt, work_stratum_ids, work_is_primaries)){
                if(pos_idx>=this->seq_annot.position_annots.size()){
                    this->seq_annot.position_annots.push_back(PositionAnnotation());
                }
                auto &annot = this->seq_annot.position_annots[pos_idx];
                annot.regions.clear();
                pos_idx++;
                annot.pos=reverse_pos;
                for(int i=0; i<cnt; i++){
                    annot.regions.push_back(RegionAnnotation(work_stratum_ids[i], prokrustean.get_stratum_size(work_stratum_ids[i]), work_is_primaries[i]));
                }
            }
            L = F;
            F = fm_index.LF(L);
            reverse_pos++;
        }
        // reuse the space as efficient as possible
        if(pos_idx!=this->seq_annot.position_annots.size()){
            // since this variable is reused, make it fit.
            this->seq_annot.position_annots.resize(pos_idx);
        }
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

    // optional<RegionIdx> first_primary(){
    //     uint32_t pidx = seq_annot.position_annots.size()-1;
    //     uint16_t sidx = 0;
    //     bool is_primary=get_region_annot(pidx, sidx)->is_primary;
    //     auto rgn = RegionIdx(pidx, sidx);
    //     if(is_primary){
    //         return rgn;
    //     } else {
    //         return next_primary(rgn);
    //     }
    // }
    // optional<RegionIdx> next_primary(RegionIdx &rgn){
    //     uint32_t pidx = rgn.pidx;
    //     uint16_t sidx = rgn.sidx;
    //     while(true){
    //         if(sidx+1<seq_annot.position_annots[pidx].regions.size()){
    //             sidx++;
    //         } else if(pidx>0){
    //             pidx--;
    //             sidx=0;
    //         } else {
    //             return nullopt;
    //         }
    //         if(get_region_annot(pidx, sidx)->is_primary){
    //             return RegionIdx(pidx, sidx);
    //         }
    //     }
    //     return nullopt;
    // }

    optional<RegionIdx> first_region(){
        uint32_t pidx = seq_annot.position_annots.size()-1;
        uint16_t sidx = 0;
        return RegionIdx(pidx, sidx);
    }

    optional<RegionIdx> next_region(RegionIdx &rgn){
        uint32_t pidx = rgn.pidx;
        uint16_t sidx = rgn.sidx;
        if(sidx+1<seq_annot.position_annots[pidx].regions.size()){
            return RegionIdx(pidx, sidx+1);
        } else if(pidx>0){
            return RegionIdx(pidx-1, 0);
        } else {
            return nullopt;
        }
    }

    bool check_inclusion(RegionIdx &a, RegionIdx &b){
        /* a includes b? */
        auto a_right=position(a.pidx)+get_region_annot(a.pidx, a.sidx)->stratum_size; 
        auto b_right=position(b.pidx)+get_region_annot(b.pidx, b.sidx)->stratum_size;
        return a_right>=b_right;
    }

    void set_sequence_output(Prokrustean &prokrustean, vector<RegionIdx>* regions=nullptr){
        if(regions==nullptr || (*regions).size()==0){
            prokrustean.set_seq_regions(this->seq_annot.seq_id, this->seq_annot.seq_size, nullptr, 0);
            return;
        }
        sort((*regions).begin(),(*regions).end(), [](const RegionIdx& lhs, const RegionIdx& rhs) { return lhs.pidx < rhs.pidx;});
        // StratifiedData arr[(*regions).size()];
        StratifiedData* arr= new StratifiedData[(*regions).size()];
        for(int i=0; i<(*regions).size(); i++){
            arr[i].pos=seq_annot.position_annots[(*regions)[i].pidx].pos;
            arr[i].stratum_id=get_stratum_id((*regions)[i]);
        }
        prokrustean.set_seq_regions(this->seq_annot.seq_id, this->seq_annot.seq_size, arr, regions->size());
    }

    void set_stratum_output(RegionIdx &primary, Prokrustean &prokrustean, vector<RegionIdx>* regions=nullptr){
        if(regions==nullptr || (*regions).size()==0){
            return;
        }
        sort((*regions).begin(),(*regions).end(), [](const RegionIdx& lhs, const RegionIdx& rhs) { return lhs.pidx < rhs.pidx;});
        // StratifiedData arr[(*regions).size()];
        StratifiedData* arr= new StratifiedData[(*regions).size()];
        for(int i=0; i<(*regions).size(); i++){
            arr[i].pos=position((*regions)[i].pidx)-position(primary.pidx);
            arr[i].stratum_id=get_stratum_id((*regions)[i]);
        }
        auto stratum_id = get_stratum_id(primary);
        prokrustean.set_stratum_regions(stratum_id, arr, regions->size());
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
    // cout << "------ seq: " << workspace.seq_annot.seq_id << " --------" << endl;
    // for(auto annot: workspace.seq_annot.position_annots){
    //     annot.print();
    // }
    if(workspace.seq_annot.position_annots.size()==0){
        workspace.set_sequence_output(prokrustean);
        return;
    }
    // 1. initialize primaries
    vector<RegionIdx> regions;
    optional<RegionIdx> focused_rgn=workspace.first_region();
    optional<uint32_t> focused_post_pidx;
    RegionIdx post_rgn;
    optional<RegionIdx> parent_of_post_rgn;
    while(focused_rgn.has_value()){
        auto is_primary=workspace.get_region_annot(focused_rgn.value().pidx, focused_rgn.value().sidx)->is_primary;
        if(is_primary){
            regions.clear();
            //leftmost case
            auto child = workspace.child(focused_rgn.value());
            if(child.has_value()){
                regions.push_back(child.value());
            }
        }
        
        //non-leftmost cases
        focused_post_pidx=workspace.consume_refs[focused_rgn.value().pidx].post_pidx;
        while(focused_post_pidx.has_value()){
            post_rgn.pidx=workspace.consume_refs[focused_post_pidx.value()].pidx;
            post_rgn.sidx=workspace.consume_refs[focused_post_pidx.value()].sidx;

            // (1) If post.rgn ∉ primary.rgn,
            if(!workspace.check_inclusion(focused_rgn.value(), post_rgn)){
                break;
            }
            // find the largest but included in the primary
            parent_of_post_rgn = workspace.parent(post_rgn);
            while(parent_of_post_rgn.has_value() && workspace.check_inclusion(focused_rgn.value(), parent_of_post_rgn.value())){
                post_rgn=parent_of_post_rgn.value();
                parent_of_post_rgn=workspace.parent(post_rgn);
            }
            if(is_primary){
                // add right region
                regions.push_back(post_rgn);
            }

            // move consuming pointer
            if(parent_of_post_rgn.has_value()){
                // workspace.consume_refs[primary_post_pidx.value()].pidx=parent_of_post_rgn.value().pidx;
                workspace.consume_refs[post_rgn.pidx].sidx=parent_of_post_rgn.value().sidx;
            } else {
                // all exhausted
                // workspace.consume_refs[primary_post_pidx.value()].exhausted();
                workspace.consume_refs[post_rgn.pidx].exhausted();
                // workspace.consume_refs[post_rgn.pidx].post_pidx=workspace.consume_refs[post_rgn.pidx].post_pidx;
                workspace.consume_refs[focused_rgn.value().pidx].post_pidx=workspace.consume_refs[post_rgn.pidx].post_pidx;
            }
            // go further to the right.
            focused_post_pidx=workspace.consume_refs[post_rgn.pidx].post_pidx;
        }
        if(is_primary){
            workspace.set_stratum_output(focused_rgn.value(), prokrustean, &regions);
        }
        // move to next region
        focused_rgn = workspace.next_region(focused_rgn.value());
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