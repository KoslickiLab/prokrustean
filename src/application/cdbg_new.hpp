#ifndef APPLICATION_CDBG_NEW_HPP_
#define APPLICATION_CDBG_NEW_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"
#include "../util/string.access.hpp"
#include "../util/data.store.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;
/* 
* (0) Data Structure
* - Prokrustean (left/right counts)
* - Sequences
* - Stratum -> ancestor map
* - Stratum/Sequence -> descendant map
* (1) Ancestor setup
* - Scan Stra -> count right reflection, assign ancestor vector index
* - fetch structure: exists bitvector, rank -> position (+right count) -> ancestors
* (2) Descendant setup
* - Seq -> Descendant
* - Stratum -> Descendant
* (3) Above two structure should allow counting maximal unitig 
*/

typedef uint32_t AncestorIdx; // 

struct Ancestor {
    uint8_t is_seq;
    StratumOrSeqId id;
    CoveringRegionIdx rgn_id;
};

struct AncestorCollection {
    bit_vector exists_bv;
    rank_support_v<> exists_rb;
    vector<AncestorIdx> ancestor_indices;
    vector<uint8_t> ancestor_counts;
    vector<uint8_t> ancestor_local_indices;

    // this is merged stratum - ancestors collection. one - to - many.
    // ancestor index -> points to some position, and then consecutive ancestor_counts[idx] is all ancestors of the stratum.
    vector<uint8_t> ancestors__is_seq;
    vector<StratumOrSeqId> ancestors__id;
    vector<CoveringRegionIdx> ancestors__rgn_id;

    AncestorIdx register_ancestor_count(StratumId stratum_id, uint8_t right_count){
        this->exists_bv[stratum_id]=true;
        AncestorIdx idx = this->ancestors__is_seq.size();
        ancestor_indices.push_back(idx);
        ancestor_counts.push_back(right_count);
        ancestor_local_indices.push_back(0);

        this->ancestors__is_seq.resize(idx+right_count);
        this->ancestors__id.resize(idx+right_count);
        this->ancestors__rgn_id.resize(idx+right_count);
        return idx;
    }

    void prepare_fetch(){
        this->exists_rb=rank_support_v<>(&this->exists_bv);
    }
    
    bool has_ancestor(StratumId stratum_id){
        return this->exists_bv[stratum_id];
    }

    AncestorIdx get_ancestor_index(StratumId stratum_id){
        assert(this->exists_bv[stratum_id]);
        return this->ancestor_indices[this->exists_rb.rank(stratum_id)];
    }

    void set_ancestor(StratumId stratum_id, uint8_t is_seq, StratumOrSeqId id, CoveringRegionIdx rgn_idx){
        assert(this->exists_bv[stratum_id]);
        auto ancestor_idx = this->ancestor_indices[this->exists_rb.rank(stratum_id)];
        assert(this->ancestor_local_indices[ancestor_idx]+1 <= this->ancestor_counts[ancestor_idx]); 
        auto ancestor_element_idx = ancestor_idx + this->ancestor_local_indices[ancestor_idx];
        this->ancestors__is_seq[ancestor_element_idx]=is_seq;
        this->ancestors__id[ancestor_element_idx]=id;
        this->ancestors__rgn_id[ancestor_element_idx]=rgn_idx;
        this->ancestor_local_indices[ancestor_idx]++;
    }

    void get_ancestors(StratumId stratum_id, vector<Ancestor> &ancestors){
        assert(this->exists_bv[stratum_id]);
        auto ancestor_idx = this->ancestor_indices[this->exists_rb.rank(stratum_id)];
        auto ancestor_count = this->ancestor_counts[ancestor_idx];
        ancestors.resize(ancestor_count);
        for(int i=0; i< ancestor_count; i++){
            ancestors[i].is_seq = this->ancestors__is_seq[ancestor_idx+i];
            ancestors[i].id = this->ancestors__id[ancestor_idx+i];
            ancestors[i].rgn_id = this->ancestors__rgn_id[ancestor_idx+i];
        }
    }
};

struct DescendantCollection {
    // vector<StratumId> sequences__left_descendant;
    vector<StratumId> strata__left_descendant;
    vector<StratumId> strata__right_descendant;
};

struct CompactedDBGWorkspace {
    AncestorCollection ancestor;
    DescendantCollection descendant;

    vector<uint8_t> is_strata_right_descendent_set;
    vector<uint8_t> is_strata_left_descendent_set;

    SpinLock ancestor_lock;

    CompactedDBGWorkspace(Prokrustean &prokrustean){
        this->reset(prokrustean);
    }

    void reset(Prokrustean &prokrustean){
        this->ancestor.exists_bv.resize(prokrustean.stratum_count);
        // this->descendant.sequences__left_descendant.resize(prokrustean.sequence_count);
        this->descendant.strata__left_descendant.resize(prokrustean.stratum_count);
        this->descendant.strata__right_descendant.resize(prokrustean.stratum_count);
        this->is_strata_right_descendent_set.resize(prokrustean.stratum_count);
        this->is_strata_left_descendent_set.resize(prokrustean.stratum_count);
    }

    // may need global lock
    AncestorIdx register_ancestor_count(StratumId stratum_id, uint8_t right_count){
        AncestorIdx idx = this->ancestor.register_ancestor_count(stratum_id, right_count);
        return idx;
    }

    // may need stratum lock
    void register_ancestor(StratumId stratum_id, uint8_t is_seq, StratumOrSeqId id, CoveringRegionIdx rgn_idx){
        this->ancestor.set_ancestor(stratum_id, is_seq, id, rgn_idx);
    }

    // void register_left_descendant(bool is_seq, StratumOrSeqId id, StratumId descendant_id){
    //     if(is_seq){
    //         this->descendant.sequences__left_descendant[id]=descendant_id;
    //     } else {
    //         this->descendant.strata__left_descendant[id]=descendant_id;
    //     }
    // }
    void register_left_descendant(StratumId id, StratumId descendant_id){
        this->descendant.strata__left_descendant[id]=descendant_id;
    }

    void register_right_descendant(StratumId id, StratumId descendant_id){
        this->descendant.strata__right_descendant[id]=descendant_id;
    }

    StratumId get_left_descendant(StratumId id){
        return this->descendant.strata__left_descendant[id];
    }

    StratumId get_right_descendant(StratumId id){
        return this->descendant.strata__right_descendant[id];
    }

    void verify_all_ancestors_set(){
        for(int i=0; i<this->ancestor.ancestor_indices.size(); i++){
            assert(this->ancestor.ancestor_local_indices[i]==this->ancestor.ancestor_counts[i]);
        }
    }
    void verify_all_descendants_set(){
        for(auto val: this->is_strata_left_descendent_set){
            assert(val==1);
        }
        for(auto val: this->is_strata_right_descendent_set){
            assert(val==1);
        }
    }
};

void _set_deepest_descendents_stra(StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work, vector<StratumId> &working_stratum_ids){
    StratifiedEdge edge;
    auto &prokrustean=ext.prokrustean;
    
    working_stratum_ids.clear();
    StratumId curr_stratum_id=stratum_id;

    while(true){
        // already visited
        if(work.is_strata_right_descendent_set[curr_stratum_id]==1 && working_stratum_ids.size()>0){
            for(auto id: working_stratum_ids){
                work.register_right_descendant(id, work.get_right_descendant(curr_stratum_id));
            }
            break;
        }
        // set visited.
        work.is_strata_right_descendent_set[curr_stratum_id]=1;
        if(!ext.is_last_reflecting(curr_stratum_id, k-1)){
            working_stratum_ids.push_back(curr_stratum_id);
            // must exists
            assert(ext.get_stratum_last_stratified(curr_stratum_id, edge));
            curr_stratum_id=edge.stratum_id;
            continue;
        } else {
            for(auto id: working_stratum_ids){
                work.register_right_descendant(id, curr_stratum_id);
            }
            work.register_right_descendant(curr_stratum_id, curr_stratum_id);
            break;
        }
    }
    // set leftmost descendent
    working_stratum_ids.clear();
    StratumId curr_stratum_id=stratum_id;

    while(true){
        // already visited
        if(work.is_strata_left_descendent_set[curr_stratum_id]==1 && working_stratum_ids.size()>0){
            for(auto id: working_stratum_ids){
                work.register_left_descendant(id, curr_stratum_id);
            }
            break;
        }
        // set visited.
        work.is_strata_left_descendent_set[curr_stratum_id]=1;
        if(!ext.is_first_reflecting(curr_stratum_id, k-1)){
            working_stratum_ids.push_back(curr_stratum_id);
            // must exists
            assert(ext.get_stratum_first_stratified(curr_stratum_id, edge));
            curr_stratum_id=edge.stratum_id;
            continue;
        } else {
            for(auto id: working_stratum_ids){
                work.register_left_descendant(id, curr_stratum_id);
            }
            work.register_left_descendant(curr_stratum_id, curr_stratum_id);
            break;
        }
    }
}

void _set_ancestor_counts(StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    if(ext.prokrustean.get_right_cnt(stratum_id)>1 && ext.is_last_reflecting(stratum_id, k-1)){
        work.register_ancestor_count(stratum_id, ext.prokrustean.get_right_cnt(stratum_id));
    }
}

// StratumId stratum_id, uint8_t is_seq, StratumOrSeqId id, CoveringRegionIdx rgn_idx
void _set_ancestors(Vertex &vertex, ProkrusteanExtension &ext, CompactedDBGWorkspace &work, stack<StratumId> &stratum_stack, StratifiedEdge &working_edge){
    for(int e=0; e<vertex.s_edges.size(); e++){
        auto &edge=vertex.s_edges[e];
        if(edge.to<vertex.size){
            auto limited_length=0;
            if(e+1<vertex.s_edges.size() && vertex.s_edges[e+1].from<edge.to){
                // intersection exists
                limited_length=edge.to-vertex.s_edges[e+1].from;
            }
            stratum_stack.push(edge.stratum_id);
            while(!stratum_stack.empty()){
                auto stratum_id=stratum_stack.top();
                stratum_stack.pop();
                
                if(work.ancestor.has_ancestor(stratum_id)){
                    work.register_ancestor(stratum_id, vertex.is_sequence, vertex.id, e);
                }
                if(ext.get_stratum_last_stratified(stratum_id, working_edge)){
                    if(working_edge.to==ext.prokrustean.get_stratum_size(stratum_id)
                    && limited_length < working_edge.size()){
                        stratum_stack.push(working_edge.stratum_id);
                    }
                }
            }
        }       
    }
}




// void _build_string_of_unitig(UnitigId uni_id, vector<Unitig> &unitigs, AbstractSequenceAccess &seq_access, int k){
//     assert(unitigs[uni_id].is_start_of_maximal);
//     // implement string
//     unitigs[uni_id].content=unitigs[uni_id].get_string(seq_access);
//     if(unitigs[uni_id].next_cnt!=1){
//         return;
//     }
//     UnitigId curr_id=uni_id;
//     // UnitigId next_id=unitigs[curr_id].nexts[0];
//     UnitigId next_id=unitigs[curr_id].nexts_new[0];
//     while(true){
//         Unitig &next_unitig=unitigs[next_id];
//         // normal case
//         if(next_unitig.is_void_k_minus_1_unitig==false){
//             if(next_unitig.is_convergence){
//                 // unitigs[uni_id].nexts.clear();
//                 // unitigs[uni_id].nexts.push_back(next_id);
//                 unitigs[uni_id].set_next_capacity(1);
//                 unitigs[uni_id].set_next(next_id);
//                 break;
//             } else {
//                 // if(unitigs[curr_id].nexts.size()<=1){
//                 if(unitigs[curr_id].next_cnt<=1){
//                     // merge
//                     unitigs[uni_id].content+=next_unitig.get_string(seq_access).substr(k-1);
//                 }

//                 // if(next_unitig.nexts.size()!=1){
//                 if(next_unitig.next_cnt!=1){
//                     // unitigs[uni_id].nexts=next_unitig.nexts;
//                     unitigs[uni_id].replace_nexts(next_unitig);
//                     break;
//                 } else {
//                     curr_id=next_id;
//                     // next_id=next_unitig.nexts[0];
//                     next_id=next_unitig.nexts_new[0];
//                 }
//             }
//         } else {
//             // void case
//             // if(next_unitig.nexts.size()!=1){
//             if(next_unitig.next_cnt!=1){
//                 // unitigs[uni_id].nexts=next_unitig.nexts;
//                 unitigs[uni_id].replace_nexts(next_unitig);
//                 break;
//             // } else if(next_unitig.nexts.size()==1){
//             } else if(next_unitig.next_cnt==1){
//                 curr_id=next_id;
//                 // next_id=next_unitig.nexts[0];
//                 next_id=next_unitig.nexts_new[0];
//             }
//         }
//     }
// }


void cdbg_stage1_set_right_extendable_structure(int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work, bool verbose=false){
    assert(k>1);
    assert(ext.prokrustean.contains_stratum_extension_count);
    int stratum_count=ext.prokrustean.stratum_count;
    int seq_count=ext.prokrustean.sequence_count;
    
    work.reset(ext.prokrustean);
    
    Vertex working_vertex;
    StratifiedEdge working_edge;
    vector<StratumId> working_stratum_ids;
    stack<StratumId> working_stratum_stack;
    

    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            _set_ancestor_counts(i, k, ext, work);
        }
    }

    work.ancestor.prepare_fetch();

    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            ext.prokrustean.get_stratum(i, working_vertex);
            _set_ancestors(working_vertex, ext, work, working_stratum_stack, working_edge);
        }
    }

    for(int i=0; i<seq_count; i++){
        if(ext.prokrustean.sequences__size[i]>=k-1){
            ext.prokrustean.get_sequence(i, working_vertex);
            _set_ancestors(working_vertex, ext, work, working_stratum_stack, working_edge);
        }
    }
    
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            _set_deepest_descendents_stra(i, k, ext, work, working_stratum_ids);
        }
    }
}

// void cdbg_stage3_construct_cdbg(vector<Unitig> &unitigs, AbstractSequenceAccess &seq_access, AbstractStringDataStore &store, int k){
//     for(UnitigId id=0; id<unitigs.size(); id++){
//         auto &unitig = unitigs[id];
//         if(unitig.is_start_of_maximal){
//             _build_string_of_unitig(id, unitigs, seq_access, k);
//             // cout << " unitig.content " << unitig.content << endl;
//             // store.store(unitig.content);
//         }
//     }
//     store.store("------------------------------------------------------------------------------------------------");
//     store.store("---------------         columns: unitig id, next unitig ids, unitig content     ----------------");
//     store.store("------------------------------------------------------------------------------------------------");
//     for(UnitigId id=0; id<unitigs.size(); id++){
//         if(unitigs[id].is_start_of_maximal){
//             string expr=to_string((int)id);
//             for(int i=0; i<unitigs[id].next_cnt; i++){
//                 expr+="  ";
//                 expr+=to_string((int)unitigs[id].nexts_new[i]);
//             }
//             expr+="    ";
//             expr+=unitigs[id].content;
//             store.store(expr);
//         }
//     }
// }


void verify_right_extendable(int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work, uint64_t calculated_maximal_unitig_cnt){
    work.verify_all_ancestors_set();
    work.verify_all_descendants_set();

    uint64_t maximal_unitig_cnt=0;
    
    Vertex working_vertex;
    vector<Ancestor> working_ancestors;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        if(ext.prokrustean.sequences__size[i]<k-1){
            continue;
        }
        if(ext.is_first_reflecting_seq(i, k-1)){
            maximal_unitig_cnt++;
        }
    }

    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        auto stratum_size=ext.prokrustean.stratums__size[i];
        if(stratum_size<k-1){
            continue;
        } else if(stratum_size>k-1){
            if(work.descendant.strata__left_descendant[i]==i && ext.prokrustean.get_left_cnt(i)>1){
                maximal_unitig_cnt++;
            }
            if(work.descendant.strata__left_descendant[i]==i && ext.prokrustean.get_left_cnt(i)==0){
                maximal_unitig_cnt++;
            }
            if(work.ancestor.has_ancestor(i)){
                work.ancestor.get_ancestors(i, working_ancestors);
                if(working_ancestors.size()>1){
                    maximal_unitig_cnt+= working_ancestors.size();
                }
            }
        } else {

        }
    }

    // int stratum_size=ext.prokrustean.stratums__size[stratum_id];
    // int cnt=0;
    // if(stratum_size<k-1){
    //     return 0;
    // } else if(stratum_size>k-1){
    //     if(ext.stratum__refracted_at_front_of_cover(stratum_id, k-1)){
    //         // convergence or tip
    //         cnt += ext.prokrustean.get_left_cnt(stratum_id)!=1? 1 : 0;
    //     }
    //     if(ext.stratum__refracted_at_back_of_cover(stratum_id, k-1)){
    //         // divergence
    //         cnt += ext.prokrustean.get_right_cnt(stratum_id)>1? ext.prokrustean.get_right_cnt(stratum_id) : 0;
    //     }
    // } else { // k-1 case
    //     if(ext.prokrustean.get_right_cnt(stratum_id)>1){
    //         // divergence implication
    //         cnt += ext.prokrustean.get_right_cnt(stratum_id);    
    //     } else if(ext.prokrustean.get_right_cnt(stratum_id)==1 &&  ext.prokrustean.get_left_cnt(stratum_id)>1){
    //         // convergence implication
    //         cnt += 1;
    //     }
    // }
    

    assert(maximal_unitig_cnt==calculated_maximal_unitig_cnt);
}



#endif