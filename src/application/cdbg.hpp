#ifndef APPLICATION_CDBG_HPP_
#define APPLICATION_CDBG_HPP_
#include <algorithm>
#include "../prokrustean.enhance.hpp"

/* 
* (0) Data Structure
* - Prokrustean (left/right counts)
* - Sequences
* - Stratum Additionals: First/Last unitig, Leftmost/Rightmost.
* - Sample occurrences
* - Unitig (marked maximal starting)
* (1) Iterate sequence/stratums (setup)
* - Set leftmost/rightmost
* - For each stratum.
*   - Make first and last unitig.
*   - The first unitig is marked maximal if converged(left_cnt>1) or tip(left_cnt>0).
*   - The last unitig is marked maximal if diverged(right_cnt)
*   - first/last unitig can be null.
* (2) Iterate sequence/stratums (explore)
* - For each reflectum.
*   - Register unitig if non-(first or last). (use registered one if first/last)
*   - Go both directions and make strings (cannot be over the starting sequence or stratum) 
*   - For both directions if converge/diverge or bottom then stop. 
*   - If diverge on left, mark maximal. itself next (lock) and has right context if not last.
*   - If not diverged on left, itself safely next and has the right context next if not last.
* consequence: 
* - Every reflectum is explored, connected by next.  
* - Drops
*   - Drop leftmost/rightmost
*   - Drop prokrustean.
*   - Drop stratum additional. 
* (3) Sample?
* - sample occurrences
* - switch node locations to sequence based. 
* - drop samples
* (4) Build
* - load sequences.
* - Iterate maximal unitigs and generate unitigs.
* - while next is single, merge. 
*/

struct Unitig {
    // careful.. but it seems in theory only downward means no jump in top level. can be extended  
    StratumOrSeqId loc_id; 
    Pos loc_from;
    Pos loc_to;
    bool is_start_of_maximal;
    bool is_from_stratum;
    bool is_void_k_minus_1_unitig;
    vector<Unitig*> nexts; // optimize
    
    void extend_left(Pos size){
        assert(loc_from>=size);
        loc_from-=size;
    };
    void extend_right(Pos size){
        assert(size>0);
        loc_to+=size;
    };

    void print(){
        cout << "loc id: " << loc_id << ", from: " << loc_from << ", to: " << loc_to << endl;
    }
};

struct StratumExt {
    StratumId leftmost_descendent;
    StratumId rightmost_descendent;
    optional<Unitig*> first_unitig;
    optional<Unitig*> last_unitig;
};

struct Workspace {
    Vertex working_vertex;
    vector<Region> working_bands;
    Unitig working_unitig;
    vector<StratumId> working_stratum_ids;

    vector<Unitig> unitigs;
    // vector<bool> is_unitig_from_stratum; 
    // vector<bool> is_unitig_maximal; 

    vector<StratumExt> stratum_exts;
    vector<bool> is_leftmost_descendent_set;
    vector<bool> is_rightmost_descendent_set;

    SpinLock unitig_gen_lock;
    UnitigId new_unitig_id=0;

    Workspace(uint64_t stratum_count){
        this->stratum_exts=vector<StratumExt>(stratum_count);
        this->is_leftmost_descendent_set=vector<bool>(stratum_count);
        this->is_rightmost_descendent_set=vector<bool>(stratum_count);
    }

    UnitigId make_unitig(){
        this->unitig_gen_lock.lock();
        /* critical region */
        UnitigId new_id=this->new_unitig_id;
        this->new_unitig_id++;
        this->unitigs.push_back(Unitig());
        // this->is_unitig_from_stratum.push_back(false);
        /* critical region */
        this->unitig_gen_lock.unlock();
        return new_id;
    }
};

void set_first_last_unitigs(StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    // optimize
    ext.prokrustean.get_stratum(stratum_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);
    // important. unitig for first made only when exploration stop here
    // if not reflected, no.
    // if not left count 1, no
    // but if bottom that bands is 1, yes 
    bool is_first_tip_or_converge=ext.stratum_left_ext_count[stratum_id]!=1;
    bool is_first_reflected=work.working_bands[0].is_reflected;
    bool make_first_unitig=(is_first_reflected && is_first_tip_or_converge) || work.working_bands.size()==1;
    if(make_first_unitig){
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[0].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[0].to;
        work.unitigs[unitig_id].is_from_stratum=true;
        //maximal-tips or convergence
        work.unitigs[unitig_id].is_start_of_maximal=is_first_tip_or_converge;
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=ext.prokrustean.stratums__size[stratum_id]==k-1;

        work.stratum_exts[stratum_id].first_unitig=&work.unitigs[unitig_id];
        // exceptional - first is last
        if(work.working_bands.size()==1){
            work.stratum_exts[stratum_id].last_unitig=&work.unitigs[unitig_id];
        }
        // work.unitigs[unitig_id].print();
    }
    // last if maximal
    auto last_band_idx=work.working_bands.size()-1;
    bool is_last_diverge=ext.stratum_right_ext_count[stratum_id]>1;
    bool make_last_unitig=last_band_idx>0 && work.working_bands[last_band_idx].is_reflected && is_last_diverge;
    if(make_last_unitig){
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[last_band_idx].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[last_band_idx].to;
        work.unitigs[unitig_id].is_from_stratum=true;

        // work.unitigs[unitig_id].print();
    }
}
void set_leftright_descendents(StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    auto &prokrustean=ext.prokrustean;
    if(!work.is_leftmost_descendent_set[stratum_id]){
        work.working_stratum_ids.clear();
        // work.working_stratum_ids.push_back(stratum_id);
        StratumId curr_stratum_id=stratum_id;
        while(true){
            // set visited.
            work.is_leftmost_descendent_set[curr_stratum_id]=true;
            // has large stratified at the front
            if(prokrustean.stratums__region_cnt[curr_stratum_id]>0
            && prokrustean.stratums__region[curr_stratum_id][0].pos==0){
                auto rgn_stra_id= prokrustean.stratums__region[curr_stratum_id][0].stratum_id;
                if(prokrustean.stratums__size[rgn_stra_id]>=k-1){
                    work.working_stratum_ids.push_back(curr_stratum_id);
                    curr_stratum_id=rgn_stra_id;
                    continue;
                } 
            }

            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].leftmost_descendent=curr_stratum_id;
            }
            work.stratum_exts[curr_stratum_id].leftmost_descendent=curr_stratum_id;
            break;
        }
    }

    if(!work.is_rightmost_descendent_set[stratum_id]){
        work.working_stratum_ids.clear();
        // work.working_stratum_ids.push_back(stratum_id);
        StratumId curr_stratum_id=stratum_id;
        while(true){
            // set visited.
            work.is_rightmost_descendent_set[curr_stratum_id]=true;
            // has large stratified at the last
            if(prokrustean.stratums__region_cnt[curr_stratum_id]>0){
                auto rgn_idx=prokrustean.stratums__region_cnt[curr_stratum_id]-1;
                auto rgn_stra_id=prokrustean.stratums__region[curr_stratum_id][rgn_idx].stratum_id;
                if(prokrustean.stratums__size[rgn_stra_id]>=k-1){
                    work.working_stratum_ids.push_back(curr_stratum_id);
                    curr_stratum_id=rgn_stra_id;
                    continue;
                } 
            }

            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].rightmost_descendent=curr_stratum_id;
            }
            work.stratum_exts[curr_stratum_id].rightmost_descendent=curr_stratum_id;
            break;
        }
    }
}

void _dig_left_and_extend(Unitig &unitig, StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    auto stratum_id_leftmost_desc=work.stratum_exts[stratum_id].leftmost_descendent;
    if(ext.stratum_right_ext_count[stratum_id_leftmost_desc]>1){
        assert(work.stratum_exts[stratum_id_leftmost_desc].first_unitig.has_value());
        unitig.nexts.push_back(work.stratum_exts[stratum_id_leftmost_desc].first_unitig.value());
    } else if(ext.prokrustean.stratums__region_cnt[stratum_id]==0){
        // bottom
        cout << "204 stratum_id_leftmost_desc: " << stratum_id_leftmost_desc <<endl;
        assert(work.stratum_exts[stratum_id_leftmost_desc].first_unitig.has_value());
        unitig.nexts.push_back(work.stratum_exts[stratum_id_leftmost_desc].first_unitig.value());
    } else {
        //extend
        auto first_rgn_stratum_size=ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[stratum_id][0].stratum_id];
        unitig.extend_right(first_rgn_stratum_size-(k-1));
        // there must be a stratified region on left
        auto right_stratum_id_of_first=ext.prokrustean.stratums__region[stratum_id][1].stratum_id;
        _dig_left_and_extend(unitig, right_stratum_id_of_first, k, ext, work);
    }
}
void _dig_right_and_extend(Unitig &unitig, StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    // set maximal if diverged
    auto stratum_id_rightmost_desc=work.stratum_exts[stratum_id].rightmost_descendent;
    if(ext.stratum_left_ext_count[stratum_id_rightmost_desc]>1){
        unitig.is_start_of_maximal=true;
        assert(work.stratum_exts[stratum_id_rightmost_desc].last_unitig.has_value());
        work.stratum_exts[stratum_id_rightmost_desc].last_unitig.value()->nexts.push_back(&unitig);
    } else if(ext.prokrustean.stratums__region_cnt[stratum_id]==0){
        // bottom
        assert(work.stratum_exts[stratum_id_rightmost_desc].last_unitig.has_value());
        work.stratum_exts[stratum_id_rightmost_desc].last_unitig.value()->nexts.push_back(&unitig);
    } else {
        //extend
        auto last_idx=ext.prokrustean.stratums__region_cnt[stratum_id]-1;
        auto stratum_size=ext.prokrustean.stratums__size[stratum_id];
        auto last_rgn_pos=ext.prokrustean.stratums__region[stratum_id][last_idx].pos;
        unitig.extend_left(stratum_size-last_rgn_pos-(k-1));
        // there must be a stratified region on left
        auto left_stratum_id_of_last=ext.prokrustean.stratums__region[stratum_id][last_idx-1].stratum_id;
        _dig_right_and_extend(unitig, left_stratum_id_of_last, k, ext, work);
    }
}

void extend_seq_unitigs(SeqId seq_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    ext.prokrustean.get_sequence(seq_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);
    // middle cases
    for(int i=0; i<work.working_bands.size(); i++){
        auto &band=work.working_bands[i];
        if(band.is_stratified){
            continue;
        }
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=seq_id;
        work.unitigs[unitig_id].loc_from=band.from;
        work.unitigs[unitig_id].loc_to=band.to;
        work.unitigs[unitig_id].is_from_stratum=false;    

        // tip
        if(i==0){
            work.unitigs[unitig_id].is_start_of_maximal=true;    

            if(work.working_bands.size()>1){
                _dig_left_and_extend(work.unitigs[unitig_id], work.working_bands[i+1].stratum_id, k, ext, work);
            }
        } else if(i==work.working_bands.size()-1){ 
            _dig_right_and_extend(work.unitigs[unitig_id], work.working_bands[i-1].stratum_id, k, ext, work);
        } else {
            // middle
            _dig_left_and_extend(work.unitigs[unitig_id], work.working_bands[i+1].stratum_id, k, ext, work);
            _dig_right_and_extend(work.unitigs[unitig_id], work.working_bands[i-1].stratum_id, k, ext, work);
        }
    }

}
void extend_stra_unitigs(StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    ext.prokrustean.get_stratum(stratum_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);

    if(work.stratum_exts[stratum_id].first_unitig.has_value()){
        auto unitig=work.stratum_exts[stratum_id].first_unitig.value();
        if(work.working_bands.size()>1){
            auto &right_stratieid=work.working_bands[1];
            assert(right_stratieid.is_stratified);
            _dig_left_and_extend(*unitig, right_stratieid.stratum_id, k, ext, work);
        }
    }
    if(work.stratum_exts[stratum_id].last_unitig.has_value()){
        auto unitig=work.stratum_exts[stratum_id].last_unitig.value();
        // can be duplicated with first unitig
        if(unitig->loc_from!=0){
            auto &left_stratified=work.working_bands[work.working_bands.size()-2];
            assert(left_stratified.is_stratified);
            _dig_right_and_extend(*unitig, left_stratified.stratum_id, k, ext, work);
        }
    }

    if(ext.prokrustean.stratums__region_cnt[stratum_id]<2){
        return;
    }

    // middle cases
    for(int i=1; i+1<work.working_bands.size(); i++){
        auto &band=work.working_bands[i];
        if(band.is_stratified){
            continue;
        }
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=stratum_id; // the parent stratum
        work.unitigs[unitig_id].loc_from=band.from;
        work.unitigs[unitig_id].loc_to=band.to;
        work.unitigs[unitig_id].is_from_stratum=true;
        
        auto &left_stratieid=work.working_bands[i-1];
        auto &right_stratieid=work.working_bands[i+1];
        assert(left_stratieid.is_stratified);
        assert(right_stratieid.is_stratified);
        _dig_left_and_extend(work.unitigs[unitig_id], left_stratieid.stratum_id, k, ext, work);
        _dig_right_and_extend(work.unitigs[unitig_id], right_stratieid.stratum_id, k, ext, work);
    }
}

void extract_compacted_dbg(int k, ProkrusteanEnhancement &ext, vector<string> &sequences, bool verbose=false){
    assert(k>1);
    int stratum_count=ext.prokrustean.stratum_count();
    int seq_count=ext.prokrustean.sequence_count();
    Workspace work(stratum_count);
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]<k-1){
            continue;
        }
        set_first_last_unitigs(i, k, ext, work);
        set_leftright_descendents(i, k, ext, work);
    }
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]<k-1){
            continue;
        }
        extend_stra_unitigs(i, k, ext, work);
    }
    // for(auto &unitig: work.unitigs){
    //     unitig.print();
    // }
    // for(int i=0; i< work.stratum_exts.size(); i++){
    //     auto &stra=work.stratum_exts[i];
    //     if(ext.prokrustean.stratums__size[i]>=k-1){
    //         cout << "stra("<< i <<"): " << "stra.leftmost_descendent: " << stra.leftmost_descendent << ", size: " << ext.prokrustean.stratums__size[stra.leftmost_descendent] << endl;
    //         cout << "stra("<< i <<"): " << "stra.rightmost_descendent: " << stra.rightmost_descendent << ", size: " << ext.prokrustean.stratums__size[stra.rightmost_descendent] << endl;
    //     }
    // }
    for(int i=0; i<seq_count; i++){
        if(ext.prokrustean.sequences__size[i]<k-1){
            continue;
        }
        extend_seq_unitigs(i, k, ext, work);
    }
    setup_stratum_example_occ(ext);
    for(auto &unitig: work.unitigs){
        if(unitig.is_from_stratum){
            auto seq_id=ext.stratum_sample_occ_seq_id[unitig.loc_id];
            auto pos=ext.stratum_sample_occ_pos[unitig.loc_id];
            cout << "unitig: " << sequences[seq_id].substr(pos, unitig.loc_to-unitig.loc_from) << endl;
        } else {
            cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << endl;
        }
    }
}

// * (2) Iterate sequence/stratums (explore)
// * - For each reflectum.
// *   - Register unitig if non-(first or last). (use registered one if first/last)
// *   - Go both directions and make strings (cannot be over the starting sequence or stratum) 
// *   - For both directions if converge/diverge or bottom then stop. 
// *   - If diverge on left, mark maximal. itself next (lock) and has right context if not last.
// *   - If not diverged on left, itself safely next and has the right context next if not last.

// void _connect_single_unitig(int k, Workspace &work, ProkrusteanEnhancement &ext, vector<string> &sequences){
//     /* Assume workspace is set.
//     *   Vertex working_vertex;
//     *   vector<Region> working_spectrum;
//     *   Unitig working_unitig;
//     */
//     /* Should be parallelizable without conflicts */
//     SeqId seq_id;
//     Pos pos;
    
//     if(work.working_vertex.is_stratum){
//         seq_id=ext.stratum_sample_occ_seq_id[work.working_vertex.id];
//         pos=ext.stratum_sample_occ_pos[work.working_vertex.id];
//     } else {
//         seq_id=work.working_vertex.id;
//         pos=0;
//     }
//     // spectrum[0] must be reflected because being the start of a unitig means it is a reflectum.
//     work.working_unitig.content=sequences[seq_id].substr(pos, work.working_bands[0].size());

//     int rgn_idx=0;
//     Region rgn;
//     while(true){
//         // find next in a sequence vertex
//         if(work.working_vertex.is_sequence){
//             // not last
//             if(rgn_idx+1<work.working_bands.size()){
//                 auto &next_rgn=work.working_bands[rgn_idx+1];
//                 if(next_rgn.is_stratified){
//                     auto next_stratum_id=work.pointers[next_rgn.stratum_id].leftmost_descendent;
//                     if(ext.stratum_left_ext_count[next_stratum_id]>1) {
//                         // convergence
//                         auto next_unitig_id=work.unitig_id_map[make_tuple(next_stratum_id, rgn_idx+1)];
//                         work.working_unitig.nexts.push_back(next_unitig_id);
//                     } else if(ext.stratum_left_ext_count[next_stratum_id]==1){
//                         // update work-spectrum
//                         ext.prokrustean.get_stratum(next_stratum_id, work.working_vertex);
//                         ext.prokrustean.get_spectrum(work.working_vertex, k, work.working_bands);
                        
//                         seq_id=ext.stratum_sample_occ_seq_id[next_stratum_id];
//                         pos=ext.stratum_sample_occ_pos[next_stratum_id];
//                         // important: skip k-1 and add the rest
//                         work.working_unitig.content+=sequences[seq_id].substr(pos-(k-1), work.working_bands[0].size());
//                         rgn_idx=0;
//                     } else {
//                         // cannot be 0 then no connection possible
//                         assert(false);
//                     }
//                 }
                
//             } else { 
//                 // last region, next does not exist
//                 // sequence is unrelated
//             }
//         } else {
//             // not last
//             if(rgn_idx+1<work.working_bands.size()){
//                 // next exists, must be stratified
//                 assert(work.working_bands[rgn_idx+1].is_stratified);
//                 auto next_stratum_id=work.pointers[work.working_bands[rgn_idx+1].stratum_id].leftmost_descendent;
//                 if(ext.stratum_left_ext_count[next_stratum_id]>1) {
//                     // convergence, stop. stratum id=unitig when convergence
//                     work.working_unitig.nexts.push_back(next_stratum_id);
//                 } else if(ext.stratum_left_ext_count[next_stratum_id]==1){
//                     // update work-spectrum
//                     ext.prokrustean.get_stratum(next_stratum_id, work.working_vertex);
//                     ext.prokrustean.get_spectrum(work.working_vertex, k, work.working_bands);
                    
//                     seq_id=ext.stratum_sample_occ_seq_id[next_stratum_id];
//                     pos=ext.stratum_sample_occ_pos[next_stratum_id];
//                     // important: skip k-1 and add the rest
//                     work.working_unitig.content+=sequences[seq_id].substr(pos-(k-1), work.working_bands[0].size());
//                     rgn_idx=0;
//                 } else {
//                     // cannot be 0 then no connection possible
//                     assert(false);
//                 }
//             } else { 
//                 // last region of stratum, next does not exist
//                 if(ext.stratum_right_ext_count[work.working_vertex.id]>1){
//                     // divergence
//                     for(auto& [seq_stra_id, band_id]: work.pointers[work.working_vertex.id].rightmost_ancestors){
                        
//                     }
//                 } else if(ext.stratum_right_ext_count[work.working_vertex.id]==1){

//                 } else {
//                     // finished
//                 }
//             }
//         }
//     }
// }
// void connect_unitigs(int k, ProkrusteanEnhancement &ext, vector<string> &sequences, vector<StratumPointers> &pointers, vector<Unitig> &output){
    
//     for(int i=0;i<ext.prokrustean.stratum_count(); i++){
//         if(pointers[i].visited){
//             continue;
//         }
//         bool is_head_of_unitig=ext.stratum_left_ext_count[i]!=1 && pointers[i].leftmost_descendent==i;
//         if(!is_head_of_unitig){
//             continue;
//         }
//         ext.prokrustean.get_stratum(i, vertex);
//         ext.prokrustean.get_spectrum(vertex, k, spectrum);
//         // must be reflectum is head.
//         assert(spectrum[0].is_reflected);
//         // initialize unitig
//         unitig.id=i;
//         auto seq_id=ext.stratum_sample_occ_seq_id[i];
//         auto pos=ext.stratum_sample_occ_pos[i];
//         auto stratum_size=ext.prokrustean.stratums__size[i];
//         unitig.content=sequences[seq_id].substr(pos, spectrum[0].size());

//         connect_single_unitig()
//         if(start.is_stratum){
//         ext.prokrustean.get_stratum(start.id, work.working_vertex);
//         ext.prokrustean.get_spectrum(work.working_vertex, k, work.working_spectrum);

//         seq_id=ext.stratum_sample_occ_seq_id[work.working_vertex.id];
//         pos=ext.stratum_sample_occ_pos[work.working_vertex.id];
//     } else {
//         ext.prokrustean.get_sequence(start.id, work.working_vertex);
//         ext.prokrustean.get_spectrum(work.working_vertex, k, work.working_spectrum);

//         seq_id=work.working_vertex.id;
//         pos=0;
//     }
//         int rgn_idx=0;
//         StratumId next_i
//         while(true){
//             // 
//             if(spectrum[rgn_idx].is_reflected){
//                 if(rgn_idx+1<spectrum.size()){
//                     // stratified
//                     assert(spectrum[rgn_idx+1].is_stratified);
//                     pointers[spectrum[rgn_idx+1].stratum_id].leftmost_descendent
//                 }
//             }
//         }
        
//     }
// }


#endif