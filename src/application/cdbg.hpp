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
    vector<UnitigId> nexts; // optimize
    Unitig(){}
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
    void print(vector<string> sequences){
        cout << "unitig: " << sequences[loc_id].substr(loc_from, loc_to-loc_from) << " seq id: " << loc_id << endl;
    }
};

struct StratumExt {
    StratumId leftmost_descendent;
    StratumId rightmost_descendent;
    optional<UnitigId> first_unitig;
    optional<UnitigId> last_unitig;
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
    bool is_first_tip_or_converge=ext.stratum_left_ext_count[stratum_id]!=1 && /* voids are exception */ ext.prokrustean.stratums__size[stratum_id]>k-1;
    bool is_first_reflected=work.working_bands[0].is_reflected;
    bool make_first_unitig=(is_first_reflected && is_first_tip_or_converge) || work.working_bands.size()==1;
    if(make_first_unitig){
        auto unitig_id=work.make_unitig();
        // cout << "unitig make for: " << stratum_id << endl;
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[0].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[0].to;
        work.unitigs[unitig_id].is_from_stratum=true;
        //maximal-tips or convergence
        work.unitigs[unitig_id].is_start_of_maximal=is_first_tip_or_converge;
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=ext.prokrustean.stratums__size[stratum_id]==k-1;

        work.stratum_exts[stratum_id].first_unitig=unitig_id;
        // exceptional - first is last
        if(work.working_bands.size()==1){
            work.stratum_exts[stratum_id].last_unitig=unitig_id;
            // cout << "last unitig set: "<< endl;
            // work.stratum_exts[stratum_id].last_unitig.value()->print();
        }
        // work.unitigs[unitig_id].print();
    }
    // last if maximal
    auto last_band_idx=work.working_bands.size()-1;
    bool is_last_diverge=ext.stratum_right_ext_count[stratum_id]>1;
    // not single element band, reflected, diverged. Bottom case is covered above(exceptional).
    bool make_last_unitig=last_band_idx>0 && work.working_bands[last_band_idx].is_reflected && is_last_diverge;
    if(make_last_unitig){
        auto unitig_id=work.make_unitig();
        // cout << "unitig make for: " << stratum_id << endl;
        // for(auto &band: work.working_bands){ band.print(); }
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[last_band_idx].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[last_band_idx].to;
        work.unitigs[unitig_id].is_from_stratum=true;

        work.unitigs[unitig_id].is_void_k_minus_1_unitig=ext.prokrustean.stratums__size[stratum_id]==k-1;
        work.stratum_exts[stratum_id].last_unitig=unitig_id;

        // cout << "last unitig set: "<< endl;
        // work.stratum_exts[stratum_id].last_unitig.value()->print();
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

void _dig_leftmosts_and_extend(UnitigId unitig_id, StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    Unitig &unitig=work.unitigs[unitig_id];
    // cout << " _dig_leftmosts_and_extend" << endl; 
    auto left_descendent_id=work.stratum_exts[stratum_id].leftmost_descendent;
    if(ext.prokrustean.stratums__size[stratum_id]>k-1){
        if(ext.stratum_left_ext_count[left_descendent_id]>1){
            assert(work.stratum_exts[left_descendent_id].first_unitig.has_value());
            unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
        } else if(ext.prokrustean.stratums__region_cnt[left_descendent_id]==0){
            // bottom
            assert(work.stratum_exts[left_descendent_id].first_unitig.has_value());
            unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
        } else {
            //extend because left count = 1, not bottom. add reflectum.
            // read reflectum inline for efficiency. no first unitig is there by design to reduce space usage.
            auto first_stratum_pos=ext.prokrustean.stratums__region[left_descendent_id][0].pos;
            assert(first_stratum_pos>0);
            // push exactly 'gap'
            unitig.extend_right(first_stratum_pos);
            // there must be a stratified region on left
            auto right_stratum_id_of_first=ext.prokrustean.stratums__region[stratum_id][1].stratum_id;
            _dig_leftmosts_and_extend(unitig_id, right_stratum_id_of_first, k, ext, work);
        }
    } else {
        cout << "k-1 maximal case " << (int)ext.stratum_right_ext_count[left_descendent_id] << endl;
        // descendent is k-1 case --> assume it node and resolve later.
        unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
    }
}

void _dig_rightmosts_and_extend(UnitigId unitig_id, StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    Unitig &unitig=work.unitigs[unitig_id];
    // cout << " _dig_rightmosts_and_extend" << endl;
    auto right_descendent_id=work.stratum_exts[stratum_id].rightmost_descendent;
    // cout << "rightmost: " << right_descendent_id << endl;
    if(ext.prokrustean.stratums__size[stratum_id]>k-1){
        auto rgn_cnt=ext.prokrustean.stratums__region_cnt[right_descendent_id];
        // set maximal if diverged
        if(ext.stratum_right_ext_count[right_descendent_id]>1){
            unitig.is_start_of_maximal=true;
            assert(work.stratum_exts[right_descendent_id].last_unitig.has_value());
            // cout << "nexts: " << work.stratum_exts[right_descendent_id].last_unitig.value()->nexts.size() << endl;
            work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
        } else if(rgn_cnt==0){
            // bottom
            assert(work.stratum_exts[right_descendent_id].last_unitig.has_value());
            work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
        } else {
            //extend because right count = 1, not bottom. so add reflectum to left direction.
            // read reflectum inline for efficiency. no last unitig is there by design to reduce space usage.
            auto last_stratified_to=ext.prokrustean.stratums__region[right_descendent_id][rgn_cnt-1].pos+ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[right_descendent_id][rgn_cnt-1].stratum_id];
            // assert(last_stratum_pos>0);
            // push exactly 'gap'
            unitig.extend_left(ext.prokrustean.stratums__size[right_descendent_id]-last_stratified_to);
            // there must be a stratified region on left
            auto left_stratum_id_of_last=ext.prokrustean.stratums__region[stratum_id][rgn_cnt-1-1].stratum_id;
            _dig_rightmosts_and_extend(unitig_id, left_stratum_id_of_last, k, ext, work);
        }
    } else {
        // k-1 case
        work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
        // if k-1 diverges or converges then maximal.
        if(ext.stratum_left_ext_count[stratum_id]>0 || ext.stratum_right_ext_count[stratum_id]>0){
            unitig.is_start_of_maximal=true;
        }
    }
}

void extend_seq_unitigs(SeqId seq_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    ext.prokrustean.get_sequence(seq_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);
    
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
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=false; // cannot be true

        // tip
        if(i==0){
            work.unitigs[unitig_id].is_start_of_maximal=true;    
            if(work.working_bands.size()>1){
                //since k-1 is applied, move pos right by 1
                work.unitigs[unitig_id].loc_to+=1;
                _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            }
        } else if(i==work.working_bands.size()-1){
            //since k-1 is applied, move pos left by 1
            work.unitigs[unitig_id].loc_from-=1; 
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i-1].stratum_id, k, ext, work);
        } else {
            //since k-1 is applied, move pos left/right by 1
            work.unitigs[unitig_id].loc_to+=1;
            work.unitigs[unitig_id].loc_from-=1; 
            // middle
            _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i-1].stratum_id, k, ext, work);
        }
    }

}

void extend_stra_unitigs(StratumId stratum_id, int k, ProkrusteanEnhancement &ext, Workspace &work){
    ext.prokrustean.get_stratum(stratum_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);

    if(work.stratum_exts[stratum_id].first_unitig.has_value()){
        auto unitig_id=work.stratum_exts[stratum_id].first_unitig.value();
        if(work.working_bands.size()>1){
            auto &right_stratified=work.working_bands[1];
            assert(right_stratified.is_stratified);
            //since k-1 is applied, move pos right by 1
            work.unitigs[unitig_id].loc_to+=1;
            _dig_leftmosts_and_extend(unitig_id, right_stratified.stratum_id, k, ext, work);
        }
    }
    if(work.stratum_exts[stratum_id].last_unitig.has_value()){
        auto unitig_id=work.stratum_exts[stratum_id].last_unitig.value();
        // can be duplicated with first unitig
        if(work.unitigs[unitig_id].loc_from!=0){
            auto &left_stratified=work.working_bands[work.working_bands.size()-2];
            assert(left_stratified.is_stratified);
            //since k-1 is applied, move pos left by 1
            work.unitigs[unitig_id].loc_from-=1;
            // if(unitig_id==2){
            //     cout << "GCA" << endl;
            //     work.working_vertex.print();
            //     for(auto &band: work.working_bands){
            //         band.print();
            //     }
            // }
            cout << "stratum id: " << stratum_id << " left stratified id: " << left_stratified.stratum_id << endl;
            _dig_rightmosts_and_extend(unitig_id, left_stratified.stratum_id, k, ext, work);
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
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=false; // cannot be true
        
        auto &left_stratieid=work.working_bands[i-1];
        auto &right_stratieid=work.working_bands[i+1];
        assert(left_stratieid.is_stratified);
        assert(right_stratieid.is_stratified);
        //since k-1 is applied, move pos left/right by 1
        work.unitigs[unitig_id].loc_from-=1;
        work.unitigs[unitig_id].loc_to+=1;
        _dig_leftmosts_and_extend(unitig_id, left_stratieid.stratum_id, k, ext, work);
        _dig_rightmosts_and_extend(unitig_id, right_stratieid.stratum_id, k, ext, work);
    }
}

void switch_loc(Unitig &unitig, ProkrusteanEnhancement &ext){
    //line order matter
    unitig.loc_from=ext.stratum_sample_occ_pos[unitig.loc_id]+unitig.loc_from;
    unitig.loc_to=ext.stratum_sample_occ_pos[unitig.loc_id]+unitig.loc_to;
    unitig.loc_id=ext.stratum_sample_occ_seq_id[unitig.loc_id]; 
}

void extract_compacted_dbg(int k, ProkrusteanEnhancement &ext, vector<string> &sequences, bool verbose=false){
    assert(k>1);
    int stratum_count=ext.prokrustean.stratum_count();
    int seq_count=ext.prokrustean.sequence_count();
    Workspace work(stratum_count);
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            set_first_last_unitigs(i, k, ext, work);
            set_leftright_descendents(i, k, ext, work);            
        }
    }
    for(int i=0; i<stratum_count; i++){
        // skip k-1
        if(ext.prokrustean.stratums__size[i]>k-1){
            extend_stra_unitigs(i, k, ext, work);
        }
    }
    setup_stratum_example_occ(ext);
    for(int i=0; i<ext.prokrustean.stratum_count(); i++){
        if(ext.prokrustean.stratums__size[i]!=k-1){
            continue;
        }
        auto seq_id=ext.stratum_sample_occ_seq_id[i];
        auto pos=ext.stratum_sample_occ_pos[i];
        cout << "stratum("<<i<<") " << sequences[seq_id].substr(pos, ext.prokrustean.stratums__size[i]) << " left cnt: " << (int)ext.stratum_left_ext_count[i] << " right cnt: " << (int)ext.stratum_right_ext_count[i] << endl;
    }
    // for(int i=0; i< work.stratum_exts.size(); i++){
    //     auto &stra=work.stratum_exts[i];
    //     if(ext.prokrustean.stratums__size[i]>=k-1){
    //         cout << "stra("<< i <<"): " << "stra.leftmost_descendent: " << stra.leftmost_descendent << ", size: " << ext.prokrustean.stratums__size[stra.leftmost_descendent] << endl;
    //         cout << "stra("<< i <<"): " << "stra.rightmost_descendent: " << stra.rightmost_descendent << ", size: " << ext.prokrustean.stratums__size[stra.rightmost_descendent] << endl;
    //     }
    // }
    for(int i=0; i<seq_count; i++){
        if(ext.prokrustean.stratums__size[i]>k-1){
            extend_seq_unitigs(i, k, ext, work);
        }
    }

    setup_stratum_example_occ(ext);
    for(auto &unitig: work.unitigs){
        if(unitig.is_from_stratum){
            switch_loc(unitig, ext);
            cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << " seq id: " << unitig.loc_id << endl;
        }
        else {
            cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << " seq id: " << unitig.loc_id << endl;
        }
    }

    int maximal_cnt=0;
    int void_cnt=0;
    for(auto &unitig: work.unitigs){
        if(unitig.is_void_k_minus_1_unitig){
            void_cnt++;
        }
        if(unitig.is_start_of_maximal){
            maximal_cnt++;
            cout << endl << "maximal start"<< endl;
            // unitig.print(sequences);
            // if(unitig.nexts.size()==1){
            //     UnitigId next_id=unitig.nexts[0];
            //     while(true){
            //         auto &next = work.unitigs[next_id];
            //         next.print(sequences);
            //         if(next.nexts.size()==1){
            //             next_id=next.nexts[0];
            //         } else {
            //             break;
            //         }
            //     }
            // }
        }
    }cout << "void: " << void_cnt<< endl;
    cout << "maximals: " << maximal_cnt<< endl;
    // for(int i=0; i<ext.prokrustean.stratum_count(); i++){
    //     auto seq_id=ext.stratum_sample_occ_seq_id[i];
    //     auto pos=ext.stratum_sample_occ_pos[i];
    //     cout << "stratum("<<i<<") " << sequences[seq_id].substr(pos, ext.prokrustean.stratums__size[i]) << endl;
    // }
    // for(auto &unitig: work.unitigs){
    //     if(unitig.is_from_stratum){
    //         auto seq_id=ext.stratum_sample_occ_seq_id[unitig.loc_id];
    //         auto pos=ext.stratum_sample_occ_pos[unitig.loc_id];
    //         cout << "unitig: " << sequences[seq_id].substr(pos+unitig.loc_from, unitig.loc_to-unitig.loc_from) << " stratum id: " << unitig.loc_id << endl;
    //     } else {
    //         cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << " seq id: " << unitig.loc_id << endl;
    //     }
    // }
}

#endif