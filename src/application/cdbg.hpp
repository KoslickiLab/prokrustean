#ifndef APPLICATION_CDBG_HPP_
#define APPLICATION_CDBG_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"
#include "../util/string.access.hpp"
#include "../util/data.store.hpp"

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
    bool is_start_of_maximal=false;
    bool is_convergence=false;
    bool is_from_stratum=false;
    // bool is_from_void_intersection=false;
    bool is_void_k_minus_1_unitig=false;
    // vector<UnitigId> nexts; // optimize

    UnitigId* nexts_new;
    uint8_t next_cnt=0;
    uint8_t next_cnt_capacity=0;

    string content;

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
    void print(vector<string> &sequences){
        cout << "unitig: " << sequences[loc_id].substr(loc_from, loc_to-loc_from) << " seq id: " << loc_id << endl;
    }
    string get_string(AbstractSequenceAccess &seq_access){
        string seq;
        seq_access.read_seq(loc_id, seq);
        bool valid=0<=loc_from && loc_to<=seq.size();
        if(!valid){
            cout << "sequence location invalid: " << loc_id << ": " << loc_from << ": " << loc_to << "(seq size " << seq.size() << ")" << endl;
            assert(false);
        }
        auto str=seq.substr(loc_from, loc_to-loc_from);
        return str;
    }
    char get_front_letter(vector<string> &sequences, int k){
        return sequences[loc_id][loc_from+(k-2)];
    }
    char get_last_letter(vector<string> &sequences, int k){
        return sequences[loc_id][loc_to-(k-3)];
    }
    
    void set_next(UnitigId next_id){
        assert(next_cnt<next_cnt_capacity);
        if(this->nexts_new==nullptr){
            this->nexts_new = new UnitigId[this->next_cnt_capacity];    
        }
        this->nexts_new[this->next_cnt]=next_id;
        this->next_cnt++;
    }
    void set_next_capacity(uint8_t capacity){
        this->next_cnt=0;
        if(next_cnt_capacity==capacity){
            return;
        } else {
            if(this->nexts_new!=nullptr){
                delete this->nexts_new;
            }
        }
        this->nexts_new = new UnitigId[capacity];
        this->next_cnt_capacity=capacity;
    }
    void replace_nexts(Unitig &merged_unitig){
        this->next_cnt=merged_unitig.next_cnt;
        this->next_cnt_capacity=merged_unitig.next_cnt_capacity;
        if(this->nexts_new!=nullptr){
            delete this->nexts_new;
        }
        this->nexts_new=merged_unitig.nexts_new;
    }
};

struct StratumExt {
    StratumId descendent_of_leftmost_refracted;
    StratumId descendent_of_rightmost_refracted;
    optional<UnitigId> first_unitig;
    optional<UnitigId> last_unitig;
};

struct CompactedDBGWorkspace {
    Vertex working_vertex;
    vector<Edge> working_bands;
    
    vector<StratumId> working_stratum_ids;
    vector<Unitig> unitigs;
    // vector<bool> is_unitig_from_stratum; 
    // vector<bool> is_unitig_maximal; 

    vector<StratumExt> stratum_exts;
    vector<uint8_t> is_leftmost_descendent_set;
    vector<uint8_t> is_rightmost_descendent_set;

    SpinLock unitig_gen_lock;
    UnitigId new_unitig_id=0;

    void reset(int stratum_count){
        this->new_unitig_id=0;
        this->unitigs.clear();
        this->stratum_exts.clear();
        this->is_leftmost_descendent_set.clear();
        this->is_rightmost_descendent_set.clear();
        this->stratum_exts.resize(stratum_count);
        this->is_leftmost_descendent_set.resize(stratum_count, 0);
        this->is_rightmost_descendent_set.resize(stratum_count, 0);
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

void _set_first_last_unitigs(StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    // optimize
    ext.prokrustean.get_stratum(stratum_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);
    // important. unitig for first made only when exploration stop here
    // if not reflected, no.
    // if not left count 1, no
    // but if bottom that bands is 1, yes 
    bool is_first_tip=ext.prokrustean.get_left_cnt(stratum_id)==0;
    bool is_first_convergence=ext.prokrustean.get_left_cnt(stratum_id)>1;
    bool is_first_reflected=work.working_bands[0].is_reflected;
    bool make_first_unitig=(
        is_first_reflected && (is_first_tip || is_first_convergence) // maximal
        ) || work.working_bands.size()==1; // bottom
    if(make_first_unitig){
        auto unitig_id=work.make_unitig();
        // cout << "unitig make for: " << stratum_id << endl;
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[0].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[0].to;
        work.unitigs[unitig_id].is_from_stratum=true;
        work.unitigs[unitig_id].is_convergence=is_first_convergence;
        work.unitigs[unitig_id].set_next_capacity(1); 
        //maximal-tips or convergence
        if(ext.prokrustean.stratums__size[stratum_id]>k-1){
            // maximal case 3
            if(is_first_convergence) work.unitigs[unitig_id].is_start_of_maximal=true;
            // maximal case 2
            if(is_first_tip) work.unitigs[unitig_id].is_start_of_maximal=true;
        } 
        // the other place covers it (the case 6)
        // else {
        //     if(is_first_convergence && ext.prokrustean.get_right_cnt(stratum_id)==1){
        //         // cout << "convergence at stratum of k-1 " << i << endl;  
        //         cout << "hihi " << endl;
        //         work.unitigs[unitig_id].is_start_of_maximal=true;    
        //     }
        //     work.unitigs[unitig_id].is_void_k_minus_1_unitig=true;    
        // }

        work.stratum_exts[stratum_id].first_unitig=unitig_id;
        // exceptional - first is last
        if(work.working_bands.size()==1){
            work.stratum_exts[stratum_id].last_unitig=unitig_id;
            work.unitigs[unitig_id].is_void_k_minus_1_unitig=ext.prokrustean.stratums__size[stratum_id]==k-1;
            work.unitigs[unitig_id].set_next_capacity(ext.prokrustean.get_right_cnt(stratum_id));
        }
    }
    // last if maximal
    auto last_band_idx=work.working_bands.size()-1;
    bool is_last_diverge=ext.prokrustean.get_right_cnt(stratum_id)>1;
    bool is_last_tip=ext.prokrustean.get_right_cnt(stratum_id)==0;
    // not single element band, reflected, diverged. Bottom case is covered above(exceptional).
    bool make_last_unitig=last_band_idx>0 && work.working_bands[last_band_idx].is_reflected && (is_last_diverge || is_last_tip );
    if(make_last_unitig){
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=stratum_id;
        work.unitigs[unitig_id].loc_from=work.working_bands[last_band_idx].from;
        work.unitigs[unitig_id].loc_to=work.working_bands[last_band_idx].to;
        work.unitigs[unitig_id].is_from_stratum=true;

        work.unitigs[unitig_id].is_void_k_minus_1_unitig=ext.prokrustean.stratums__size[stratum_id]==k-1;
        work.stratum_exts[stratum_id].last_unitig=unitig_id;
        work.unitigs[unitig_id].set_next_capacity(ext.prokrustean.get_right_cnt(stratum_id));
    }
}

void _set_deepest_descendents(StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    StratifiedEdge edge;
    auto &prokrustean=ext.prokrustean;
    
    work.working_stratum_ids.clear();
    StratumId curr_stratum_id=stratum_id;
    while(true){
        // already visited
        if(work.is_rightmost_descendent_set[curr_stratum_id]==1 && work.working_stratum_ids.size()>0){
            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].descendent_of_rightmost_refracted=work.stratum_exts[curr_stratum_id].descendent_of_rightmost_refracted;
            }
            break;
        }
        // set visited.
        work.is_rightmost_descendent_set[curr_stratum_id]=1;
        // has large stratified at the last
        if(ext.get_stratum_last_stratified(curr_stratum_id, edge)
        && edge.to==ext.prokrustean.get_stratum_size(curr_stratum_id)
        && edge.size()>=k-1){
            work.working_stratum_ids.push_back(curr_stratum_id);
            curr_stratum_id=edge.stratum_id;
            continue;
        } else {
            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].descendent_of_rightmost_refracted=curr_stratum_id;
            }
            work.stratum_exts[curr_stratum_id].descendent_of_rightmost_refracted=curr_stratum_id;
            break;
        }
    }
    // set leftmost descendent
    work.working_stratum_ids.clear();
    curr_stratum_id=stratum_id;
    while(true){
        // already visited
        if(work.is_leftmost_descendent_set[curr_stratum_id]==1 && work.working_stratum_ids.size()>0){
            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].descendent_of_leftmost_refracted=work.stratum_exts[curr_stratum_id].descendent_of_leftmost_refracted;
            }
            break;
        }
        // set visited.
        work.is_leftmost_descendent_set[curr_stratum_id]=1;
        // has large stratified at the front
        if(ext.get_stratum_first_stratified(curr_stratum_id, edge)
        && edge.from==0
        && edge.size()>=k-1){
            work.working_stratum_ids.push_back(curr_stratum_id);
            curr_stratum_id=edge.stratum_id;
            continue;
        } else {
            for(auto id: work.working_stratum_ids){
                work.stratum_exts[id].descendent_of_leftmost_refracted=curr_stratum_id;
            }
            work.stratum_exts[curr_stratum_id].descendent_of_leftmost_refracted=curr_stratum_id;
        }
        break;
    }
}

void _dig_leftmosts_and_extend(UnitigId unitig_id, StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    Unitig &unitig=work.unitigs[unitig_id];
    // cout << " _dig_leftmosts_and_extend" << endl; 
    auto left_descendent_id=work.stratum_exts[stratum_id].descendent_of_leftmost_refracted;
    if(ext.prokrustean.stratums__size[left_descendent_id]>k-1){
        if(ext.prokrustean.get_left_cnt(left_descendent_id)>1){
            // assert(work.stratum_exts[left_descendent_id].first_unitig.has_value());
            unitig.extend_right(1);
            // unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
            unitig.set_next(work.stratum_exts[left_descendent_id].first_unitig.value());
            // assert(work.unitigs[work.unitigs[unitig_id].nexts[work.unitigs[unitig_id].nexts.size()-1]].loc_from==0);
        } else if(no_stratified_region_in_stra(left_descendent_id, k-1, ext)){
            // bottom
            assert(work.stratum_exts[left_descendent_id].first_unitig.has_value());
            unitig.extend_right(1);
            // unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
            unitig.set_next(work.stratum_exts[left_descendent_id].first_unitig.value());
            // assert(work.unitigs[work.unitigs[unitig_id].nexts[work.unitigs[unitig_id].nexts.size()-1]].loc_from==0);
        } else {
            //extend because left count = 1, not bottom. add reflectum.
            // read reflectum inline for efficiency. no first unitig is there by design to reduce space usage.
            auto idx = first_stratified_region_idx_in_stra(left_descendent_id, k-1, ext);
            auto &first_stratum_data=ext.prokrustean.stratums__region[left_descendent_id][idx];
            // auto first_stratum_pos=ext.prokrustean.stratums__region[left_descendent_id][idx].pos;
            // assert(first_stratum_data.pos>0);
            // push exactly 'gap'
            unitig.extend_right(first_stratum_data.pos);
            // try more to right
            auto right_stratum_id_of_first=first_stratum_data.stratum_id;
            _dig_leftmosts_and_extend(unitig_id, right_stratum_id_of_first, k, ext, work);
        }
    } else {
        // descendent is k-1 case --> assume it node and resolve later.
        unitig.extend_right(1);
        // unitig.nexts.push_back(work.stratum_exts[left_descendent_id].first_unitig.value());
        unitig.set_next(work.stratum_exts[left_descendent_id].first_unitig.value());
    }
}

void _dig_rightmosts_and_extend(UnitigId unitig_id, StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    if(ext.prokrustean.stratums__size[stratum_id]<k-1){
        return;
    }
    Unitig &unitig=work.unitigs[unitig_id];
    // cout << " _dig_rightmosts_and_extend" << endl;
    auto right_descendent_id=work.stratum_exts[stratum_id].descendent_of_rightmost_refracted;
    // cout << "rightmost: " << right_descendent_id << endl;
    if(ext.prokrustean.stratums__size[right_descendent_id]>k-1){
        // set maximal if diverged
        if(ext.prokrustean.get_right_cnt(right_descendent_id)>1){
            // maximal case 4
            unitig.is_start_of_maximal=true;
            assert(work.stratum_exts[right_descendent_id].last_unitig.has_value());
            unitig.extend_left(1);
            // work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
            auto right_descendent_last_id=work.stratum_exts[right_descendent_id].last_unitig.value();
            work.unitigs[right_descendent_last_id].set_next(unitig_id);
            // assert(ext.prokrustean.get_right_cnt(right_descendent_id)>=work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.size());
            
        } else if(no_stratified_region_in_stra(right_descendent_id, k-1, ext)){
            // bottom
            assert(work.stratum_exts[right_descendent_id].last_unitig.has_value());
            unitig.extend_left(1);
            // work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
            auto right_descendent_last_id=work.stratum_exts[right_descendent_id].last_unitig.value();
            work.unitigs[right_descendent_last_id].set_next(unitig_id);
            // work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].set_next(unitig_id);
        } else {
            //extend because right count = 1, not bottom. so add reflectum to left direction.
            // read reflectum inline for efficiency. no last unitig is there by design to reduce space usage.
            auto last_idx=last_stratified_region_idx_in_stra(right_descendent_id, k-1, ext);
            auto &last_stratified_data=ext.prokrustean.stratums__region[right_descendent_id][last_idx];
            auto last_stratified_to=last_stratified_data.pos+ext.prokrustean.stratums__size[last_stratified_data.stratum_id];
            // assert(last_stratum_pos>0);
            // push exactly 'gap'
            unitig.extend_left(ext.prokrustean.stratums__size[right_descendent_id]-last_stratified_to);
            // there must be a stratified region on left
            auto left_stratum_id_of_last=last_stratified_data.stratum_id;
            _dig_rightmosts_and_extend(unitig_id, left_stratum_id_of_last, k, ext, work);
        }
    } else {
        // k-1 case
        unitig.extend_left(1);
        // work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].nexts.push_back(unitig_id);
        
        auto right_descendent_last_id=work.stratum_exts[right_descendent_id].last_unitig.value();
        work.unitigs[right_descendent_last_id].set_next(unitig_id);
        // work.unitigs[work.stratum_exts[right_descendent_id].last_unitig.value()].set_next(unitig_id);
        
        if(ext.prokrustean.get_right_cnt(right_descendent_id)>1){
            // divergence multiple -> convergence does not matter
            // maximal case 5
            unitig.is_start_of_maximal=true;
        } else if(ext.prokrustean.get_right_cnt(right_descendent_id)==1){
            // convergence
            if(ext.prokrustean.get_left_cnt(right_descendent_id)>1){
                // maximal case 6
                unitig.is_start_of_maximal=true;
                unitig.is_convergence=true;
                // the others be literally void: 
                // * if right cnt multiple -> then each branch gets maximal. 
                // * if right cnt single/zero -> impossible
            } else if(ext.prokrustean.get_left_cnt(right_descendent_id)==0){
                assert(false); // tip stratum - theoretically seems impossible but let's check.
            }
        }
    }
}

void _extend_seq_unitigs(SeqId seq_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    ext.prokrustean.get_sequence(seq_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);
    // if(seq_id==0){
    //     for(auto &band: work.working_bands){
    //         band.print();
    //     }
    // }
    for(int i=0; i<work.working_bands.size(); i++){
        auto &band=work.working_bands[i];
        if(!band.is_reflected){
            continue;
        }
        assert(band.size()>=k-1);
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=seq_id;
        work.unitigs[unitig_id].loc_from=band.from;
        work.unitigs[unitig_id].loc_to=band.to;
        work.unitigs[unitig_id].is_from_stratum=false;    
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=false; // cannot be true
        work.unitigs[unitig_id].set_next_capacity(1);
        // tip
        if(i==0){
            // maximal case 1
            work.unitigs[unitig_id].is_start_of_maximal=true; 
            if(work.working_bands.size()>1){
                //since k-1 is applied, move pos right by 1
                _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            }
        } else if(i==work.working_bands.size()-1){
            //since k-1 is applied, move pos left by 1
            // work.unitigs[unitig_id].loc_from-=1; 
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i-1].stratum_id, k, ext, work);
            
        } else {
            //since k-1 is applied, move pos left/right by 1
            // work.unitigs[unitig_id].loc_to+=1;
            // work.unitigs[unitig_id].loc_from-=1; 
            // middle
            _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i-1].stratum_id, k, ext, work);
        }
    }

    // connect void intersections
    for(int i=0; i+1<work.working_bands.size(); i++){
        if(work.working_bands[i].is_stratified 
        && work.working_bands[i+1].is_stratified
        // && work.stratum_exts[work.working_bands[i].stratum_id].rightmost_reflectum_descendent==work.working_bands[i].stratum_id
        // && work.stratum_exts[work.working_bands[i+1].stratum_id].leftmost_reflectum_descendent==work.working_bands[i+1].stratum_id
        && work.working_bands[i].to - work.working_bands[i+1].from==k-2){
            // populates length k unitig.
            auto unitig_id=work.make_unitig();
            work.unitigs[unitig_id].loc_id=seq_id; 
            work.unitigs[unitig_id].loc_from=work.working_bands[i+1].from;
            work.unitigs[unitig_id].loc_to=work.working_bands[i].to;
            work.unitigs[unitig_id].is_from_stratum=false;
            work.unitigs[unitig_id].is_void_k_minus_1_unitig=false;
            work.unitigs[unitig_id].set_next_capacity(std::max<uint8_t>(1, ext.prokrustean.get_right_cnt(work.working_bands[i].stratum_id)));
            // work.unitigs[unitig_id].is_from_void_intersection=true;

            // work.unitigs[unitig_id].loc_from-=1;
            // work.unitigs[unitig_id].loc_to+=1;
            _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i].stratum_id, k, ext, work);
        }
    }
}

void _extend_stra_unitigs(StratumId stratum_id, int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    ext.prokrustean.get_stratum(stratum_id, work.working_vertex);
    ext.prokrustean.get_spectrum(work.working_vertex, k-1, work.working_bands);

    if(work.stratum_exts[stratum_id].first_unitig.has_value()){
        auto unitig_id=work.stratum_exts[stratum_id].first_unitig.value();
        if(work.working_bands.size()>1){
            auto &right_stratified=work.working_bands[1];
            assert(right_stratified.is_stratified);
            //since k-1 is applied, move pos right by 1
            // work.unitigs[unitig_id].loc_to+=1;
            _dig_leftmosts_and_extend(unitig_id, right_stratified.stratum_id, k, ext, work);
        }
    }
    if(work.stratum_exts[stratum_id].last_unitig.has_value()){
        auto unitig_id=work.stratum_exts[stratum_id].last_unitig.value();
        // can be duplicated with first unitig
        if(work.unitigs[unitig_id].loc_from!=0){
            auto &left_stratified=work.working_bands[work.working_bands.size()-2];
            assert(left_stratified.is_stratified);
            // cout << "stratum id: " << stratum_id << " left stratified id: " << left_stratified.stratum_id << endl;
            //since k-1 is applied, move pos left by 1
            // work.unitigs[unitig_id].loc_from-=1;
            _dig_rightmosts_and_extend(unitig_id, left_stratified.stratum_id, k, ext, work);
        }
    }

    if(ext.prokrustean.stratums__region_cnt[stratum_id]<2){
        return;
    }

    // middle cases
    for(int i=1; i+1<work.working_bands.size(); i++){
        auto &band=work.working_bands[i];
        if(!band.is_reflected){
            continue;
        }
        assert(band.size()>=k-1);
        auto unitig_id=work.make_unitig();
        work.unitigs[unitig_id].loc_id=stratum_id; // the parent stratum
        work.unitigs[unitig_id].loc_from=band.from;
        work.unitigs[unitig_id].loc_to=band.to;
        work.unitigs[unitig_id].is_from_stratum=true;
        work.unitigs[unitig_id].is_void_k_minus_1_unitig=false; // cannot be true
        work.unitigs[unitig_id].set_next_capacity(1);

        auto &left_stratieid=work.working_bands[i-1];
        auto &right_stratieid=work.working_bands[i+1];
        assert(left_stratieid.is_stratified);
        assert(right_stratieid.is_stratified);
        //since k-1 is applied, move pos left/right by 1
        // work.unitigs[unitig_id].loc_from-=1;
        // work.unitigs[unitig_id].loc_to+=1;
        _dig_rightmosts_and_extend(unitig_id, left_stratieid.stratum_id, k, ext, work);
        _dig_leftmosts_and_extend(unitig_id, right_stratieid.stratum_id, k, ext, work);
    }

    // connect void intersections
    for(int i=0; i+1<work.working_bands.size(); i++){
        if(work.working_bands[i].is_stratified 
        && work.working_bands[i+1].is_stratified
        // && work.stratum_exts[work.working_bands[i].stratum_id].rightmost_reflectum_descendent==work.working_bands[i].stratum_id
        // && work.stratum_exts[work.working_bands[i+1].stratum_id].leftmost_reflectum_descendent==work.working_bands[i+1].stratum_id
        && work.working_bands[i].to - work.working_bands[i+1].from==k-2){
            // populates length k unitig.
            auto unitig_id=work.make_unitig();
            work.unitigs[unitig_id].loc_id=stratum_id; // the parent stratum
            work.unitigs[unitig_id].loc_from=work.working_bands[i+1].from;
            work.unitigs[unitig_id].loc_to=work.working_bands[i].to;
            work.unitigs[unitig_id].is_from_stratum=true;
            work.unitigs[unitig_id].is_void_k_minus_1_unitig=false;
            // work.unitigs[unitig_id].is_from_void_intersection=true;
            work.unitigs[unitig_id].set_next_capacity(std::max<uint8_t>(1, ext.prokrustean.get_right_cnt(work.working_bands[i].stratum_id)));

            _dig_leftmosts_and_extend(unitig_id, work.working_bands[i+1].stratum_id, k, ext, work);
            _dig_rightmosts_and_extend(unitig_id, work.working_bands[i].stratum_id, k, ext, work);
        }
    }
}

void switch_loc(Unitig &unitig, ProkrusteanExtension &ext){
    //line order matter
    unitig.loc_from=ext.stratum_sample_occ_pos[unitig.loc_id]+unitig.loc_from;
    unitig.loc_to=ext.stratum_sample_occ_pos[unitig.loc_id]+unitig.loc_to;
    unitig.loc_id=ext.stratum_sample_occ_seq_id[unitig.loc_id]; 
}

void extract_paritial_unitigs(int k, ProkrusteanExtension &ext, CompactedDBGWorkspace &work, bool verbose=false){
    assert(k>1);
    assert(ext.prokrustean.contains_stratum_extension_count);
    int stratum_count=ext.prokrustean.stratum_count;
    int seq_count=ext.prokrustean.sequence_count;
    
    work.reset(stratum_count);
    
    cout << "first/last " << endl;
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            _set_first_last_unitigs(i, k, ext, work);
        }
    }
    cout << "stratum " << ext.prokrustean.stratums__size.size() << endl;
    int first_cnt=0;
    int last_cnt=0;
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        if(work.stratum_exts[i].first_unitig.has_value()){
            first_cnt++;
        }
        if(work.stratum_exts[i].last_unitig.has_value()){
            last_cnt++;
        }
    }
    cout << "first_cnt " <<first_cnt << " last_cnt  " << last_cnt << endl;
    cout << "deepest descendents " << endl;
    for(int i=0; i<stratum_count; i++){
        if(ext.prokrustean.stratums__size[i]>=k-1){
            _set_deepest_descendents(i, k, ext, work);            
        }
    }
    work.is_leftmost_descendent_set.clear();
    work.is_leftmost_descendent_set.shrink_to_fit();
    work.is_rightmost_descendent_set.clear();
    work.is_rightmost_descendent_set.shrink_to_fit();
    
    cout << "extend stra" << endl;
    for(int i=0; i<stratum_count; i++){
        // skip <k-2
        if(ext.prokrustean.stratums__size[i]>k-1){
            _extend_stra_unitigs(i, k, ext, work);
        }
    }
    
    cout << "extend seq" << endl;
    for(int i=0; i<seq_count; i++){
        if(ext.prokrustean.sequences__size[i]>k-1){
            _extend_seq_unitigs(i, k, ext, work);
        }
    }
}

void update_stratum_based_loc_to_seq_based_loc(ProkrusteanExtension &ext, CompactedDBGWorkspace &work){
    assert(ext.is_sample_occ_available);
    for(auto &unitig: work.unitigs){
        if(unitig.is_from_stratum){
            switch_loc(unitig, ext);
            // cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << " seq id: " << unitig.loc_id << endl;
        }
        else {
            // cout << "unitig: " << sequences[unitig.loc_id].substr(unitig.loc_from, unitig.loc_to-unitig.loc_from) << " seq id: " << unitig.loc_id << endl;
        }
    }
}

struct CdbgInvariants {
    int total_left_count=0;
    int total_pointers_pointing_down=0;

    int total_right_count=0;
    int total_pointers_pointing_up=0;

    int maximal_starting_unitig_count=0;
    int multiple_branches_include_a_void_branch=0;

    vector<optional<int>> each_left_extension_count_of_stratum; // reflectum at pos 0
    vector<optional<int>> each_unitig_attached_on_first_of_stratum_referenced_count; // reflectum at pos 0

    vector<optional<int>> each_unitig_attached_on_last_of_stratum_next_count; // reflectum at pos last
    vector<optional<int>> each_right_extension_count_of_stratum; // reflectum at pos last

    void set(CompactedDBGWorkspace &work, ProkrusteanExtension &ext){
        this->each_left_extension_count_of_stratum.resize(ext.prokrustean.stratum_count);
        this->each_right_extension_count_of_stratum.resize(ext.prokrustean.stratum_count);
        this->each_unitig_attached_on_first_of_stratum_referenced_count.resize(ext.prokrustean.stratum_count);
        this->each_unitig_attached_on_last_of_stratum_next_count.resize(ext.prokrustean.stratum_count);

        for(int i=0; i<work.stratum_exts.size(); i++){
            if(work.stratum_exts[i].last_unitig.has_value()){
                this->each_right_extension_count_of_stratum[i]=ext.prokrustean.get_right_cnt(i);
                this->total_right_count+=ext.prokrustean.get_right_cnt(i);

                this->each_unitig_attached_on_last_of_stratum_next_count[i]=work.unitigs[work.stratum_exts[i].last_unitig.value()].next_cnt;
                this->total_pointers_pointing_up+=work.unitigs[work.stratum_exts[i].last_unitig.value()].next_cnt;


            }
            if(work.stratum_exts[i].first_unitig.has_value()){
                this->each_left_extension_count_of_stratum[i]=ext.prokrustean.get_left_cnt(i);
                this->total_left_count+=ext.prokrustean.get_left_cnt(i);
            }
        }

        for(auto &unitig: work.unitigs){
            for(int i=0; i<unitig.next_cnt; i++){
                // auto next_id=unitig.nexts[i];
                auto next_id=unitig.nexts_new[i];
                // what is this next?
                if(work.unitigs[next_id].is_from_stratum){
                    auto &first_unitig=work.stratum_exts[work.unitigs[next_id].loc_id].first_unitig;
                    // check if next points to a first or not
                    if(first_unitig.has_value() && next_id==first_unitig.value()){
                        this->total_pointers_pointing_down++;
                        auto &referenced_cnt=this->each_unitig_attached_on_first_of_stratum_referenced_count[work.unitigs[next_id].loc_id];
                        referenced_cnt = referenced_cnt.has_value()? referenced_cnt.value()+1: 1;
                    }
                } else {
                    // let's not do for seq yet.
                }
            }
        }
        
        int void_cnt=0;
        // int void_intersection_cnt=0;
        for(auto &unitig: work.unitigs){
            // if(unitig.nexts.size()>1){
            if(unitig.next_cnt>1){
                // for(auto next: unitig.nexts){
                for(int i=0; i<unitig.next_cnt; i++){
                    if(work.unitigs[unitig.nexts_new[i]].is_void_k_minus_1_unitig){
                    // if(work.unitigs[unitig.nexts[i]].is_void_k_minus_1_unitig){
                        this->multiple_branches_include_a_void_branch++;
                    }
                }
            }
            if(unitig.is_void_k_minus_1_unitig){
                void_cnt++;
            }
            if(unitig.is_start_of_maximal){
                this->maximal_starting_unitig_count++;
            }
            // if(unitig.is_from_void_intersection){
            //     void_intersection_cnt++;
            // }
        }
    }

    void assert_result(){
        assert(total_right_count==total_pointers_pointing_up);
        assert(total_left_count==total_pointers_pointing_down);
        assert(multiple_branches_include_a_void_branch==0);
    }
    void print_result(){
        cout << "total left count: " << total_left_count << ", total pointing down: " << total_pointers_pointing_down << endl;
        cout << "total right count: " << total_right_count << ", total pointing down: " << total_pointers_pointing_up << endl;
    }
};


void build_strings(UnitigId uni_id, vector<Unitig> &unitigs, AbstractSequenceAccess &seq_access, int k){
    assert(unitigs[uni_id].is_start_of_maximal);
    // implement string
    unitigs[uni_id].content=unitigs[uni_id].get_string(seq_access);
    if(unitigs[uni_id].next_cnt!=1){
        return;
    }
    UnitigId curr_id=uni_id;
    // UnitigId next_id=unitigs[curr_id].nexts[0];
    UnitigId next_id=unitigs[curr_id].nexts_new[0];
    while(true){
        Unitig &next_unitig=unitigs[next_id];
        // normal case
        if(next_unitig.is_void_k_minus_1_unitig==false){
            if(next_unitig.is_convergence){
                // unitigs[uni_id].nexts.clear();
                // unitigs[uni_id].nexts.push_back(next_id);
                unitigs[uni_id].set_next_capacity(1);
                unitigs[uni_id].set_next(next_id);
                break;
            } else {
                // if(unitigs[curr_id].nexts.size()<=1){
                if(unitigs[curr_id].next_cnt<=1){
                    // merge
                    unitigs[uni_id].content+=next_unitig.get_string(seq_access).substr(k-1);
                }

                // if(next_unitig.nexts.size()!=1){
                if(next_unitig.next_cnt!=1){
                    // unitigs[uni_id].nexts=next_unitig.nexts;
                    unitigs[uni_id].replace_nexts(next_unitig);
                    break;
                } else {
                    curr_id=next_id;
                    // next_id=next_unitig.nexts[0];
                    next_id=next_unitig.nexts_new[0];
                }
            }
        } else {
            // void case
            // if(next_unitig.nexts.size()!=1){
            if(next_unitig.next_cnt!=1){
                // unitigs[uni_id].nexts=next_unitig.nexts;
                unitigs[uni_id].replace_nexts(next_unitig);
                break;
            // } else if(next_unitig.nexts.size()==1){
            } else if(next_unitig.next_cnt==1){
                curr_id=next_id;
                // next_id=next_unitig.nexts[0];
                next_id=next_unitig.nexts_new[0];
            }
        }
    }
}

void construct_cdbg(vector<Unitig> &unitigs, AbstractSequenceAccess &seq_access, AbstractStringDataStore &store, int k){
    for(UnitigId id=0; id<unitigs.size(); id++){
        auto &unitig = unitigs[id];
        if(unitig.is_start_of_maximal){
            build_strings(id, unitigs, seq_access, k);
            // cout << " unitig.content " << unitig.content << endl;
            store.store(unitig.content);
        }
    }
    // for(UnitigId id=0; id<unitigs.size(); id++){
    //     string expr=unitigs[id].content;
    //     for(int i=0; i<unitigs[id].next_cnt; i++){
    //         expr+=" ";
    //         expr+=unitigs[unitigs[id].nexts_new[i]].content;
    //     }
    //     store.store(expr);
    // }
}

#endif