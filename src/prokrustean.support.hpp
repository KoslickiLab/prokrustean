#ifndef PROKRUSTEAN_SUPPORT_HPP_
#define PROKRUSTEAN_SUPPORT_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "prokrustean.hpp"

using namespace std;

struct ProkrusteanSupport {
    Prokrustean& prokrustean;
    ProkrusteanOptional& prokrustean_optional;

    vector<tuple<SeqId, Pos>> stratum_occ_samples;

    // CDBG
    // for counting maximal unitigs of single k
    vector<optional<bool>> stratum__left_bound_stratum_has_single_right_extension;
    vector<bool> stratum__is_start_of_maximal_unitig;
    vector<bool> sequence__is_start_of_maximal_unitig;

    vector<StratumId> stratum_leftmost_descendent_id;
    vector<StratumId> stratum_rightmost_ancestor_id;
    vector<uint8_t> stratum_rightmost_ancestor_rgn_idx;

    bool ran__fill_stratum_left_bound_single_right_extension;
    
    ProkrusteanSupport(Prokrustean& prokrustean, ProkrusteanOptional& prokrustean_optional):prokrustean(prokrustean), prokrustean_optional(prokrustean_optional) {}

    void fill_stratum_left_bound_single_right_extension(int k){
        // for(int i=0; i<prokrustean.sequence_count(); i++){
        //     Sequence seq = prokrustean.get_sequence(i);
        //     seq.print();
        // }
        // for(int i=0; i<prokrustean.stratum_count(); i++){
        //     Stratum stra = prokrustean.get_stratum(i);
        //     stra.print();
        // }
        this->stratum__left_bound_stratum_has_single_right_extension=vector<optional<bool>>(prokrustean.stratum_count(), 0);
        vector<uint8_t> indices;
        
        /* single thread version of filling extensions */
        for(int i=0; i<prokrustean.sequence_count(); i++){
            Sequence seq = prokrustean.get_sequence(i);
            seq.get_valid_indices(k, indices);
            // has to pair two
            for(int r=0; r+1<indices.size(); r++){
                fill_for_two_(seq.s_regions[indices[r]], seq.s_regions[indices[r+1]], k);
            }
        }
        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            stra.get_valid_indices(k, indices);
            // has to pair two
            for(int r=0; r+1<indices.size(); r++){
                fill_for_two_(stra.s_regions[indices[r]], stra.s_regions[indices[r+1]], k);
            }
        }

        this->stratum__is_start_of_maximal_unitig=vector<bool>(prokrustean.stratum_count(), false);
        this->sequence__is_start_of_maximal_unitig=vector<bool>(prokrustean.sequence_count(), false);
        // count
        int cnt=0;
        for(int i=0; i<prokrustean.sequence_count(); i++){
            Sequence seq = prokrustean.get_sequence(i);
            
            if(seq.s_regions.size()>0 && seq.s_regions[0].size()>=k){
                this->sequence__is_start_of_maximal_unitig[i]=false;
            } else {
                this->sequence__is_start_of_maximal_unitig[i]=true;
                cnt++;
            }
        }

        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            
            if(stra.s_regions.size()>0 && stra.s_regions[0].size()>=k){
                this->stratum__is_start_of_maximal_unitig[i]=false;
            } else if(prokrustean_optional.stratum_left_ext_count[i]!=1) {
                this->stratum__is_start_of_maximal_unitig[i]=false;
            } else if(stratum__left_bound_stratum_has_single_right_extension[i].value()){
                this->stratum__is_start_of_maximal_unitig[i]=false;
            } else {
                this->stratum__is_start_of_maximal_unitig[i]=true;
                cnt++;
            }
        }

        cout << "maximal unitig counted: " << cnt << endl;
        ran__fill_stratum_left_bound_single_right_extension=true;
    }
    
    void fill_for_two_(StratifiedRegion& l_rgn, StratifiedRegion& r_rgn, int k){
        // go deepest to rightmost of left and leftmost of right.
        // if k is fixed, there is no way but go into the deepest stratum such that the rightmost/leftmost is not (stratified and >=k)
        StratumId l_id= get_deepest_stratum(l_rgn, k, false);
        bool left_desendant_is_single_extending=prokrustean_optional.stratum_right_ext_count[l_id]==1;
        // StratumId r_id= get_deepest_stratum(r_rgn, k, true);
        // bool left_descendant_is_diverging=prokrustean_optional.stratum_right_ext_count[l_id]!=1;
        auto r_id=r_rgn.stratum_id;
        bool has_to_go_deeper_recursion;
        while(true){
            //already computed 
            if(stratum__left_bound_stratum_has_single_right_extension[r_id].has_value()){
                break;
            }

            has_to_go_deeper_recursion=false;
            auto child_stratum=prokrustean.get_stratum(r_id);
            //leftmost shares edge
            if(child_stratum.s_regions.size()>0){
                auto& leftmost_rgn = child_stratum.s_regions[0]; 
                if(leftmost_rgn.from==0 && leftmost_rgn.size()>=k){
                    // has 1st stratified region that is large
                    has_to_go_deeper_recursion=true;
                    stratum__left_bound_stratum_has_single_right_extension[r_id]=false;
                    r_id=leftmost_rgn.stratum_id;
                }
            }
                
            if(!has_to_go_deeper_recursion){
                stratum__left_bound_stratum_has_single_right_extension[r_id]=left_desendant_is_single_extending;
                break;
            }
        }

    }

    StratumId get_deepest_stratum(StratifiedRegion& rgn, int k, bool left_false_right_true){
        auto stratum_id=rgn.stratum_id;
        while(true){
            auto child_stratum=prokrustean.get_stratum(stratum_id);
            if(left_false_right_true){
                //rightmost shares edge
                if(child_stratum.s_regions.size()>0){
                    auto& rightmost_rgn = child_stratum.s_regions[child_stratum.s_regions.size()-1]; 
                    if(rightmost_rgn.to==child_stratum.size && rightmost_rgn.size()>=k){
                        stratum_id=rightmost_rgn.stratum_id;
                        continue;
                    }
                }
                // stop
                break;
            } else {
                //leftmost shares edge
                if(child_stratum.s_regions.size()>0){
                    auto& leftmost_rgn = child_stratum.s_regions[0]; 
                    if(leftmost_rgn.from==0 && leftmost_rgn.size()>=k){
                        stratum_id=leftmost_rgn.stratum_id;
                        continue;
                    }
                }
                // stop
                break;
            }
        }
        return stratum_id;
    }

    void setup_stratum_example_occ(){
        stratum_occ_samples.resize(prokrustean.stratum_count());
        vector<bool> visits(prokrustean.stratum_count());
        std::stack<Stratum> stratum_stack;
        for(int i=0; i<prokrustean.sequence_count(); i++){
            for(auto &rgn: prokrustean.get_sequence(i).s_regions){
                auto stratum = prokrustean.get_stratum(rgn.stratum_id);
                stratum.set_occ(i, rgn.from);
                stratum_occ_samples[rgn.stratum_id]=stratum.example_occ.value();
                visits[rgn.stratum_id]=true;

                stratum_stack.push(stratum);
                while(!stratum_stack.empty()){
                    auto stratum=stratum_stack.top();
                    stratum_stack.pop();
                    for(auto &c_rgn: stratum.s_regions){
                        if(visits[c_rgn.stratum_id]) continue;
                        SeqId seq_id = get<0>(stratum.example_occ.value());
                        Pos rel_pos = get<1>(stratum.example_occ.value())+c_rgn.from;
                        auto c_stratum = prokrustean.get_stratum(c_rgn.stratum_id);
                        c_stratum.set_occ(seq_id, rel_pos);
                        stratum_occ_samples[c_rgn.stratum_id]=c_stratum.example_occ.value();
                        visits[c_rgn.stratum_id]=true;
                        stratum_stack.push(c_stratum);
                    }
                }
            }
        }
    }
};

#endif