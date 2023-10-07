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
    vector<optional<bool>> stratum__left_bound_stratum_is_branching;
    vector<bool> stratum__is_start_of_maximal_unitig;
    vector<bool> sequence__is_start_of_maximal_unitig;

    vector<StratumId> stratum_leftmost_descendent_id;
    vector<StratumId> stratum_rightmost_ancestor_id;
    vector<uint8_t> stratum_rightmost_ancestor_rgn_idx;

    bool ran__fill_stratum_left_bound_single_right_extension;
    
    ProkrusteanSupport(Prokrustean& prokrustean, ProkrusteanOptional& prokrustean_optional):prokrustean(prokrustean), prokrustean_optional(prokrustean_optional) {}

    int compute_maximal_unitigs(int k, vector<string> *sequences=nullptr){
        if(sequences!=nullptr){
            setup_stratum_example_occ();
            int idx=0;
            for(auto &pair: stratum_occ_samples){
                cout << "sample: " << (*sequences)[get<0>(pair)].substr(get<1>(pair), prokrustean.stratums__size[idx]) << endl;
                idx++;
            }
        }
        
        int cnt=0;
        vector<int> stats(7);
        vector<vector<string>> strings(7);
        // vector<uint8_t> indices;
        vector<Region> spectrum;
        /* single thread version of filling extensions */
        for(int i=0; i<prokrustean.sequence_count(); i++){
            Sequence seq = prokrustean.get_sequence(i);
            if(seq.size<k){
                continue;
            }
            // seq.print();
            // seq.get_valid_indices(k-1, indices);
            prokrustean.get_spectrum(seq, k-1, spectrum);
            if(spectrum[0].is_reflected){
                cnt++;
                stats[0]++;
                if(sequences!=nullptr){
                    auto txt = (*sequences)[i].substr(spectrum[0].from, spectrum[0].size());
                    strings[0].push_back(txt);
                }
            }
            // if(indices.size()>0 && seq.s_regions[indices[0]].from==0 && seq.s_regions[indices[0]].size()>=k-1){
            //     not tip
            // } else {
            //     cnt++;
            //     stats[0]++;
            // }
        }
        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            if(stra.size<k-1){
                continue;
            }
            // stra.print();
            // stra.get_valid_indices(k-1, indices);
            prokrustean.get_spectrum(stra, k-1, spectrum);
            if(stra.size>k-1){
                if(spectrum[0].is_reflected){
                    if(prokrustean_optional.stratum_left_ext_count[i]==0){
                        // tip
                        // cout << "tip at stratum " << i << endl; 
                        cnt++;
                        stats[1]++;
                        if(sequences!=nullptr){
                            auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(spectrum[0].from, spectrum[0].size());
                            strings[1].push_back(txt);
                        }
                    } else if(prokrustean_optional.stratum_left_ext_count[i]>1){
                        // convergence
                        // cout << "convergence at stratum " << i << endl; 
                        cnt++;
                        stats[2]++;
                        if(sequences!=nullptr){
                            auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(spectrum[0].from, spectrum[0].size());
                            strings[2].push_back(txt);
                        }
                    }
                } 
                // if(indices.size()>0 && stra.s_regions[indices[0]].from==0 && stra.s_regions[indices[0]].size()>=k-1){
                // } else {
                //     if(prokrustean_optional.stratum_left_ext_count[i]==0){
                //         // tip
                //         // cout << "tip at stratum " << i << endl; 
                //         cnt++;
                //         stats[1]++;
                //     } else if(prokrustean_optional.stratum_left_ext_count[i]>1){
                //         // convergence
                //         // cout << "convergence at stratum " << i << endl; 
                //         cnt++;
                //         stats[2]++;
                //     }
                // }
                if(spectrum[spectrum.size()-1].is_reflected){
                        if(prokrustean_optional.stratum_right_ext_count[i]>1){
                        // divergence
                        // cout << "divergence at stratum " << i << " (" << (int)prokrustean_optional.stratum_right_ext_count[i] << ")" << endl;  
                        cnt+=prokrustean_optional.stratum_right_ext_count[i];
                        stats[3]+=prokrustean_optional.stratum_right_ext_count[i];
                        if(sequences!=nullptr){
                            // auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(spectrum[0].from, spectrum[0].size());
                            auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(get<1>(stratum_occ_samples[i])+spectrum[spectrum.size()-1].from, spectrum[spectrum.size()-1].size());
                            txt+="("+to_string(prokrustean_optional.stratum_right_ext_count[i])+")";
                            strings[3].push_back(txt);
                            // for(int i=0; i< prokrustean_optional.stratum_right_ext_count[i]; i++)
                            // strings[3].push_back("div");
                        }
                    }
                }
                // if(indices.size()>0 && stra.s_regions[indices[indices.size()-1]].to==stra.size && stra.s_regions[indices[indices.size()-1]].size()>=k-1){
                // } else {
                //     if(prokrustean_optional.stratum_right_ext_count[i]>1){
                //         // divergence
                //         // cout << "divergence at stratum " << i << " (" << (int)prokrustean_optional.stratum_right_ext_count[i] << ")" << endl;  
                //         cnt+=prokrustean_optional.stratum_right_ext_count[i];
                //         stats[3]+=prokrustean_optional.stratum_right_ext_count[i];
                //     }
                // }
            } else {
                // special case
                if(prokrustean_optional.stratum_right_ext_count[i]>1){
                    // divergence multiple -> convergence does not matter
                    // cout << "divergence at stratum of k-1 " << i << " (" << (int)prokrustean_optional.stratum_right_ext_count[i] << ")" << endl;  
                    cnt+=prokrustean_optional.stratum_right_ext_count[i];
                    stats[4]+=prokrustean_optional.stratum_right_ext_count[i];
                    if(sequences!=nullptr){
                        cout << "stratum xxxx " << i << endl;
                        auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(get<1>(stratum_occ_samples[i])+spectrum[spectrum.size()-1].from, spectrum[spectrum.size()-1].size());
                        txt+="("+to_string(prokrustean_optional.stratum_right_ext_count[i])+")";
                        strings[4].push_back(txt);
                        // for(int i=0; i< prokrustean_optional.stratum_right_ext_count[i]; i++){
                        //     // strings[4].push_back("div");
                        // }
                    }
                } else if(prokrustean_optional.stratum_right_ext_count[i]==1){
                    if(prokrustean_optional.stratum_left_ext_count[i]==0){
                        // divergence single -> convergence can work
                        // cout << "tip at stratum of k-1 " << i << endl;  
                        cnt++;
                        stats[5]++;
                        if(sequences!=nullptr){
                            auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(get<1>(stratum_occ_samples[i])+spectrum[spectrum.size()-1].from, spectrum[spectrum.size()-1].size());
                            strings[5].push_back(txt);
                        }
                    } else if(prokrustean_optional.stratum_left_ext_count[i]>1){
                        // cout << "convergence at stratum of k-1 " << i << endl;  
                        cnt++;
                        stats[6]++;
                        if(sequences!=nullptr){
                            auto txt = (*sequences)[get<0>(stratum_occ_samples[i])].substr(spectrum[spectrum.size()-1].from, spectrum[spectrum.size()-1].size());
                            strings[6].push_back(txt);
                        }
                    }
                }
            }
        }
        cout << "maximal unitigs : " << cnt << endl;
        if(sequences==nullptr){
            for(int i=0; i<stats.size(); i++){
                if(i==0) cout << "tip at seq " << stats[i] << endl; 
                if(i==1) cout << "tip at stratum " << stats[i] << endl; 
                if(i==2) cout << "convergence at stratum " << stats[i] << endl; 
                if(i==3) cout << "divergnece at stratum " << stats[i] << endl; 
                if(i==4) cout << "divergnece at stratum of k-1 " << stats[i] << endl; 
                if(i==5) cout << "tip at stratum of k-1 " << stats[i] << endl; 
                if(i==6) cout << "convergence at stratum of k-1 " << stats[i] << endl; 
            }
        } else {
            for(int i=0; i<stats.size(); i++){
                if(i==0){
                    cout << "tip at seq " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                } 
                if(i==1) {
                    cout << "tip at stratum " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                }
                if(i==2){
                    cout << "convergence at stratum " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                } 
                if(i==3) {
                     cout << "divergnece at stratum " << stats[i] << endl; 
                     for(auto &s: strings[i]) cout << s << endl;
                }
                if(i==4) {
                    cout << "divergnece at stratum of k-1 " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                }
                if(i==5) {
                    cout << "tip at stratum of k-1 " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                }
                if(i==6) {
                    cout << "convergence at stratum of k-1 " << stats[i] << endl; 
                    for(auto &s: strings[i]) cout << s << endl;
                }
            }
        }
        return cnt;
    }
    void fill_stratum_left_bound_single_right_extension(int k){
        // for(int i=0; i<prokrustean.sequence_count(); i++){
        //     Sequence seq = prokrustean.get_sequence(i);
        //     seq.print();
        // }
        // for(int i=0; i<prokrustean.stratum_count(); i++){
        //     Stratum stra = prokrustean.get_stratum(i);
        //     stra.print();
        // }
        this->stratum__left_bound_stratum_is_branching=vector<optional<bool>>(prokrustean.stratum_count());
        vector<uint8_t> indices;
        /* single thread version of filling extensions */
        for(int i=0; i<prokrustean.sequence_count(); i++){
            Sequence seq = prokrustean.get_sequence(i);
            // seq.print();
            seq.get_valid_indices(k, indices);
            // cout<< "seq " << i << " indices " << indices.size() << endl;
            // has to pair two
            for(int r=0; r+1<indices.size(); r++){
                fill_for_two_(seq.s_regions[indices[r]], seq.s_regions[indices[r+1]], k);
            }
        }
        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            stra.print();
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
            if(seq.size<k){
                continue;
            }
            
            if(seq.s_regions.size()>0 && seq.s_regions[0].size()>=k){
                this->sequence__is_start_of_maximal_unitig[i]=false;
            } else {
                // tip location: k-reflectum will have single occurrence (in sequence) 
                cout << "seq tip " << i << " found" << endl;
                this->sequence__is_start_of_maximal_unitig[i]=true;
                cnt++;
            }
        }

        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            if(stra.size<k){
                continue;
            }
            if(stra.s_regions.size()>0 && stra.s_regions[0].size()>=k){
                this->stratum__is_start_of_maximal_unitig[i]=false;
            } else if(prokrustean_optional.stratum_left_ext_count[i]!=1) {
                // if count is 0, then tip location.
                // if count is >1, then branch convergence.
                this->stratum__is_start_of_maximal_unitig[i]=true;
                cnt++;
            } else if(!stratum__left_bound_stratum_is_branching[i].has_value()){
                // 1. maybe stratified region is on pos 0
                // 2. maybe ancestors of every occurrence is leftmost but not pos 0 -> always reflectum defines the maximal
                this->stratum__is_start_of_maximal_unitig[i]=false;
            } else if(stratum__left_bound_stratum_is_branching[i].value()){
                this->stratum__is_start_of_maximal_unitig[i]=true;
                cnt++;
            } else {
                this->stratum__is_start_of_maximal_unitig[i]=false;
            }
        }

        vector<Region> spectrum;
        for(int i=0; i<prokrustean.sequence_count(); i++){
            Sequence seq = prokrustean.get_sequence(i);
            prokrustean.get_spectrum(seq, k, spectrum);
            for(int r=0; r< spectrum.size(); r++){
                if(r>0 && spectrum[r].is_reflected){
                    // must be stratified
                    if(1<prokrustean_optional.stratum_right_ext_count[spectrum[r-1].stratum_id]){
                        cnt++;
                        cout << "found reflectum after branching stratum id: " << spectrum[r-1].stratum_id << " rgn: " << r << endl;
                    }
                }
            }
        }

        for(int i=0; i<prokrustean.stratum_count(); i++){
            Stratum stra = prokrustean.get_stratum(i);
            prokrustean.get_spectrum(stra, k, spectrum);
            for(int r=0; r< spectrum.size(); r++){
                if(r>0 && spectrum[r].is_reflected){
                    // must be stratified
                    if(1<prokrustean_optional.stratum_right_ext_count[spectrum[r-1].stratum_id]){
                        cnt++;
                        cout << "found reflectum after branching stratum id: " << spectrum[r-1].stratum_id << " rgn: " << r << endl;
                    }
                }
            }
        }

        cout << "maximal unitig counted: " << cnt << endl;
        ran__fill_stratum_left_bound_single_right_extension=true;
    }
    
    void fill_for_two_(StratifiedRegion& l_rgn, StratifiedRegion& r_rgn, int k){
        // go deepest to rightmost of left and leftmost of right.
        // if k is fixed, there is no way but go into the deepest stratum such that the rightmost/leftmost is not (stratified and >=k)
        StratumId l_id= get_deepest_stratum(l_rgn, k, true);
        // note: since it is left side region of a pair, right ext is always 1 or more
        bool left_desendant_is_branching=prokrustean_optional.stratum_right_ext_count[l_id]>1;
        cout << "left descendant: " << left_desendant_is_branching << endl;
        // StratumId r_id= get_deepest_stratum(r_rgn, k, true);
        // bool left_descendant_is_diverging=prokrustean_optional.stratum_right_ext_count[l_id]!=1;
        auto r_id=r_rgn.stratum_id;
        while(true){
            //already computed 
            if(stratum__left_bound_stratum_is_branching[r_id].has_value()){
                break;
            }
            auto child_stratum=prokrustean.get_stratum(r_id);
            //leftmost shares edge
            if(child_stratum.s_regions.size()>0){
                auto& leftmost_rgn = child_stratum.s_regions[0]; 
                if(leftmost_rgn.from==0 && leftmost_rgn.size()>=k){
                    r_id=leftmost_rgn.stratum_id;
                    continue;
                }
            }    
            stratum__left_bound_stratum_is_branching[r_id]=left_desendant_is_branching;
            break;
        }

    }

    StratumId get_deepest_stratum(StratifiedRegion& rgn, int k, bool is_direction_left){
        auto stratum_id=rgn.stratum_id;
        if(is_direction_left){
            while(true){
                auto child_stratum=prokrustean.get_stratum(stratum_id);
                //rightmost of the left shares edge
                if(child_stratum.s_regions.size()>0){
                    auto& rightmost_rgn = child_stratum.s_regions[child_stratum.s_regions.size()-1]; 
                    if(rightmost_rgn.to==child_stratum.size && rightmost_rgn.size()>=k){
                        stratum_id=rightmost_rgn.stratum_id;
                        continue;
                    }
                }
                // stop
                break;
            }
        } else {
            while(true){
                auto child_stratum=prokrustean.get_stratum(stratum_id);
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