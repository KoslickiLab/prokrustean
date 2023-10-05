#ifndef PROKRUSTEAN_HPP_
#define PROKRUSTEAN_HPP_

#include <stack>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "base.hpp"

using namespace std;

struct Prokrustean {
    /* data structure is succinct as possible*/
    vector<SequenceSize> sequences__size;
    vector<StratifiedData*> sequences__region;
    vector<uint8_t> sequences__region_cnt;
    
    vector<StratumSize> stratums__size;
    vector<StratifiedData*> stratums__region;
    vector<uint8_t> stratums__region_cnt;

    vector<tuple<SeqId, Pos>> stratum_occ_samples;

    uint64_t sequence_count(){
        return this->sequences__size.size();
    }

    uint64_t stratum_count(){
        return this->stratums__size.size();
    }

    void setup_stratum_example_occ(){
        stratum_occ_samples.resize(stratum_count());
        vector<bool> visits(stratum_count());
        std::stack<Stratum> stratum_stack;
        for(int i=0; i<sequence_count(); i++){
            for(auto &rgn: get_sequence(i).s_regions){
                auto stratum = get_stratum(rgn.stratum_id);
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
                        auto c_stratum = get_stratum(c_rgn.stratum_id);
                        c_stratum.set_occ(seq_id, rel_pos);
                        stratum_occ_samples[c_rgn.stratum_id]=c_stratum.example_occ.value();
                        visits[c_rgn.stratum_id]=true;
                        stratum_stack.push(c_stratum);
                    }
                }
            }
        }
        // for(auto &el: stratum_pos){
        //     cout << get<0>(el) << ", " << get<1>(el) << endl;
        // }

    }
    Sequence get_sequence(SeqId id){
        auto sequence=Sequence(id, sequences__size[id], sequences__region[id], sequences__region_cnt[id], stratums__size);
        sequence.set_occ(id, 0);
        assert(sequence.is_sequence);
        return sequence;
    }

    Stratum get_stratum(StratumId id){
        auto stratum=Stratum(id, stratums__size[id], stratums__region[id], stratums__region_cnt[id], stratums__size);
        if(stratum_occ_samples.size()>0){
            auto occ=stratum_occ_samples[id];
            stratum.example_occ=occ;
        }
        assert(stratum.is_stratum);
        return stratum;
    }

    void get_spectrum(Vertex v, int k,vector<Region> &output){
        output.clear();
        if(v.size<k){
            return;
        }
        // cout << "---- spectrum of -----" << endl;
        // v.print();
        vector<StratifiedRegion> regions_least_k;
        for(auto &r:v.s_regions){
            if(r.size()>=k){
                regions_least_k.push_back(r);
            }
        }
        int cnt=regions_least_k.size();
        // single reflected 
        if(cnt==0){
            // cout << "single reflected" << endl;
            auto rgn=ReflectedRegion(0, v.size);
            output.push_back(rgn);
            assert(rgn.is_reflected&&!rgn.is_stratified);
            return;
        }

        // regions are already sorted
        for(int i=0; i<cnt; i++){
            if(i==0 && regions_least_k[i].from>0){
                output.push_back(ReflectedRegion(0, regions_least_k[i].from+(k-1)));  
            }

            output.push_back(regions_least_k[i]);

            if(i<cnt-1){
                if(regions_least_k[i].to - regions_least_k[i+1].from >= k-1){
                } else {
                    output.push_back(ReflectedRegion(regions_least_k[i].to-(k-1), regions_least_k[i+1].from+(k-1)));          
                }
            }

            if(i==cnt-1 && regions_least_k[i].to<v.size){
                output.push_back(ReflectedRegion(regions_least_k[i].to-(k-1), v.size));  
            }
        }
    } 
    
    uint64_t get_cardinality(){
        uint64_t cnt = 0;
        for(auto rcnt: sequences__region_cnt){
            cnt+=rcnt;
        }
        for(auto rcnt: stratums__region_cnt){
            cnt+=rcnt;
        }
        return cnt;
    }


    void print_stratum(uint64_t id, vector<string> &seq_texts){
        auto stratum = get_stratum(id);
        auto seq_id=get<0>(stratum_occ_samples[id]);
        auto pos=get<1>(stratum_occ_samples[id]);
        auto str = seq_texts[seq_id];
        if(str.substr(pos, stratums__size[id]).find("GGC") == std::string::npos){
            return;
        }
        cout << "stratum("<< id <<"): " << str.substr(pos, stratums__size[id]) << endl; 
        for(auto &rgn: stratum.s_regions){
            cout << "stratified: " << str.substr(pos+rgn.from, rgn.size()) << endl; 
        }
    }
};

struct ProkrusteanOptional{
    bool collect_ext=false;
    bool collect_occ_samples=false;

    vector<uint8_t> stratum_left_ext_count;
    vector<uint8_t> stratum_right_ext_count;

    vector<tuple<SeqId, Pos>> stratum_occ_samples;

    ProkrusteanOptional(bool collect_ext){
        this->collect_ext=collect_ext;
    }
};
#endif