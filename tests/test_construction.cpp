#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"

using namespace std;

int C_CNT = 4;

optional<SuffixArrayIdx> _get_repr_right(CharId c, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
    optional<SuffixArrayIdx> repr_idx = nullopt;
    for(auto t: distinct_extensions){
        auto r = get<1>(t);
        if(r!=c) continue;
        
        auto sa_idx = get<2>(t);
        repr_idx = repr_idx.has_value()? min(repr_idx.value(), sa_idx): sa_idx; 
    }
    if(repr_idx.has_value()){
        vector<SuffixArrayIdx> sa_indexes = decide_repr_sa_extensions(C_CNT, distinct_extensions);
        bool exists = false;
        for(auto sa_idx: sa_indexes){
            exists = exists? exists : repr_idx.value() == sa_idx;
        }
        if(!exists){
            repr_idx = nullopt;
        }
    }
    return repr_idx;
}

optional<SuffixArrayIdx> _get_repr_left(CharId c, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
    optional<SuffixArrayIdx> repr_idx = nullopt;
    for(auto t: distinct_extensions){
        auto l = get<0>(t);
        if(l!=c) continue;
        
        auto sa_idx = get<2>(t);
        repr_idx = repr_idx.has_value()? min(repr_idx.value(), sa_idx): sa_idx; 
    }
    if(repr_idx.has_value()){
        vector<SuffixArrayIdx> sa_indexes = decide_repr_sa_extensions(C_CNT, distinct_extensions);
        bool exists = false;
        for(auto sa_idx: sa_indexes){
            exists = exists? exists : repr_idx.value() == sa_idx;
        }
        if(!exists){
            repr_idx = nullopt;
        }
    }
    return repr_idx;
}

optional<SuffixArrayIdx> _get_repr_naive_right(CharId c, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
    optional<SuffixArrayIdx> repr_idx = nullopt;
    vector<tuple<CharId, CharId, SuffixArrayIdx>> c_extensions;
    for(auto t: distinct_extensions){
        auto r = get<1>(t);
        if(r!=c) continue;
        
        c_extensions.push_back(t);
        auto sa_idx = get<2>(t);
        repr_idx = repr_idx.has_value()? min(repr_idx.value(), sa_idx): sa_idx;
    }
    //exclusive
    if(c_extensions.size()==1){
        auto l = get<0>(c_extensions[0]);
        if(l>0){ //l is not term 
            int r_cnt = 0;
            for(auto t: distinct_extensions){
                if(l==get<0>(t)) r_cnt++;
            }
            // not exclusive to each other    
            if(r_cnt>1) repr_idx = nullopt;
        }
    }
    return repr_idx;
}

optional<SuffixArrayIdx> _get_repr_naive_left(CharId c, vector<tuple<CharId, CharId, SuffixArrayIdx>> distinct_extensions){
    optional<SuffixArrayIdx> repr_idx = nullopt;
    vector<tuple<CharId, CharId, SuffixArrayIdx>> c_extensions;
    for(auto t: distinct_extensions){
        auto l = get<0>(t);
        if(l!=c) continue;

        c_extensions.push_back(t);
        auto sa_idx = get<2>(t);
        repr_idx = repr_idx.has_value()? min(repr_idx.value(), sa_idx): sa_idx;
    }

    //some r is exclusive to l
    if(c_extensions.size()==1){
        auto r = get<1>(c_extensions[0]);
        if(r>0){ //r is not term 
            int l_cnt = 0;
            for(auto t: distinct_extensions){
                if(r==get<1>(t)) l_cnt++;
            }
            // not exclusive to each other    
            if(l_cnt>1) repr_idx = nullopt;
        }
    }
    return repr_idx;
}

vector<tuple<CharId, CharId, SuffixArrayIdx>> _sample_exts(vector<bool> choices){
    vector<tuple<CharId, CharId, SuffixArrayIdx>> exts;
    int sa_idx = 0;
    for(CharId l=0; l<C_CNT; l++){
        for(CharId r=0; r<C_CNT; r++){
            exts.push_back(make_tuple(l, r, sa_idx));
            sa_idx++;
        }
    }
    vector<tuple<CharId, CharId, SuffixArrayIdx>> filtered;
    for(int i=0; i<choices.size(); i++){
        if(choices[i]){
            filtered.push_back(exts[i]);
        }
    }
    return filtered;
}

void test_repr_rule(){
    int all_choices_cnt = C_CNT*C_CNT;
    int all_choices_permut = pow(2, all_choices_cnt);
    for(int i=0; i<all_choices_permut; i++){
        vector<bool> sampling;
        for(int p=0; p<all_choices_cnt; p++){
            int bit_pos = pow(2,p);
            sampling.push_back(i%bit_pos==0);
        }
        auto exts = _sample_exts(sampling);
        // cout << "-- ext --" << endl;
        // for(auto ext: exts){
        //     cout << (int)get<0>(ext) << " , " << (int)get<1>(ext) << " , " << get<2>(ext) << endl;
        // }

        for(int c = 1; c < C_CNT; c++){
            // cout << "c: " << c << ", "; 
            // if(_get_repr_left(c, exts).has_value()) cout << _get_repr_left(c, exts).value();
            // cout << " , ";
            // if(_get_repr_naive_left(c, exts).has_value()) cout << _get_repr_naive_left(c, exts).value();
            // cout << endl;
            assert(_get_repr_left(c, exts)==_get_repr_naive_left(c, exts));
            assert(_get_repr_right(c, exts)==_get_repr_naive_right(c, exts));
        }
    }
}


void main_construction() {
    test_repr_rule();
}
