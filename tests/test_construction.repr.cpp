#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"

using namespace std;

int C_CNT = 4;

bool _decide_repr_left(CharId c, vector<tuple<CharId, CharId>> distinct_extensions){
    auto decision = decide_repr_sa_extensions(C_CNT, distinct_extensions);
    auto left = get<0>(decision);
    for(auto l: left){
        if(l==c){
            return true;
        }
    }
    return false;
}

bool _decide_repr_right(CharId c, vector<tuple<CharId, CharId>> distinct_extensions){
    auto decision = decide_repr_sa_extensions(C_CNT, distinct_extensions);
    auto right = get<1>(decision);
    for(auto r: right){
        if(r==c){
            return true;
        }
    }
    return false;
}


bool _decide_repr_naive_left(CharId c, vector<tuple<CharId, CharId>> distinct_extensions){
    vector<tuple<CharId, CharId>> c_extensions;
    for(auto t: distinct_extensions){
        auto l = get<0>(t);
        if(l!=c) continue;
        
        c_extensions.push_back(t);
    }
    if(c_extensions.size()==0){
        return false;
    }
    //exclusive
    if(c_extensions.size()==1){
        auto r = get<1>(c_extensions[0]);
        if(r>0){ //r is not term 
            int l_cnt = 0;
            for(auto t: distinct_extensions){
                if(r==get<1>(t)) 
                l_cnt++;
            }
            // not exclusive to each other    
            if(l_cnt>1)
            return false;
        }
    }
    return true;
}

bool _decide_repr_naive_right(CharId c, vector<tuple<CharId, CharId>> distinct_extensions){
    vector<tuple<CharId, CharId>> c_extensions;
    for(auto t: distinct_extensions){
        auto r = get<1>(t);
        if(r!=c) continue;
        
        c_extensions.push_back(t);
    }
    if(c_extensions.size()==0){
        return false;
    }
    //exclusive
    if(c_extensions.size()==1){
        auto l = get<0>(c_extensions[0]);
        if(l>0){ //l is not term 
            int r_cnt = 0;
            for(auto t: distinct_extensions){
                if(l==get<0>(t)) 
                r_cnt++;
            }
            // not exclusive to each other    
            if(r_cnt>1)
            return false;
        }
    }
    return true;
}

vector<tuple<CharId, CharId>> _sample_exts(vector<bool> choices){
    vector<tuple<CharId, CharId>> exts;
    for(CharId l=0; l<C_CNT; l++){
        for(CharId r=0; r<C_CNT; r++){
            exts.push_back(make_tuple(l, r));
        }
    }
    vector<tuple<CharId, CharId>> filtered;
    for(int i=0; i<choices.size(); i++){
        if(choices[i]){
            filtered.push_back(exts[i]);
        }
    }
    return filtered;
}

vector<bool> convertToBinary(unsigned int n, int digits)
{
    vector<bool> bits;
    unsigned int no = n;
    while(digits>0){
        bits.push_back(no/2);
        no /= 2;
        digits--;
    }
    return bits;
}

void test_repr_rule(){
    int digits = C_CNT * C_CNT;
    int all_choices_permut = pow(2, digits);
    vector<vector<bool>> exustive_samples;
    for(int i=0; i<all_choices_permut; i++){
        vector<bool> choices = convertToBinary(i, digits);
        exustive_samples.push_back(choices);
    }
    for(auto sampling: exustive_samples){
        auto exts = _sample_exts(sampling);
        // cout << "-- ext --" << endl;
        // for(auto ext: exts){
        //     cout << (int)get<0>(ext) << " , " << (int)get<1>(ext) << endl;
        // }
        
        for(int c = 1; c < C_CNT; c++){
            // cout << "c: " << c << ", "; 
            // cout << (int)_decide_repr_left(c, exts) << " , " << (int)_decide_repr_naive_left(c, exts) << endl;
            assert(_decide_repr_left(c, exts)==_decide_repr_naive_left(c, exts));
            assert(_decide_repr_right(c, exts)==_decide_repr_naive_right(c, exts));
        }
    }
}


void main_construction_repr() {
    test_repr_rule();
}
