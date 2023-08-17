#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/fm_index/tree.hpp"
#include "../src/fm_index/locate.hpp"

using namespace std;

vector<char> characters = {'A', 'C', 'G', 'T'};

int _occurrence(vector<string> sequences, string str){
    int occ = 0;
    for(auto seq: sequences){
        int pos = 0;
        while(pos+str.size()<seq.size()){
            if(str == seq.substr(pos, str.size())){
                occ++;
            }
            pos++;
        }
    }
    return occ;
}

set<string> _distinct(vector<string> sequences, int L){
    set<string> strings;
    for(auto seq: sequences){
        if(seq.size()<L){
            continue;
        }
        int pos = 0;
        while(pos < seq.size()-L){
            strings.insert(seq.substr(pos, L));
            pos++;    
        }
    }
    return strings;
}

vector<string> _find_maximal_repeats_naive(vector<string> sequences, int Lmin){
    int Lmax=0;
    for(auto seq:sequences){
        if(Lmax<seq.size()) 
        Lmax=seq.size();
    }
    
    vector<string> repeats;
    for(int L=Lmin; L<=Lmax; L++){
        auto dist_strings = _distinct(sequences, L);
        for(auto str: dist_strings){
            bool repeat = true;
            bool l_maximal = true;
            bool r_maximal = true;
            for(auto c: characters){
                int occ = _occurrence(sequences, str);
                int l_occ = _occurrence(sequences, c+str);
                int r_occ = _occurrence(sequences, str+c);
                // extension
                if(occ < 2) repeat = false;
                if(l_occ == occ) l_maximal = false;
                if(r_occ == occ) r_maximal = false;
            }
            if(repeat && l_maximal && r_maximal){
                repeats.push_back(str);
            }
        }
    }
    sort(repeats.begin(), repeats.end());
    return repeats;
}

void test_maximal_repeat(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    auto sequences = recover_text(fm_idx);
    auto repeats_naive = _find_maximal_repeats_naive(sequences, Lmin);

    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> rep_annot = navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx);
    auto sa = recover_suffix_array(fm_idx);
    set<string> uniq_repeats;
    for(auto r: rep_annot){
        // cout << "r size: " << r.size << endl;
        for(auto sa_idx: r.repr_indexes){
            auto str = sa[sa_idx].substr(0, r.size);
            // cout << "idx: " << sa_idx << " str: " << str << endl;
            uniq_repeats.insert(str);
        }
    }
    vector<string> repeats(uniq_repeats.begin(), uniq_repeats.end());
    sort(repeats.begin(), repeats.end());

    assert(repeats_naive.size() == repeats.size());
    for(int i=0; i<repeats.size(); i++){
        assert(repeats_naive[i]==repeats[i]);
    }
}


void main_construction_max_rep() {
    test_maximal_repeat();
}
