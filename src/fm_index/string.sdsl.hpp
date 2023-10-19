// tbd licence

/*
 * string.hpp
 *
 *  Created on: Aug 15, 2023
 *      Author: Adam Park
 * 
 *  Originally worked by sdsl
 *  
 * A simple 
 */

#ifndef FM_INDEX_RANK_HPP_
#define FM_INDEX_RANK_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "index.hpp"
#include "../sdsl/wavelet_trees.hpp"

using namespace std;
using namespace sdsl;

class WaveletString: public AbstractString{
private:
    vector<char> characters;
    map<char, uint8_t> characters_inv;
    int characters_cnt;
    // wt_blcd<> wt;
public:
    wt_blcd<> wt;
	WaveletString(){}

	/*
	 * constructor from ASCII file
	 */
	WaveletString(string path, char term='$'){
        construct(wt, path, std::is_rvalue_reference<string &&>::value);
        
        for(uint8_t c=0; c < std::numeric_limits<uint8_t>::max(); c++){
            if(wt.rank(wt.size(),c)>0) characters.push_back(c);
        }
        characters_cnt=characters.size();
        for(int i=0; i< characters.size(); i++){
            characters_inv[characters[i]] = i;
        }
        
        bool only_ACGTN = true;
        bool has_term = false;
        for(auto c: characters){
            if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N' && c!=term){
                cout << "warning: A symbol (" << c <<") not in AGCTN and term("<< term <<")." << endl;
                only_ACGTN = false;
            }
            if(c==term){
                has_term=true;
            }
        }
        if(!has_term){
            cout << "the termination ("<< term <<") does not exist." << endl;
            assert(has_term);
        }
    }

    vector<char> get_characters(){
        return characters;
    };
	
    char operator[](uint64_t i){
        return wt[i];
    }

	CharId access(uint64_t i){
        return characters_inv[wt[i]];
    }

	uint64_t rank(SuffixArrayIdx i, CharId c){
        return wt.rank(i, characters[c]);
    }

    uint64_t rank(SuffixArrayIdx i, char c){
        return wt.rank(i, c);
    }

    void ranks(CharId c, vector<SuffixArrayIdx> &firsts, vector<uint64_t> &ext_ranks){
        /* check interval first, so that if empty, just skip it*/
        ext_ranks[0]=wt.rank(firsts[0], characters[c]);
        ext_ranks[characters_cnt]=wt.rank(firsts[characters_cnt], characters[c]);
        if(ext_ranks[0] == ext_ranks[characters_cnt]){
            for(int i=1; i< characters_cnt; i++){
                ext_ranks[i]=ext_ranks[0];
            }
            return;
        }
        // do the rest
        for(int i=1; i< characters_cnt; i++){
            if(firsts[i-1]!=firsts[i]){
                ext_ranks[i]=wt.rank(firsts[i], characters[c]);
            } else {
                ext_ranks[i]=ext_ranks[i-1];
            }
        }
    }

    RankArray ranks(SuffixArrayIdx i){
        RankArray ranks;
        for(auto c: characters){
            ranks.push_back(wt.rank(i, c));
        }
        return ranks;
    }

    uint64_t select(uint64_t i, CharId c){
        return wt.select(i, characters[c]);
    }

    uint64_t select_by_char(uint64_t i, char c){
        return wt.select(i, c);
    }

	uint64_t size(){
        return wt.size();
    }

    void dispose(){
        wt=wt_blcd<>();
    }
};

#endif /* FM_INDEX_RANK_HPP_ */
