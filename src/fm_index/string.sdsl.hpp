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
    wt_blcd<> wt;
public:
	WaveletString(){}

	/*
	 * constructor from ASCII file
	 */
	WaveletString(string path, char term='#'){
        construct(wt, path, std::is_rvalue_reference<string &&>::value);

        for(uint8_t c=0; c < std::numeric_limits<uint8_t>::max(); c++){
            if(wt.rank(wt.size(),c)>0) characters.push_back(c);
        }

        for(int i=0; i< characters.size(); i++){
            characters_inv[characters[i]] = i;
        }
        
        bool only_ACGTN = true;
        for(auto c: characters){
            if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N' && c!=term){
                cout << "-- characters(" << characters.size() << ") --" << endl;
                cout << "Warning: A symbol not in AGCTN and term is included. The code still runs but please check the data validity." << endl;
                for(auto c: characters){
                    cout << c << " ";
                }
                cout << endl;
            }
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

	//the order follows the character sequence
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

	virtual uint64_t size(){
        return wt.size();
    }
};

#endif /* FM_INDEX_RANK_HPP_ */
