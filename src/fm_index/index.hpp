#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

#include <algorithm>
#include <ranges>

using namespace std;

//character ids
typedef uint8_t CharId;
//C array of fm index. Each index is CharId 
typedef vector<uint64_t> CArray;
typedef vector<uint64_t> RankArray;

class AbstractString{
public:
	CharId term_id=0;
	//the order follows the character sequence
	virtual vector<char> get_characters() = 0;
	virtual char operator[](uint64_t i) = 0;
	virtual CharId access(uint64_t i) = 0;
	virtual uint64_t rank(uint64_t i, CharId c) = 0;
	//the order follows the character sequence
    virtual RankArray ranks(uint64_t i) = 0;
	virtual uint64_t size() = 0;
};

class FmIndex{

public:
	/*
	 * constructor path of a STRING file containing the STRING in ASCII format
	 */
	FmIndex(AbstractString &string){
		this->STRING = &string;
        this->C = get_c_array(string);
    }
	char TERM='#'; //Lexicographically first
	CharId term_id=0;
    CArray C;
	AbstractString* STRING;

    uint64_t LF(uint64_t r){
        // cout << STRING[r] << endl;
		CharId cid = STRING->access(r);
		uint64_t f = C[cid]+STRING->rank(r, cid);
        return f;
    }
    
    char get_character(uint64_t pos){
        return (*STRING)[pos];
    }

    uint64_t size(){
        return STRING->size();
    }

	uint64_t seq_cnt(){
        return C[1];
    }

private:
	CArray get_c_array(AbstractString &string){
		RankArray r = string.ranks(string.size());
		CArray c_array(r.size());
		for(int i=0; i<c_array.size(); i++){
			for(int j=0; j<i; j++){
				c_array[i] += r[j];
			}
		}
		return c_array;
	}
};

#endif /* FM_INDEX_INDEX_HPP_ */
