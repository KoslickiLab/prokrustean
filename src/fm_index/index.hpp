

#include "rank.hpp"
#include <algorithm>

using namespace std;

#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

/*
 * factory. For potential optimization later.
 */
// FmIndex create_fm_index(string path, char TERM = '#'){
// 	FmIndex idx = FmIndex(path, TERM);
// 	return idx;
// }

struct CArray {
	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;
};

class FmIndex{

public:
	/*
	 * constructor path of a STRING file containing the STRING in ASCII format
	 */
	FmIndex(string path, char TERM = '#'){

		this->TERM = TERM;

        STRING = SuccintString(path, TERM);
		auto r = STRING.parallel_rank(STRING.size());
        C = {
			r.TERM, 
			r.TERM+r.A, 
			r.TERM+r.A+r.C, 
			r.TERM+r.A+r.C+r.G
		};
    }
	char TERM = '#';
    CArray C;
	SuccintString STRING;

    uint64_t LF(uint64_t r){
        uint64_t f;
        // cout << STRING[r] << endl;
        switch (STRING[r])
        {
            case 'A': f=C.A+STRING.rank(r, 'A'); break;
            case 'C': f=C.C+STRING.rank(r, 'C'); break;
            case 'G': f=C.G+STRING.rank(r, 'G'); break;
            case 'T': f=C.T+STRING.rank(r, 'T'); break;
            default: //TERM, otherwise invalidated already.
                f=STRING.rank(r, TERM);
                break;
        }
        return f;
    }
    
    char get_character(uint64_t pos){
        return STRING[pos];
    }

    uint64_t size(){
        return STRING.size();
    }

	uint64_t seq_cnt(){
        return C.A;
    }
};

#endif /* FM_INDEX_INDEX_HPP_ */
