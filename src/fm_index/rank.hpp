// tbd licence

/*
 * rank.hpp
 *
 *  Created on: Jul 29, 2023
 *      Author: Adam Park
 * 
 *  Originally worked by Nicola Prezza
 *  https://github.com/nicolaprezza/bwt2lcp/blob/master/internal/dna_string.hpp
 * 
 *  Optimized string with rank on DNA alphabet: {A,C,G,T,TERM}
 *
 *  One access or a parallel rank for all 4 letters A,C,G,T causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^64
 *
 *  Supports very efficient (1 cache miss) parallel rank for (A,C,G,T), and (1 cache miss) single rank for TERM
 *
 *  Data is stored and cache-aligned in blocks of 512 bits (64 bytes)
 *
 *  Size of the string: 4n bits, where n = string length
 *
 *  memory layout:
 *
 *  | 32-bit rank A | 32-bit rank C | 32-bit rank G | 32-bit rank T|
 *  | 128-bit 1st bits | 128-bit 2nd bits | 128-bit 3rd bits |
 *
 *  careful: __uint128_t inverts the first and last 64 bits of the number!
 */

#ifndef FM_INDEX_RANK_HPP_
#define FM_INDEX_RANK_HPP_

#define SUPERBLOCK_SIZE 0x100000000 	//number of characters in a superblock = 2^32 characters
#define BLOCKS_PER_SUPERBLOCK 33554432	//blocks in a superblock
#define BYTES_PER_SUPERBLOCK 2147483648	//bytes in a superblock
#define BLOCK_SIZE 128 					//number of characters inside a block
#define BYTES_PER_BLOCK 64				//bytes in a block of 512 bits
#define ALN 64							//alignment

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>

using namespace std;

struct ParallelRank{

public:

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;
	uint64_t TERM = 0; // filled later

	void fill_term(uint64_t i){
		TERM = i - (A + C + G + T);
	}
	
	ParallelRank operator+(const ParallelRank& a) const{

		return {
			a.A + A,
			a.C + C,
			a.G + G,
			a.T + T
		};

	}

	bool operator==(const ParallelRank& a) const{

		return a.A == A and a.C == C and a.G == G and a.T == T;

	}

	bool operator!=(const ParallelRank& a) const{

		return a.A != A or a.C != C or a.G != G or a.T != T;

	}

	bool operator<=(const ParallelRank& a) const{

		return A <= a.A and C <= a.C and G <= a.G and T <= a.T;

	}

};

struct CArray {
	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;
};

class SuccintString{

public:

	SuccintString(){}

	/*
	 * constructor from ASCII file
	 */
	SuccintString(string path, char TERM = '#'){

		this->TERM = TERM;
		c_array = {0,0,0,0};

		n = uint64_t(filesize(path));

		n_superblocks = (n+1)/SUPERBLOCK_SIZE + ((n+1)%SUPERBLOCK_SIZE != 0);
		n_blocks = (n+1)/BLOCK_SIZE + ((n+1)%BLOCK_SIZE != 0);
		nbytes = (n_blocks * BYTES_PER_BLOCK);//number of bytes effectively filled with data

		superblock_ranks = vector<ParallelRank>(n_superblocks);

		/*
		 * this block of code ensures that data is aligned by 64 bytes = 512 bits
		 */
		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;

		//cout << "alignment of data: " << (void*)data << endl;
		
		ifstream ifs(path);

		string BUF(BLOCK_SIZE,'A');

		for(uint64_t i = 0; i<n; ++i){

			char c;

			ifs.read((char*)&c, sizeof(char));

			BUF[i%BLOCK_SIZE] = c;

			//buffer is full
			if((i%BLOCK_SIZE) == BLOCK_SIZE-1) set(i/BLOCK_SIZE, BUF);

			switch (c)
			{
				case 'A': c_array.C++; break;
				case 'C': c_array.G++; break;
				case 'G': c_array.T++; break;
				case 'T': break;
				default:
				if(c==TERM){
					c_array.A++;
				} else {
					cout << "Error while reading file: read forbidden character '" <<  c << "' (ASCII code " << int(c) << ")." << endl;
					cout << "Only A,C,G,T, and " << TERM << " are admitted in the input BWT!" << endl;
					// if(c=='N'){
					// 	cout << "Possible solution: it seems that your file contains 'N' characters. Please, re-run with option -n." << endl;
					// }else{
					// 	cout << "Possible solution: if the unknown character is the terminator, you can solve the problem by adding option \"-t " << int(c) << "\"." << endl;
					// }
					exit(1);
				}
				break;
			}
		}
		c_array.C += c_array.A;
		c_array.G += c_array.C;
		c_array.T += c_array.G;

		if(n % BLOCK_SIZE != 0) set(n/BLOCK_SIZE, BUF);

		

		// assert(check_content(path));
		build_rank_support();

	}

	std::ifstream::pos_type filesize(string filename){
		std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
		return in.tellg();
	}

	//return i-th character
	char operator[](uint64_t i){

		assert(i<n);

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		//chars[0..3] contains 1st, 2nd, 3rd least significant bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		uint64_t b =	((chars[0]>>(BLOCK_SIZE-(block_off+1)))&0x1) +
						(((chars[1]>>(BLOCK_SIZE-(block_off+1)))&0x1)<<1) +
						(((chars[2]>>(BLOCK_SIZE-(block_off+1)))&0x1)<<2);

		return 	(b == 0)*'A' +
				(b == 1)*'C' +
				(b == 2)*'G' +
				(b == 3)*'T' +
				(b == 4)*TERM;

	}

	/*
	 * Parallel rank of (A,C,T,G) at position i.
	 */
	ParallelRank parallel_rank(uint64_t i){

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		ParallelRank superblock_r = superblock_ranks[superblock_number];
		ParallelRank block_r = get_counters(superblock_number,block_number);

		ParallelRank r = superblock_r + block_r + block_rank(superblock_number, block_number, block_off);
		r.fill_term(i);
		return r;
	}

	/*
	 * standard rank. c can be A,C,G,T, or TERM
	 * todo: optimize?
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		if(c==TERM) return rank_non_dna(i);
		
		ParallelRank pr = parallel_rank(i);

		switch(c){
			case 'A' : return pr.A; break;
			case 'C' : return pr.C; break;
			case 'G' : return pr.G; break;
			case 'T' : return pr.T; break;
		}

		return 0;
	}

	/*
	 * return number of non-dna symbols in the prefix of length i of the text. At most 1 cache miss!
	 */
	uint64_t rank_non_dna(uint64_t i){

		assert(i<=n);
		auto r = parallel_rank(i);

		assert(r.A + r.C + r.G + r.T <= i);

		return i - (r.A + r.C + r.G + r.T);

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&nbytes,sizeof(nbytes));
		out.write((char*)&n_superblocks,sizeof(n_superblocks));
		out.write((char*)&n_blocks,sizeof(n_blocks));

		w_bytes += sizeof(n) + sizeof(nbytes) + sizeof(n_superblocks) + sizeof(n_blocks);

		out.write((char*)superblock_ranks.data(),n_superblocks*sizeof(ParallelRank));
		w_bytes += n_superblocks*sizeof(ParallelRank);

		out.write((char*)data,nbytes*sizeof(uint8_t));
		w_bytes += nbytes*sizeof(uint8_t);

		return w_bytes;

	}

	CArray get_count_array(){
		return c_array;
	}
	// void load(std::istream& in) {

	// 	in.read((char*)&n,sizeof(n));
	// 	in.read((char*)&nbytes,sizeof(nbytes));
	// 	in.read((char*)&n_superblocks,sizeof(n_superblocks));
	// 	in.read((char*)&n_blocks,sizeof(n_blocks));

	// 	superblock_ranks = vector<p_rank>(n_superblocks);
	// 	in.read((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank));

	// 	memory = vector<uint8_t>(nbytes+ALN,0);
	// 	data = memory.data();
	// 	while(uint64_t(data) % ALN != 0) data++;
	// 	in.read((char*)data,nbytes*sizeof(uint8_t));

	// 	assert(check_rank());

	// }

	uint64_t size(){
		return n;
	}

private:

	void build_rank_support(){

		ParallelRank superblock_r = {};
		ParallelRank block_r = {};

		for(uint64_t bl = 0; bl < n_blocks-1; ++bl){

			uint64_t superblock_number = bl/BLOCKS_PER_SUPERBLOCK;
			uint64_t block_number = bl%BLOCKS_PER_SUPERBLOCK;

			if(block_number == 0){

				superblock_ranks[superblock_number]=superblock_r;
				block_r = {};

			}

			set_counters(superblock_number, block_number,block_r);

			ParallelRank local_rank = block_rank(superblock_number, block_number);

			block_r = block_r + local_rank;
			superblock_r = superblock_r + local_rank;

		}

		uint64_t superblock_number = (n_blocks-1)/BLOCKS_PER_SUPERBLOCK;
		uint64_t block_number = (n_blocks-1)%BLOCKS_PER_SUPERBLOCK;

		if(block_number == 0){

			superblock_ranks[superblock_number]=superblock_r;
			block_r = {};

		}

		set_counters(superblock_number, block_number,block_r);

		// assert(check_rank());

	}

	/*
	 * set i-th block to s. Assumption: s.length() == BLOCK_SIZE
	 */
	void set(uint64_t i, string & s){

		assert(s.length()==BLOCK_SIZE);
		assert(i<n_blocks);

		uint64_t superblock_number = i / BLOCKS_PER_SUPERBLOCK;
		uint64_t block_number = i % BLOCKS_PER_SUPERBLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		chars[2] = 0;
		chars[1] = 0;
		chars[0] = 0;

		for(auto c : s){

			if(c==TERM){

				chars[2] = (chars[2] << 1) | __uint128_t(0x1);
				chars[1] = chars[1] << 1;
				chars[0] = chars[0] << 1;

			}else{

				switch(c){

				case 'A' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = chars[0] << 1;break;

				case 'C' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				case 'G' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = chars[0] << 1;break;

				case 'T' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				}

			}

		}

	}


	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline ParallelRank block_rank(uint64_t superblock_number, uint64_t block_number, uint64_t block_off){

		assert(block_off<BLOCK_SIZE);

		if(block_off < 64) return block_rank64(superblock_number, block_number, block_off);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(start);

		__uint128_t PAD = ((~__uint128_t(0))>>block_off);

		__uint128_t b = ~(chars[2] | PAD); //most significant bit, padded and negated

		__uint128_t b0 = b & (~chars[1]);
		__uint128_t b1 = b & chars[1];

		return {

			popcount128(b0 & (~chars[0])),
			popcount128(b0 & (chars[0])),
			popcount128(b1 & (~chars[0])),
			popcount128(b1 & (chars[0]))

		};

	}

    inline uint64_t popcount128(__uint128_t x){

	    return __builtin_popcountll(uint64_t(x>>64)) + __builtin_popcountll( x & 0xFFFFFFFFFFFFFFFF );

    }

	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline ParallelRank block_rank64(uint64_t superblock_number, uint64_t block_number, uint64_t block_off){

		assert(block_off < 64);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0,2,4] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		uint64_t* chars = (uint64_t*)(start);

		uint64_t PAD = ((~uint64_t(0))>>block_off);

		uint64_t b = ~(chars[5] | PAD); //most significant bit, padded and negated

		uint64_t b0 = b & (~chars[3]);
		uint64_t b1 = b & chars[3];

		return {

			(uint64_t)__builtin_popcountll(b0 & (~chars[1])),
			(uint64_t)__builtin_popcountll(b0 & (chars[1])),
			(uint64_t)__builtin_popcountll(b1 & (~chars[1])),
			(uint64_t)__builtin_popcountll(b1 & (chars[1]))

		};

	}

	/*
	 * rank in whole block given as coordinates: superblock, block
	 */
	inline ParallelRank block_rank(uint64_t superblock_number, uint64_t block_number){

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(start);

		__uint128_t b2 = ~chars[2]; //most significant bit

		ParallelRank res = {

			popcount128(b2 & (~chars[1]) & (~chars[0])),
			popcount128(b2 & (~chars[1]) & (chars[0])),
			popcount128(b2 & (chars[1]) & (~chars[0])),
			popcount128(b2 & (chars[1]) & (chars[0]))

		};

		//assert(check_rank_local(res, superblock_number*SUPERBLOCK_SIZE + block_number*BLOCK_SIZE, block_off));

		return res;

	}

	/*
	 * set counters in the i-th block to r
	 */
	void set_counters(uint64_t superblock_number, uint64_t block_number, ParallelRank r){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		block_ranks[0] = r.A;
		block_ranks[1] = r.C;
		block_ranks[2] = r.G;
		block_ranks[3] = r.T;

		assert(get_counters(superblock_number,block_number) == r);

	}

	/*
	 * get counters of the i-th block
	 */
	inline ParallelRank get_counters(uint64_t superblock_number, uint64_t superblock_off){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + superblock_off*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		return {
			block_ranks[0],
			block_ranks[1],
			block_ranks[2],
			block_ranks[3]
		};

	}

	char TERM = '#';

	uint64_t n_superblocks = 0;
	uint64_t n_blocks = 0;

	vector<uint8_t> memory; //allocated memory

	//data aligned with blocks of 64 bytes = 512 bits
	uint8_t * data = NULL;

	vector<ParallelRank> superblock_ranks;

	uint64_t nbytes = 0; //bytes used in data
	uint64_t n = 0;

	CArray c_array;

};


#endif /* FM_INDEX_RANK_HPP_ */
