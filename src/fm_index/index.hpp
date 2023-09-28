#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

#include <algorithm>
#include <ranges>
#include "../prokrustean.hpp"

using namespace std;

//character ids -> lexicographical order of characters 
typedef uint8_t CharId;
typedef uint64_t SuffixArrayIdx;
//C array of fm index. Each index is CharId 
typedef vector<uint64_t> RankArray;

class AbstractString{
public:
	CharId term_id=0;
	//the order follows the character sequence
	virtual vector<char> get_characters() = 0;
	virtual char operator[](uint64_t i) = 0;
	virtual CharId access(uint64_t i) = 0;
	virtual uint64_t rank(SuffixArrayIdx i, CharId c) = 0;
	//the order follows the character sequence
    virtual RankArray ranks(SuffixArrayIdx i) = 0;
	virtual void ranks(CharId c, vector<SuffixArrayIdx> &firsts, vector<uint64_t> &ext_ranks)=0;
	virtual uint64_t select(uint64_t i, CharId c) = 0;
	virtual uint64_t select_by_char(uint64_t i, char c) = 0;
	virtual uint64_t size() = 0;
	virtual void dispose()=0;
};

class AbstractLocator{
public:
	virtual tuple<SeqId, Pos> get_location(SuffixArrayIdx suffix_idx)=0;
};

class FmIndex{

public:
	/*
	 * constructor path of a STRING file containing the STRING in ASCII format
	 */
	FmIndex(AbstractString &string, char TERM='$'){
		this->STRING = &string;
		this->characters_cnt = STRING->get_characters().size();
		this->TERM = TERM;

        this->C = get_c_array(string);
		this->char_abundances = get_abundences(this->C);
		this->characters_ranked_by_abundance = get_char_ids_by_abundence_desc(this->char_abundances);
		cout << "total length: " << this->STRING->size() << endl;
		cout << "abundences:";
		for(int c=0; c<characters_cnt; c++){
			cout << " " << STRING->get_characters()[c]<< ":" << this->char_abundances[c];
		}
		cout << endl;
    }
	char TERM='$'; //Lexicographically first
	CharId term_id=0;
    vector<uint64_t> C;
	AbstractString* STRING;
	AbstractLocator* locator=nullptr;
	// vector<char> characters;
	int characters_cnt;
	vector<uint64_t> char_abundances;
	vector<CharId> characters_ranked_by_abundance; // rank by abundances

    uint64_t LF(uint64_t r){
        // cout << STRING[r] << endl;
		CharId cid = STRING->access(r);
		uint64_t f = C[cid]+STRING->rank(r, cid);
        return f;
    }
    
	void set_sampled_suffix_array(AbstractLocator &locator){
		this->locator=&locator;
	}

    uint64_t size(){
        return STRING->size();
    }

	uint64_t seq_cnt(){
        return C[1];
    }

	void dispose(){
		this->STRING->dispose();
		this->STRING=nullptr;
	}

	string recover_text(int seq_id){
		uint64_t L = seq_id;
		uint64_t F = this->LF(L);

		string seq;
		while(F >= this->seq_cnt()){
			seq = (*this->STRING)[L] + seq;
			L = F;
			F = this->LF(L);
		}
		// must be terminator symbol
		// seq += (*this->STRING)[L]; 
		return seq;
	}

	void recover_all_texts(vector<string> &seqs){
		for(int i=0; i<this->seq_cnt(); i++){
			seqs.push_back(recover_text(i));
		}
	}

	/*
	debugging purpose
	index is position in the sequence,
	first: sa index
	second: actual suffix
	*/
	vector<pair<uint64_t, string>> recover_suffix_array(FmIndex &fm_idx, int seq_no, bool with_term=true){
		uint64_t L = seq_no;
		uint64_t F = this->LF(L);

		string seq(1, this->TERM);
		vector<pair<uint64_t, string>> sa;
		while(F >= this->seq_cnt()){
			seq = (*this->STRING)[L] + seq;
			sa.push_back(make_pair(F, seq));
			L = F;
			F = this->LF(L);
		}
		// important: this can be misleading because F is randomly (in lexicographical order) chosen
		if(with_term){
			string term(1, this->TERM);
			sa.insert(sa.begin(), make_pair(F, term));
		}
		reverse(sa.begin(), sa.end());

		return sa;
	}

	/*
	debugging purpose
	index is sa index
	*/
	vector<string> recover_suffix_array(FmIndex &fm_idx, bool with_term=true){
		vector<pair<uint64_t, string>> sa;
		for(int i=0; i<this->seq_cnt(); i++){
			for(auto pair: recover_suffix_array(fm_idx, i, with_term)){
				sa.push_back(pair);
			}
		}
		std::sort(sa.begin(), sa.end(), 
			[](tuple<int, string> const &t1, tuple<int, string> const &t2) {
				return get<0>(t1) < get<0>(t2); 
			});
		
		vector<string> suffixes;
		for(auto pair: sa){
			suffixes.push_back(pair.second);
		}
		return suffixes;
	}

private:
	vector<uint64_t> get_c_array(AbstractString &string){
		RankArray r = string.ranks(string.size());
		vector<uint64_t> c_array(r.size());
		for(int i=0; i<c_array.size(); i++){
			for(int j=0; j<i; j++){
				c_array[i] += r[j];
			}
		}
		return c_array;
	}

	vector<uint64_t> get_abundences(vector<uint64_t> C){
		vector<uint64_t> abundances;
		for(int c=0; c<this->characters_cnt-1; c++){
			abundances.push_back(C[c+1]-C[c]);
		}
		abundances.push_back(this->STRING->size()-C[this->characters_cnt-1]);
		return abundances;
	}

	vector<CharId> get_char_ids_by_abundence_desc(vector<uint64_t> abundances){
		vector<CharId> char_id_by_abundance;
		for(int c=0; c<this->characters_cnt; c++){
			char_id_by_abundance.push_back(c);
		}
		// Sort the indices vector based on the values in vector C
		std::sort(char_id_by_abundance.begin(), char_id_by_abundance.end(), [&abundances](CharId a, CharId b) {
			return abundances[a] > abundances[b];
		});
		return char_id_by_abundance;
	}

};

#endif /* FM_INDEX_INDEX_HPP_ */
