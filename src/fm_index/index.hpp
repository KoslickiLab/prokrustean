#ifndef FM_INDEX_INDEX_HPP_
#define FM_INDEX_INDEX_HPP_

#include <algorithm>
#include <ranges>
#include <future>
#include "../data_types.hpp"

using namespace std;

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
		this->char_abundances = get_abundances(this->C);
		this->characters_ranked_by_abundance = get_char_ids_by_abundance_desc(this->char_abundances);
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

	string recover_text(int seq_id, bool include_term=false){
		uint64_t L = seq_id;
		uint64_t F = this->LF(L);

		string seq;
		while(F >= this->seq_cnt()){
			seq = (*this->STRING)[L] + seq;
			L = F;
			F = this->LF(L);
		}
		
		if(include_term){
			seq += (*this->STRING)[L]; 
		}
		return seq;
	}

	void recover_all_texts(vector<string> &seqs, bool include_term=false){
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
	vector<pair<uint64_t, string>> recover_suffix_array(int seq_no, bool with_term=false){
		uint64_t L = seq_no;
		uint64_t F = this->LF(L);

		string seq;
		if(with_term){
			seq.push_back(this->TERM);
		}
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
	vector<string> recover_suffix_array(bool with_term=true){
		vector<pair<uint64_t, string>> sa;
		for(int i=0; i<this->seq_cnt(); i++){
			for(auto pair: recover_suffix_array(i, with_term)){
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

	vector<uint64_t> get_abundances(vector<uint64_t> C){
		vector<uint64_t> abundances;
		for(int c=0; c<this->characters_cnt-1; c++){
			abundances.push_back(C[c+1]-C[c]);
		}
		abundances.push_back(this->STRING->size()-C[this->characters_cnt-1]);
		return abundances;
	}

	vector<CharId> get_char_ids_by_abundance_desc(vector<uint64_t> abundances){
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

/********************************************************************************************************/
/*                              parallel text recovery                                                  */
/********************************************************************************************************/
auto func__recover_texts = [](FmIndex &fm_index, vector<string> &output, atomic<int> &idx_gen) {
    while(true){
        auto idx = idx_gen.fetch_add(1);
        if(idx>=fm_index.seq_cnt()){
            break;
        }
        output[idx]=fm_index.recover_text(idx);
    }
};

void recover_sequences_parallel(FmIndex &fm_idx, vector<string> &sequences, int num_threads){
    sequences.resize(fm_idx.seq_cnt());
    vector<future<void>> futures;
    atomic<int> seq_id_iter;
    for(int i=0; i<num_threads; i++){
        futures.push_back(std::async(std::launch::async, func__recover_texts, ref(fm_idx), ref(sequences), ref(seq_id_iter)));
    }
    for (auto &f : futures) {
        f.wait();
    }
}
#endif /* FM_INDEX_INDEX_HPP_ */
