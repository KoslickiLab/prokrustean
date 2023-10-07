#ifndef FM_INDEX_SSA_HPP_
#define FM_INDEX_SSA_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "../prokrustean.hpp"
#include "index.hpp"
#include "../sdsl/int_vector.hpp"
#include "../sdsl/rank_support_v.hpp"
#include "../sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

struct SampleBlock{
    sdsl::bit_vector sampled_bv=sdsl::bit_vector(std::numeric_limits<SuffixArrayIdx_InBlock>::max());
    sdsl::rank_support_v<> sampled_rank;
    vector<SeqId> seq_id_samples;
    vector<Pos> reverse_position_samples;

    std::atomic_flag flag = ATOMIC_FLAG_INIT;
    vector<tuple<SeqId, Pos, SuffixArrayIdx_InBlock>> raw_data;

    SampleBlock(){
        this->sampled_bv=sdsl::bit_vector(std::numeric_limits<SuffixArrayIdx_InBlock>::max());
    }

    void add_raw(SuffixArrayIdx_InBlock local_idx, SeqId seq_id, Pos pos){
        while (this->flag.test_and_set(std::memory_order_acquire)) {}
        // lock this critical section
        raw_data.push_back(make_tuple(local_idx, seq_id, pos));
        this->flag.clear(std::memory_order_release);
    }

    void consume_raw(){
        for(auto &r: raw_data){
            this->sampled_bv[get<0>(r)]=true;
        }
        this->sampled_rank=sdsl::rank_support_v<>(&this->sampled_bv);
        int total_cnt=this->sampled_rank.rank(this->sampled_bv.size());
        this->seq_id_samples.resize(total_cnt);
        this->reverse_position_samples.resize(total_cnt);
        for(auto &r: raw_data){
            int idx=this->sampled_rank(get<0>(r));
            this->seq_id_samples[idx]=get<1>(r);
            this->reverse_position_samples[idx]=get<2>(r);
        }
        this->raw_data.clear();
        this->raw_data.shrink_to_fit();
    }
};

class SampledSuffixArray: public AbstractLocator{
    /* Sample suffix array for each position 
    * Note: seq id: pos pairs should be identified.
    * Then just sampling by regular intervals is not enough because for example first index of a sequence can be missed.
    */
    int sampling_factor;
    FmIndex* fm_index;
    vector<SampleBlock> sample_blocks;
    int block_size;
    atomic<int> block_idx_generator;
    atomic<int> sequence_idx_generator;
public:
    int block_count;
    vector<Pos> seq_lengths;
    SampledSuffixArray(){}
    SampledSuffixArray(FmIndex &fm_index, int sampling_factor){
        this->sampling_factor=sampling_factor;
        this->fm_index=&fm_index;
        this->block_size=numeric_limits<SuffixArrayIdx_InBlock>::max();
        this->block_count=fm_index.size()/this->block_size+1;
        this->sample_blocks=vector<SampleBlock>(this->block_count);
        this->seq_lengths.resize(fm_index.seq_cnt());
    }

    bool sample_one_sequence_and_store(){
        /* assuming parallel workers outside
        */
        int seq_id=sequence_idx_generator.fetch_add(1);
        if(seq_id>=fm_index->seq_cnt()){
            return false;
        }
        uint64_t seq_cnt=this->fm_index->seq_cnt();
        SuffixArrayIdx L = seq_id;
	    SuffixArrayIdx F = this->fm_index->LF(L);
        uint64_t idx;
        Pos reverse_pos=0;
        vector<Pos> reverse_positions;
        vector<char> string;
        auto characters = fm_index->STRING->get_characters();
        while(F >= seq_cnt){
            L = F;
            F = fm_index->LF(L);
            // sample suffix indices of (mod factor=0) 
            if(reverse_pos % this->sampling_factor==0){
                // store raw data with local idx
                this->sample_blocks[L/this->block_size].add_raw(L%this->block_size, seq_id, reverse_pos);
            }
            string.push_back(characters[fm_index->STRING->access(L)]);
            reverse_pos++;
        }
        // sample each first position.
        if((reverse_pos-1) % this->sampling_factor!=0){
            this->sample_blocks[L/this->block_size].add_raw(L%this->block_size, seq_id, reverse_pos-1);
        }
        this->seq_lengths[seq_id]=reverse_pos;
        return true;
    }

    bool consume_one_block_and_release(){
        int block_no=block_idx_generator.fetch_add(1);
        if(block_no>=this->block_count){
            return false;
        }
        this->sample_blocks[block_no].consume_raw();
        return true;
    }

    tuple<SeqId, Pos> get_location(SuffixArrayIdx suffix_idx){
        /* if input is a termination suffix index, the answer will be invalid (seq length) pos */
        int original=suffix_idx;
        int miss=0;
        while(!this->sample_blocks[suffix_idx/this->block_size].sampled_bv[suffix_idx%this->block_size]){
            suffix_idx=this->fm_index->LF(suffix_idx);
            miss++;
            if(miss>this->sampling_factor){
                cout << "cycling detected: " << endl;
                cout << "idx: " << original << endl;
                for(int i=0; i<miss; i++){
                    original=this->fm_index->LF(original);
                    cout << "idx: " << original << endl;
                }
                assert(false);
            }
        }
        auto &block=this->sample_blocks[suffix_idx/this->block_size];
        auto sample_idx=block.sampled_rank.rank(suffix_idx%this->block_size);
        
        SeqId seq_id=block.seq_id_samples[sample_idx];
        Pos sample_pos=this->seq_lengths[seq_id]-block.reverse_position_samples[sample_idx]-1;
        
        return make_tuple(seq_id, sample_pos + miss);
    }

    void set_sequence_lengths(Prokrustean &prokrustean){
        auto seq_cnt= this->seq_lengths.size();
        prokrustean.sequences__size.resize(seq_cnt);
        prokrustean.sequences__region.resize(seq_cnt);
        prokrustean.sequences__region_cnt.resize(seq_cnt);
        for(uint64_t i=0; i<seq_cnt; i++){
            prokrustean.sequences__size[i]=this->seq_lengths[i];
        }
    }
    void dispose(){
        this->seq_lengths.clear();
        this->seq_lengths.shrink_to_fit();
    }

    void validate(){
        int total_sequence_length_except_term=0;
        int id;
        for(auto &length: this->seq_lengths){
            total_sequence_length_except_term+=length;
        }
        if(total_sequence_length_except_term!=fm_index->size()-fm_index->seq_cnt()){
            cout << "scanned length: " << total_sequence_length_except_term << ", total length - sequence cnt: " << fm_index->size()-fm_index->seq_cnt() << endl;
            cout << "Possible error in bwt. " << endl;
            cout << "Gathering sampled suffix array, the program tries to scan all possible locations of the sequence. " << endl;
            cout << "However, the string length stored in the fm_index (subtracting the terminal symbols - sequence count) does not match with the scanned symbols" << endl;
            cout << "This has happened with ropebwt2 for example, where running a few fast.gz files result in a wrong result(the original string is not recovered)" << endl;
            assert(total_sequence_length_except_term==fm_index->size()-fm_index->seq_cnt());
        }
        
    }
};

// /*
// debugging purpose
// */
// string recover_text(FmIndex &fm_idx, int seq_no){
// 	uint64_t L = seq_no;
// 	uint64_t F = fm_idx.LF(L);

// 	string seq;
// 	while(F >= fm_idx.seq_cnt()){
// 		seq = fm_idx.get_character(L) + seq;
// 		L = F;
// 		F = fm_idx.LF(L);
// 	}
// 	// must be terminator symbol
// 	seq += fm_idx.get_character(L); 
// 	return seq;
// }

// vector<string> recover_text(FmIndex &fm_idx){
// 	vector<string> seqs;
// 	for(int i=0; i<fm_idx.seq_cnt(); i++){
// 		seqs.push_back(recover_text(fm_idx, i));
// 	}
// 	return seqs;
// }

// /*
// debugging purpose
// index is position in the sequence,
// first: sa index
// second: actual suffix
// */
// vector<pair<uint64_t, string>> recover_suffix_array(FmIndex &fm_idx, int seq_no, bool with_term=true){
// 	uint64_t L = seq_no;
// 	uint64_t F = fm_idx.LF(L);

// 	string seq(1, fm_idx.TERM);
// 	vector<pair<uint64_t, string>> sa;
// 	while(F >= fm_idx.seq_cnt()){
// 		seq = fm_idx.get_character(L) + seq;
// 		sa.push_back(make_pair(F, seq));
// 		L = F;
// 		F = fm_idx.LF(L);
// 	}
// 	// important: this can be misleading because F is randomly (in lexicographical order) chosen
// 	if(with_term){
// 		string term(1, fm_idx.TERM);
// 		sa.insert(sa.begin(), make_pair(F, term));
// 	}
// 	reverse(sa.begin(), sa.end());

// 	return sa;
// }

// /*
// debugging purpose
// index is sa index
// */
// vector<string> recover_suffix_array(FmIndex &fm_idx, bool with_term=true){
// 	vector<pair<uint64_t, string>> sa;
// 	for(int i=0; i<fm_idx.seq_cnt(); i++){
// 		for(auto pair: recover_suffix_array(fm_idx, i, with_term)){
// 			sa.push_back(pair);
// 		}
// 	}
// 	std::sort(sa.begin(), sa.end(), 
//         [](tuple<int, string> const &t1, tuple<int, string> const &t2) {
//             return get<0>(t1) < get<0>(t2); 
//         });
	
// 	vector<string> suffixes;
// 	for(auto pair: sa){
// 		suffixes.push_back(pair.second);
// 	}
// 	return suffixes;
// }

// vector<int> recover_lcp(FmIndex &fm_idx){
// 	vector<int> lcp;
// 	auto sa = recover_suffix_array(fm_idx);
// 	for(int i=1; i<sa.size(); i++){
// 		int pre = 0;
// 		for (int j = 0; j < min(sa[i-1].size(), sa[i].size()); j++) {
// 			if (sa[i-1][j] != sa[i][j] || sa[i][j] == '#')
// 				break;
// 			pre++;
// 		}
// 		lcp.push_back(pre);
// 	}
// 	return lcp;
// }

#endif