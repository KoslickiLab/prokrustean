#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "util.cpp"	
#include "../src/prokrustean.hpp"
// #include "../src/construction/algorithms.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/fm_index/tree_new.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;

struct SampledSuffixArray{
    /* Sample suffix array for each position 
    * Note: seq id: pos pairs should be identified.
    * Then just sampling by regular intervals is not enough because for example first index of a sequence can be missed.
    */
    int sampling_factor;
    FmIndex* fm_index;
    sdsl::bit_vector is_sampled_bv;
    sdsl::rank_support_v<> is_sampled_rank;
    vector<Pos> seq_lengths;
    vector<SeqId> seq_id_samples;
    vector<Pos> reverse_position_samples;

    SampledSuffixArray(FmIndex &fm_index, int sampling_factor){
        this->sampling_factor=sampling_factor;
        this->fm_index=&fm_index;
    }

    void setup_sample_indices(){
        /* may be expensive filling the bit vector*/
        uint64_t seq_cnt=this->fm_index->seq_cnt();
        this->is_sampled_bv=sdsl::bit_vector(this->fm_index->size());
        uint64_t idx=0;
        uint64_t total_length=fm_index->size();
        // mark all idx s.t. idx mod sampling_factor is 0 
        while(idx<total_length){
            this->is_sampled_bv[idx]=true;
            idx+=this->sampling_factor;
        }
        // mark all idx that is the first position of each seq
        // of course, can be duplicated with sampling factor
        for(int i=0; i<seq_cnt; i++){
            int sa_idx = fm_index->STRING->select_by_char(i+1, fm_index->TERM);
            this->is_sampled_bv[sa_idx]=true;

            // remove the termination from samples because non-zero 'position' is not defined
            sa_idx=this->fm_index->LF(sa_idx);
            this->is_sampled_bv[sa_idx]=false;
        }
        
        this->is_sampled_rank=sdsl::rank_support_v<>(&this->is_sampled_bv);
        
        //prepare vector spaces
        uint64_t total_sample_cnt=is_sampled_rank.rank(total_length);

        this->seq_id_samples.resize(total_sample_cnt, -1);
        this->reverse_position_samples.resize(total_sample_cnt);
        this->seq_lengths.resize(seq_cnt);
        
        cout << "suffix samples count: " << total_sample_cnt << endl;
    }

    // for parallelism
    void sample_suffices(int seq_id){
        /* note, the termination($) will not be sampled. so there can be a few void samples. 
        * I leave it just because the logic will not pick the sample.
        */
        uint64_t seq_cnt=this->fm_index->seq_cnt();
        SuffixArrayIdx L = seq_id;
	    SuffixArrayIdx F = this->fm_index->LF(L);
        uint64_t idx;
        Pos reverse_pos=0;
        // until the sequence part is done
        while(F >= seq_cnt){
            L = F;
            F = fm_index->LF(L);
            // sample suffix indices of (mod factor=0) 
            if(L % this->sampling_factor==0){
                assert(this->is_sampled_bv[L]); // TODO: remove later
                idx=this->is_sampled_rank.rank(L);
                this->seq_id_samples[idx]=seq_id;
                this->reverse_position_samples[idx]=reverse_pos;
            }
            reverse_pos++;
        }
        // sample each first position.
        assert(this->is_sampled_bv[L]); // TODO: remove later
        idx=this->is_sampled_rank.rank(L);
        this->seq_id_samples[idx]=seq_id;
        this->reverse_position_samples[idx]=reverse_pos-1;
        this->seq_lengths[seq_id]=reverse_pos;
    }

    tuple<SeqId, Pos> get_location(SuffixArrayIdx suffix_idx){
        /* if input is a termination suffix index, the answer will be invalid (seq length) pos */
        int miss=0;
        while(!this->is_sampled_bv[suffix_idx]){
            suffix_idx=this->fm_index->LF(suffix_idx);
            miss++;
        }
        
        uint64_t sample_idx=is_sampled_rank.rank(suffix_idx);
        SeqId seq_id=this->seq_id_samples[sample_idx];
        Pos sample_pos=this->seq_lengths[seq_id]-this->reverse_position_samples[sample_idx]-1;
        return make_tuple(seq_id, sample_pos+miss);
    }

    void validate(){
        int cnt=0;
        for(auto id: this->seq_id_samples){
            if(id!=-1)
            cnt ++;
            // assert(id!=-1);
        }
        cout << "cout" << cnt << endl;
    }
};


void test_sampling_basic(){
    // auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto index = FmIndex(str);
    auto sample=SampledSuffixArray(index, 3);
    sample.setup_sample_indices();
    for(int i=0; i<index.seq_cnt(); i++){
        sample.sample_suffices(i);
    }
    sample.validate();
    
    auto start = std::chrono::steady_clock::now();
}

void test_sampling_to_fill_all_positions(){
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto index = FmIndex(str);
    sdsl::bit_vector is_sampled_bv(str.size());
    auto sample=SampledSuffixArray(index, 3);
    sample.setup_sample_indices();
    for(int i=0; i<index.seq_cnt(); i++){
        sample.sample_suffices(i);
    }

    vector<vector<int>> suffix_indices_by_pos_and_seq(index.seq_cnt());
    for(int i=0; i<index.seq_cnt(); i++){
        suffix_indices_by_pos_and_seq[i].resize(sample.seq_lengths[i],-1);
    }
    
    for(int i=0; i<index.size(); i++){
        auto pairs=sample.get_location(i);
        SeqId id=get<0>(pairs);
        Pos pos=get<1>(pairs);
        // termination case
        if(pos>=sample.seq_lengths[id]){
            continue;
        }
        assert(suffix_indices_by_pos_and_seq[id][pos]==-1);
        suffix_indices_by_pos_and_seq[id][pos]=i;
    }
    
    int cnt=0;
    for(auto& indices: suffix_indices_by_pos_and_seq){
        for(auto idx: indices){
            assert(idx !=-1);
            cnt++;
        }
    }
    int total_seq_length=0;
    for(auto length: sample.seq_lengths){
        total_seq_length+=length;
    }
    assert(cnt == total_seq_length);
    // auto start = std::chrono::steady_clock::now();
}

void test_sampling_still_works_if_factor_exceeds_seq_length(){
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto index = FmIndex(str);
    sdsl::bit_vector is_sampled_bv(str.size());
    auto sampling_factor=200;
    auto sample=SampledSuffixArray(index, sampling_factor);
    sample.setup_sample_indices();
    for(int i=0; i<index.seq_cnt(); i++){
        sample.sample_suffices(i);
    }
    
    vector<vector<int>> suffix_indices_by_pos_and_seq(index.seq_cnt());
    for(int i=0; i<index.seq_cnt(); i++){
        suffix_indices_by_pos_and_seq[i].resize(sample.seq_lengths[i],-1);
    }
    
    for(int i=0; i<index.size(); i++){
        auto pairs=sample.get_location(i);
        SeqId id=get<0>(pairs);
        Pos pos=get<1>(pairs);
        // termination case
        if(pos>=sample.seq_lengths[id]){
            continue;
        }
        assert(suffix_indices_by_pos_and_seq[id][pos]==-1);
        suffix_indices_by_pos_and_seq[id][pos]=i;
    }
    
    int cnt=0;
    for(auto& indices: suffix_indices_by_pos_and_seq){
        for(auto idx: indices){
            assert(idx !=-1);
            cnt++;
        }
    }
    int total_seq_length=0;
    for(auto length: sample.seq_lengths){
        total_seq_length+=length;
    }
    assert(cnt == total_seq_length);
}

void test_sampling_works_the_same_for_sampling_factors(){
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto index = FmIndex(str);
    sdsl::bit_vector is_sampled_bv(str.size());
    auto sample1=SampledSuffixArray(index, 3);
    auto sample2=SampledSuffixArray(index, 10);
    auto sample3=SampledSuffixArray(index, 1000);
    sample1.setup_sample_indices();
    sample2.setup_sample_indices();
    sample3.setup_sample_indices();
    for(int i=0; i<index.seq_cnt(); i++){
        sample1.sample_suffices(i);
        sample2.sample_suffices(i);
        sample3.sample_suffices(i);
    }
    
    for(int i=0; i<index.size(); i++){
        auto pairs1=sample1.get_location(i);
        auto pairs2=sample2.get_location(i);
        auto pairs3=sample3.get_location(i);
        assert(pairs1==pairs2);
        assert(pairs2==pairs3);
    }
}

void main_performance_step1_new() {
    // test_basic_step3_operation();
    test_sampling_basic();
    test_sampling_to_fill_all_positions();
    test_sampling_still_works_if_factor_exceeds_seq_length();
    test_sampling_works_the_same_for_sampling_factors();
}
