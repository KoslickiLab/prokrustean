#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <set>
#include "const.cpp"
#include "naive_impl.cpp"	
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/braycurtis.hpp"
#include "../src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;


void test_single_dataset_returns_frequency(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);

    Prokrustean prokrustean;
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    vector<DatasetId> ids_by_sequence(prokrustean.sequence_count, 0);
    int dataset_count=1;
    vector<BrayCurtisOutput> outputs;
    BrayCurtisIntermediate intermediate;
    intermediate.from=5;
    intermediate.to=10;
    intermediate.dataset_count=dataset_count;
    intermediate.dataset_ids_by_sequence=ids_by_sequence;
    intermediate.initialize(ext.prokrustean);

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, intermediate);

    // verification
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        auto frequency = ext.prokrustean.stratums__frequency_cnt[i];
        auto computed_frequency = intermediate.stratum_frequencies_by_dataset[0][i];
        assert(frequency==computed_frequency);
    }
}

void test_two_datasets_return_correct_frequencies(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);

    Prokrustean prokrustean;
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    vector<DatasetId> dataset_ids;
    for(int i=0; i<prokrustean.sequence_count; i++){
        if(i%2==0) dataset_ids.push_back(0);
        else dataset_ids.push_back(1);
    }
    int dataset_count=2;
    vector<BrayCurtisOutput> outputs;
    BrayCurtisIntermediate intermediate;
    intermediate.from=5;
    intermediate.to=10;
    intermediate.dataset_count=dataset_count;
    intermediate.dataset_ids_by_sequence=dataset_ids;
    intermediate.initialize(ext.prokrustean);

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, intermediate);

    // verification
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        auto frequency = ext.prokrustean.stratums__frequency_cnt[i];
        auto computed_frequency = intermediate.stratum_frequencies_by_dataset[0][i]+intermediate.stratum_frequencies_by_dataset[1][i];
        assert(frequency==computed_frequency);
    }
}


void test_single_dataset_returns_frequency_sum_denominator(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    vector<string> seqs; 
    fm_idx.recover_all_texts(seqs);

    Prokrustean prokrustean;
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    int dataset_count=1;
    int from=5;
    int to=10;
    vector<BrayCurtisOutput> outputs;
    compute_incoming_degrees(prokrustean, ext.stratum_incoming_degrees);
    // compute_braycurtis_k_range(5, 150, ext, ids_by_sequence, dataset_count, outputs);;

    BrayCurtisIntermediate intermediate;
    intermediate.from=from;
    intermediate.to=to;
    intermediate.dataset_count=dataset_count;
    intermediate.dataset_ids_by_sequence=vector<DatasetId>(prokrustean.sequence_count, 0);
    intermediate.initialize(ext.prokrustean);
    
    compute_frequencies_by_datasets(ext, intermediate);
    
    Vertex vertex;
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        for(int k=from; k<to+1; k++){
            ext.prokrustean.get_stratum(i, vertex, k);
            _count_k_mers_per_datasets(vertex, k, intermediate);
        }
    }

    int least_frequency=2;
    for(auto output: outputs){
        if(output.k%5==0){
            uint64_t k_mer_sum=get_kmer_frequency_sum_naive(seqs, output.k, least_frequency);
            cout << "k: " << output.k << " denominator: " << output.denominator << " kmer sum: " << k_mer_sum << endl;
            assert(k_mer_sum==output.denominator);
        }
    }
}

void test_two_datasets_return_correct_nominator(){
    int Lmin = 1;
    // frequency=1 usually means sequencing error
    int least_frequency=2;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    vector<string> seqs; 
    fm_idx.recover_all_texts(seqs);

    Prokrustean prokrustean;
    prokrustean.contains_stratum_frequency=true;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    int dataset_count=2;
    vector<DatasetId> dataset_ids;
    for(int i=0; i<prokrustean.sequence_count; i++){
        if(i%2==0) dataset_ids.push_back(0);
        else dataset_ids.push_back(1);
    }
    
    vector<BrayCurtisOutput> outputs;
    compute_incoming_degrees(prokrustean, ext.stratum_incoming_degrees);
    compute_braycurtis_k_range(3, 10, ext, dataset_ids, dataset_count, outputs);

    for(auto output: outputs){
        if(output.k%5==0){
            uint64_t nominator_naive=get_braycurtis_nominator_naive(seqs, dataset_ids, output.k, least_frequency, 0, 1);
            // cout << "k: " << output.k << " denominator: " << output.denominator << " nominator: " << output.nominator << " nominator naive: " << nominator_sum << endl;
            // cout << "kmer: " << get_distinct_kmer_count_naive(seqs, output.k, least_frequency)<< endl;;
            assert(nominator_naive==output.nominator);
        }
    }

    // three datasets
    dataset_count=3;
    for(int i=0; i<prokrustean.sequence_count; i++){
        if(i%3==0) dataset_ids.push_back(0);
        else if(i%3==1) dataset_ids.push_back(1);
        else dataset_ids.push_back(2);
    }
    
    compute_incoming_degrees(prokrustean, ext.stratum_incoming_degrees);
    compute_braycurtis_k_range(3, 50, ext, dataset_ids, dataset_count, outputs);
    
    for(auto output: outputs){
        // to make test concise
        if(output.k%5==0 && output.id1==1 && output.id2==2){
            uint64_t nominator_naive=get_braycurtis_nominator_naive(seqs, dataset_ids, output.k, least_frequency, output.id1, output.id2);
            assert(nominator_naive==output.nominator);
        }
    }
}

void main_application_braycurtis() {
    test_single_dataset_returns_frequency();
    test_two_datasets_return_correct_frequencies();
    test_single_dataset_returns_frequency_sum_denominator();
    test_two_datasets_return_correct_nominator();
}
