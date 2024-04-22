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
    
    int dataset_count=1;
    vector<DatasetId> dataset_ids_by_sequence(prokrustean.sequence_count, 0);
    vector<vector<FrequencyCount>> stratum_frequencies_by_dataset(dataset_count, vector<FrequencyCount>(ext.prokrustean.stratum_count, 0));
    vector<BrayCurtisOutput> outputs;

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, dataset_count, dataset_ids_by_sequence, stratum_frequencies_by_dataset);

    // verification
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        auto frequency = ext.prokrustean.stratums__frequency_cnt[i];
        auto computed_frequency = stratum_frequencies_by_dataset[0][i];
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
    vector<DatasetId> dataset_ids_by_sequence(prokrustean.sequence_count, 0);
    vector<vector<FrequencyCount>> stratum_frequencies_by_dataset(dataset_count, vector<FrequencyCount>(ext.prokrustean.stratum_count, 0));

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, dataset_count, dataset_ids_by_sequence, stratum_frequencies_by_dataset);

    // verification
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        auto frequency = ext.prokrustean.stratums__frequency_cnt[i];
        auto computed_frequency = stratum_frequencies_by_dataset[0][i]+stratum_frequencies_by_dataset[1][i];
        assert(frequency==computed_frequency);
    }
}

void test_two_datasets_frequencies_parallel(){
    int Lmin = 1;
    int thread_cnt = 4;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);

    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    
    vector<DatasetId> dataset_ids;
    for(int i=0; i<prokrustean.sequence_count; i++){
        if(i%2==0) dataset_ids.push_back(0);
        else dataset_ids.push_back(1);
    }
    int dataset_count=2;
    vector<BrayCurtisOutput> outputs;
    vector<vector<FrequencyCount>> frequencies(dataset_count, vector<FrequencyCount>(ext.prokrustean.stratum_count, 0));

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, dataset_count, dataset_ids, frequencies);

    ext.set_stratum_locks(thread_cnt);
    vector<vector<FrequencyCount>> frequencies_parallel(dataset_count, vector<FrequencyCount>(ext.prokrustean.stratum_count, 0));
    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets_parallel(ext, thread_cnt, dataset_count, dataset_ids, frequencies_parallel);

    // verification
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        assert(frequencies[0][i]+frequencies[1][i]==frequencies_parallel[0][i]+frequencies_parallel[1][i]);
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
    vector<DatasetId> dataset_ids_by_sequence(prokrustean.sequence_count, 0);
    vector<vector<FrequencyCount>> stratum_frequencies_by_dataset(dataset_count, vector<FrequencyCount>(ext.prokrustean.stratum_count, 0));

    compute_incoming_degrees(ext.prokrustean, ext.stratum_incoming_degrees);
    compute_frequencies_by_datasets(ext, dataset_count, dataset_ids_by_sequence, stratum_frequencies_by_dataset);

    BrayCurtisIntermediate intermediate(from, to, dataset_ids_by_sequence, dataset_count, stratum_frequencies_by_dataset);
    
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

void test_two_datasets_nominator_parallel(){
    int Lmin = 1;
    int thread_cnt = 4;
    // usually means sequencing error
    int least_frequency=2;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);

    Prokrustean prokrustean;
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

    ext.set_stratum_locks(thread_cnt);
    vector<BrayCurtisOutput> outputs_parallel;
    compute_incoming_degrees(prokrustean, ext.stratum_incoming_degrees);
    compute_braycurtis_k_range_parallel(3, 10, ext, thread_cnt, dataset_ids, dataset_count, outputs_parallel);

    for(int i=0; i<outputs.size(); i++){
        assert(outputs[i].nominator==outputs_parallel[i].nominator);
    }
}

void main_application_braycurtis() {
    test_single_dataset_returns_frequency();
    test_two_datasets_return_correct_frequencies();
    test_two_datasets_frequencies_parallel();
    test_single_dataset_returns_frequency_sum_denominator();
    test_two_datasets_return_correct_nominator();
    test_two_datasets_nominator_parallel();
}
