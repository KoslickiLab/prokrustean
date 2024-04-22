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
#include "../src/application/kmers.hpp"
#include "../src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;


void test_counting_distinct_kmers_single_k_naive(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    for(int k=1; k<180; k++){
        if(k%15==0){
            auto naive_count = get_distinct_kmer_count_naive(seq_texts, k, 1);
            auto computed_count = count_distinct_kmers(k, prokrustean);
            // cout << " naive count: " << naive_count << " count_distinct_kmers(k, prokrustean): " << computed_count << endl;
            assert(naive_count==computed_count);
        }
    }
}

void test_counting_distinct_kmers_single_k(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    MemorySequenceAccess sequence_access(seq_texts);
    MemoryStringDataStore string_store;

    for(int k=1; k<180; k++){
        string_store.strings.clear();
        get_distinct_kmers(k, ext, sequence_access, string_store);
        assert(string_store.strings.size()==count_distinct_kmers(k, prokrustean));
    }
}

// if prokrustean if correct, the kmers will be perfectly collected
void test_counting_distinct_kmers_k_range(){
    int Lmin = 20;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);

    vector<uint64_t> counts;
    count_distinct_kmers_of_range(30, 180, prokrustean, counts);

    for(int k=30; k<180; k++){
        auto single_k_count=count_distinct_kmers(k, prokrustean);
        assert(counts[k]==single_k_count);
    }
}


void test_counting_distinct_kmers_k_range_parallel(){
    int Lmin = 5;
    int thread_cnt=4;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);

    vector<uint64_t> counts;
    count_distinct_kmers_of_range(30, 180, prokrustean, counts);
    vector<uint64_t> counts_parallel;
    count_distinct_kmers_of_range_parallel(30, 180, thread_cnt, prokrustean, counts_parallel);

    for(int k=30; k<180; k++){
        assert(counts[k]==counts_parallel[k]);
    }
}
void main_application_kmer_count() {
    test_counting_distinct_kmers_single_k_naive();
    test_counting_distinct_kmers_single_k();
    test_counting_distinct_kmers_k_range();
    test_counting_distinct_kmers_k_range_parallel();
}
