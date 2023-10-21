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
        
        cout << "counts[k]: " << counts[k] << " count by single k "<< single_k_count << endl;
        
        assert(counts[k]==single_k_count);
    }
}

void main_application_kmer_count() {
    // test_counting_distinct_kmers_single_k();
    test_counting_distinct_kmers_k_range();
}
