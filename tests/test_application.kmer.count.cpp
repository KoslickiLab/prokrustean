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


// if prokrustean if correct, the kmers will be perfectly collected
void test_counting_distinct_kmers(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<uint64_t> counts;
    count_distinct_kmers_of_range(1, 150, prokrustean, counts);

    vector<string> output;
    for(int k=1; k<40; k++){
        get_distinct_kmers(k, ext, seq_texts, output);
        assert(counts[k]==output.size());
    }
}

void main_application_kmer_count() {
    test_counting_distinct_kmers();
}
