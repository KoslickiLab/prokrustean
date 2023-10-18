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
void test_distinct_kmers(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<string> output;
    for(int k=2; k<10; k++){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        sort(output.begin(), output.end());
        auto output_naive = get_distinct_kmers_naive(seq_texts, k);
        if(output!=output_naive){
            bool has_duplicates = false;
            for (size_t i = 1; i < output.size(); ++i) {
                if (output[i] == output[i-1]) {  // Compare current element with previous element
                    // Print the duplicate element, only if it hasn't been printed before
                    if(i == 1 || output[i] != output[i-2]) {
                        std::cout << "Duplicate element: " << output[i] << std::endl;
                        has_duplicates = true;
                    }
                }
            }
        }
        // cout << "output " << output.size() << " output_naive " << output_naive.size() << endl;
        assert(output==output_naive);
    }
}

void main_application_kmer() {
    test_distinct_kmers(); 
}
