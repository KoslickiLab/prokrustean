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
#include "../src/util/data.store.hpp"
#include "../src/util/string.access.hpp"

using namespace std;
using namespace sdsl;


// if prokrustean if correct, the kmers will be perfectly collected
void test_distinct_kmers(){
    int Kmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Kmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    MemorySequenceAccess sequence_access(seq_texts);
    MemoryStringDataStore string_store;

    auto &output= string_store.strings;
    for(int k=2; k<10; k++){
        // get_distinct_kmers(k, ext, seq_texts, output);
        output.clear();
        get_distinct_kmers(k, ext, sequence_access, string_store);
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

void test_distinct_kmers_ondisk(){
    int Kmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean_single_thread(fm_idx, prokrustean, Kmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);

    vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}

    DiskSequenceAccess sequence_access("seq_texts.dat");
    sequence_access.write_open();
    sequence_access.write_metadata(seq_sizes);
    // for(int i=0; i<seq_texts.size();i++){
    //     sequence_access.write_sequence(i, seq_texts[i]);
    // }
    sequence_access.write_strings(seq_texts);
    sequence_access.write_close();
    DiskStringDataStore string_store("kmer.txt");
    sequence_access.read_open();

    int k=3;
    get_distinct_kmers(k, ext, sequence_access, string_store);
}

void main_application_kmer() {
    test_distinct_kmers(); 
    // test_distinct_kmers_ondisk();
}
