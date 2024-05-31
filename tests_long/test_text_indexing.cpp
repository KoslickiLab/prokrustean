#ifndef TEST_TEXT_INDEXING_HPP_
#define TEST_TEXT_INDEXING_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>	
#include "const.cpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/util/string.access.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/cdbg.hpp"

using namespace std;

void test_bwt_prokrustean_indexing_update(){
    WaveletString str(PATH_SREAD_001001_GRLBWT_BWT, '$');
    auto fm_idx = FmIndex(str);
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<tuple<SeqId, int, int>> substring_locs={
        make_tuple(0, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42),
        make_tuple(500, 31, 10)
    };
    Prokrustean prokrustean;
    prokrustean.set_seq_count(seq_texts.size());
    for(int i=0; i<seq_texts.size(); i++){
        prokrustean.sequences__size[i]=seq_texts[i].size();
    }
    prokrustean.stratum_count=1;
    prokrustean.kmin=4;

    vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}
    // construct_prokrustean_parallel(fm_idx, prokrustean, 12);
    auto start = std::chrono::steady_clock::now();
    DiskSequenceAccess sequence_access("test_bwt_prokrustean_indexing.dat");
    sequence_access.write_open();
    sequence_access.write_metadata(seq_sizes, prokrustean);
    sequence_access.write_strings(seq_texts);
    sequence_access.write_close();
    start = std::chrono::steady_clock::now();
    cout << "save meta: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // sequence_access.update_open();
    

    // for(int i=0; i< seq_texts.size(); i++){
    //     sequence_access.update_single_sequence(i, seq_texts[i]);
    // }
    // cout << "save sequences " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // sequence_access.update_close();

    DiskSequenceAccess sequence_access2("test_bwt_prokrustean_indexing.dat");
    sequence_access2.load_metadata();
    // sequence_access2.load_all_strings();
    sequence_access2.metadata.print();
    sequence_access2.read_open();
    string result;
    for(auto [idx, from, size]: substring_locs){
        sequence_access2.read_seq(idx, result);
        assert(seq_texts[idx]==result);
        sequence_access2.read_seq_substr(idx, from, size, result);
        assert(seq_texts[idx].substr(from, size)==result);
    }
    sequence_access2.read_close(); 
    
    DiskSequenceAccess sequence_access3("test_bwt_prokrustean_indexing.dat");
    sequence_access3.load_metadata();

    // load this time!
    sequence_access3.load_all_strings();

    sequence_access3.metadata.print();
    sequence_access3.read_open();
    for(auto [idx, from, size]: substring_locs){
        sequence_access3.read_seq(idx, result);
        // cout << "seq_texts[idx] " << seq_texts[idx] <<endl;
        // cout << "result " << result <<endl; 
        assert(seq_texts[idx]==result);
        sequence_access3.read_seq_substr(idx, from, size, result);
        assert(seq_texts[idx].substr(from, size)==result);
    }
    sequence_access3.read_close(); 
}

// void test_temp_____(){
//     DiskSequenceAccess sequence_access("../../../prokrustean_data/experiments/ERR3450203_1.bwt.prokrustean.sequences");
//     // DiskSequenceAccess sequence_access("../../../prokrustean_data/experiments/SRR20044276.bwt.prokrustean.sequences");
    
// 	sequence_access.load_metadata();
// 	sequence_access.read_open();
//     sequence_access.metadata.print();
//     cout << "seq0: " << sequence_access.sequence_start_positions[27510303] << endl;
//     string seq;
//     sequence_access.read_seq(27510303, seq);
//     cout << "seq0: " << seq << endl;
// }

void main_text_indexing() {
    test_bwt_prokrustean_indexing_update();
    // test_temp();

}

#endif