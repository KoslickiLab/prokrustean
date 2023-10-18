#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "const.cpp"	
#include "naive_impl.cpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;

void test_step1_binary_storage(){
    auto annot_cnt=8;
    std::vector<StratumId> stratum_ids(annot_cnt, 0);
    stratum_ids[0]=2;
    stratum_ids[1]=6;
    stratum_ids[3]=7;
    stratum_ids[7]=9;
    std::vector<bool> is_primary(annot_cnt, false);
    is_primary[0]=true;
    is_primary[4]=false;
    is_primary[2]=true;

    SuccinctStratifiedData data;
    data.set_data(annot_cnt, stratum_ids, is_primary);
    std::vector<StratumId> stratum_ids2(annot_cnt, 0);
    std::vector<bool> is_primary2(annot_cnt, false);
    data.get_data(annot_cnt,stratum_ids2, is_primary2 );

    for(int i=0; i< annot_cnt; i++){
        assert(stratum_ids[i]==stratum_ids2[i]);
        assert(is_primary[i]==is_primary2[i]);
    }

    annot_cnt=20;
    stratum_ids.clear();
    stratum_ids.resize(annot_cnt, 0);
    stratum_ids[0]=3;
    stratum_ids[1]=6;
    stratum_ids[2]=9;
    stratum_ids[11]=30;
    is_primary.clear();
    is_primary.resize(annot_cnt);
    is_primary[0]=true;
    is_primary[1]=false;
    is_primary[2]=true;
    is_primary[17]=true;

    SuccinctStratifiedData data2;
    data.set_data(annot_cnt, stratum_ids, is_primary);
    std::vector<StratumId> stratum_ids3(annot_cnt, 0);
    std::vector<bool> is_primary3(annot_cnt, false);
    data.get_data(annot_cnt,stratum_ids3, is_primary3);

    for(int i=0; i< annot_cnt; i++){
        assert(stratum_ids[i]==stratum_ids3[i]);
        assert(is_primary[i]==is_primary3[i]);
    }
}


void test_basic_construction(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
}

// if prokrustean if correct, the kmers will be perfectly collected
void test_basic_construction_w_kmers(){
    int Lmin = 1;
    auto str = WaveletString(PATH4_SREAD_PARTITIONED);
    // auto str = WaveletString(PATH3_BWT);
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    construct_prokrustean(fm_idx, prokrustean, Lmin);
    ProkrusteanExtension ext(prokrustean);
    setup_stratum_example_occ(ext);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    
    vector<string> output;
    vector<int> testing_ks={1,5,10,20};
    for(auto k: testing_ks){
        get_distinct_kmers(k, prokrustean, seq_texts, output);
        sort(output.begin(), output.end());
        assert(output==get_distinct_kmers_naive(seq_texts, k));
    }
}

void main_construction() {
    test_step1_binary_storage();
    test_basic_construction();
    test_basic_construction_w_kmers();
}
