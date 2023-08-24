#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

void test_real_data_ropebwt2(){
    int threads = 4;
    int Lmin = 30;
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;
    vector<string> recovered_sequences = recover_text(fm_idx);
    // cout << "1st seq: " << recover_text(fm_idx, 2) << endl;

    std::ifstream file(PATH1_PERFORMANCE_SREAD_SEQ);
    string line;
    vector<string> sequences;

    // Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    Prokrustean pk = build_prokrustean_parallel(fm_idx, threads, Lmin, true);
    // vector<string> mers = collect_distinct_kmers(pk, 40);
    // cout << "mers count: " << mers.size() << " ex: " << mers[0] <<endl;
	// while(getline(file, line))
	// {
	// 	// cout<<line<<endl;
    //     sequences.push_back(line);
	// }
    // file.close();
    // int i = 1000;
    // // cout << "seq " << i << " : " << sequences[i] << endl;
    // i = 0;
    // int found_cnt = 0;
    // int not_found_cnt = 0;
    // for(auto r_seq: recovered_sequences){
    //     bool found = false;
    //     auto r_seq_no_term = r_seq.substr(0, r_seq.size()-1);
    //     for(auto seq: sequences){
    //         if(r_seq_no_term==seq)
    //         found = true;
    //     }
    //     if(found) found_cnt++;
    //     else not_found_cnt++;
    //     if(i>=1000){
    //         cout << "seq " << i << ": " << sequences[i] << endl;
    //     }
    //     i++;
    // }

    // file.open()); //open a file to perform read operation using file object
    // if (file.is_open()){   //checking whether the file is open
    //     string tp;
    //     while(getline(file, tp)){ //read data from file object and put it into string.
    //         //  cout << tp << "\n"; //print the data of the string
    //         sequences.push_back(tp);
    //     }
    //     file.close(); //close the file object.
    // }
    // cout << "found: "<< found_cnt << "not found: " << not_found_cnt << endl;
    
}

void test_real_data_gut_ropebwt2(){
    int threads = 60;
    int Lmin = 30;
    auto start = std::chrono::steady_clock::now();
    auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');
    cout << "bwt reading takes " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;
    // vector<string> recovered_sequences = recover_text(fm_idx);
    // cout << "1st seq: " << recover_text(fm_idx, 2) << endl;

    start = std::chrono::steady_clock::now();
    // Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    vector<string> mers = collect_distinct_kmers(pk, 40);
    cout << "mers count: " << mers.size() << " ex: " << mers[0] <<endl;
    // cout << "results count: " << results.size() <<endl;
    cout << "collection completed, " << (std::chrono::steady_clock::now()-start).count()/1000000/1000 << "s" << endl;
}

void test_real_data_ful_ropebwt2(){
    int threads = 16;
    int Lmin = 30;
    auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << endl;

    auto start = std::chrono::steady_clock::now();
    // Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    Prokrustean pk = build_prokrustean_parallel(fm_idx, threads, Lmin, true);
    vector<string> mers = collect_distinct_kmers(pk, 40);
    cout << "mers count: " << mers.size() << " ex: " << mers[0] <<endl;
}

void main_performance_construction() {
    test_real_data_ropebwt2();
    test_real_data_gut_ropebwt2();
    test_real_data_ful_ropebwt2();
}
