#ifndef TEST_CONST_CPP
#define TEST_CONST_CPP
#include <iostream>

using namespace std;

/*relative path from tests/build/*/
string PATH1_SEQ = "../data/1_sequences.txt";
string PATH1_BWT = "../data/1_ebwt.txt";
string PATH2_SEQ = "../data/2_sequences_unsorted.txt";
string PATH2_BWT = "../data/2_ebwt.txt";
string PATH3_SEQ = "../data/3_sequences_unsorted_tied.txt";
string PATH3_BWT = "../data/3_ebwt.txt";
string PATH4_SREAD_PARTITIONED="../../../prokrustean_data/SRR20044276.001001.bwt";
string PATH4_SREAD_PARTITIONED_BI="../../../prokrustean_data/SRR20044276.001001.bidirection.bwt";
string PATH5_CDBG_SAMPLE="../../../prokrustean_data/cdbg_sample.bwt";
// If parameter is not true, test fails
// This check function would be provided by the test framework
#define IS_TRUE(x) { if (!(x)) std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; }

#endif