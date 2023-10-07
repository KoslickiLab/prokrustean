#include <iostream>
#include <vector>
#include <set>
#include <tuple>

#ifndef TEST_UTIL_CPP
#define TEST_UTIL_CPP

using namespace std;

/*relative path from tests/build/*/
string PATH1_SEQ = "../data/1_sequences.txt";
string PATH1_BWT = "../data/1_ebwt.txt";
string PATH2_SEQ = "../data/2_sequences_unsorted.txt";
string PATH2_BWT = "../data/2_ebwt.txt";
string PATH3_SEQ = "../data/3_sequences_unsorted_tied.txt";
string PATH3_BWT = "../data/3_ebwt.txt";
string PERFORMANCE_DATA_FOLDER = "../../../prokrustean_data";
string PATH1_PERFORMANCE_SREAD_SEQ = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.txt";
// string PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.bwt";
string PATH_SREAD_001001_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.bwt";
string PATH_SREAD_002_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.002.bwt";
string PATH_SREAD_FULL_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.bwt";
string PATH_SREAD_GUT_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/ERR3450203_1.bwt";
string PATH_GENOME_GCA_019_ROPEBWT2_BWT = PERFORMANCE_DATA_FOLDER+"/GCA_019.bwt";
string PATH_SREAD_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR19392960.bwt";

#endif