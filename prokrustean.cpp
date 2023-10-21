// tbd: license comment
/*
 * prokrustean.cpp
 *  Created on: Jul 29, 2023
 *      Author: Adam Park
 */

#include <iostream>
#include <unistd.h>
#include "src/fm_index/index.hpp"
#include "src/fm_index/string.sdsl.hpp"
#include "src/construction/algorithms.hpp"
#include "src/prokrustean.hpp"
#include "src/prokrustean.support.hpp"
#include "src/util/string.access.hpp"

using namespace std;
using namespace sdsl;

int lmin=-1;
string input_bwt;
string output_file;
string output_seq_file;
bool output_txt=false;
int num_threads=12;
bool contains_char_ext = false;
bool contains_frequency = false;
bool recover_seuqneces = false;
char TERM = '$';

void help(){
	cout << "prokrustean [options]" << endl <<
	"Input: ebwt of a collection of sequences." << endl <<
	"Output: A data structure representing Prokrustean Graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-l <arg>    (REQUIRED) lmin - minimum length of a stratum(maximal repeat)" << endl <<
	"-q <arg>    terminator. default:" << TERM << "($) ASCII code. Cannot be the code for A,C,G,T,N." << endl <<
	"-r     	 recover sequences of the bwt and store in {output}.sequences" << endl <<
	"-o <arg>    output file. default: {input}.prokrustean" << endl <<
	"-t <arg>    thread count. default: " << num_threads << endl <<
	"-c          prokrustean is a txt file. Cannot be reused for applications. default: no" << endl <<
	"-e          include stratum character extensions. default: no" << endl <<
	"-f          include stratum frequency. default: no" << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "i:o:l:q:t:cefhr")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_bwt = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'l':
				lmin = stoi(string(optarg));
			break;
			case 't':
				num_threads = stoi(string(optarg));
			break;
			case 'q':
				TERM = atoi(optarg);
			break;
			case 'c':
				output_txt=true;
			break;
			case 'e':
				contains_char_ext = true;
			break;
			case 'f':
				contains_frequency = true;
			case 'r':
				recover_seuqneces = true;
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_bwt.size()==0){
		cout << "input empty" << endl;
		help();
	}
	if(output_file.size()==0) {
		if(!output_txt){
			output_file=input_bwt+".prokrustean";
		} else {
			output_file=input_bwt+".prokrustean.txt";
		}
	};
	output_seq_file = output_file+".sequences";

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output prokrustean file: " << output_file << endl;
	cout << "Lmin: " << lmin << endl;
	cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
    cout << "stage0 wavelet tree ... " ;
    WaveletString str(input_bwt);
    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

	FmIndex fm_idx(str);
	Prokrustean prokrustean;
	prokrustean.contains_stratum_extension_count=contains_char_ext;
	prokrustean.contains_stratum_frequency=contains_frequency;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, lmin);

	prokrustean.print_abstract();
	
	if(!output_txt){
		store_prokrustean(prokrustean, output_file);
	} else {
		store_prokrustean_text(prokrustean, output_file);
	}
	// dispose prokrustean
	// prokrustean=Prokrustean();

	if(recover_seuqneces){
		auto start = std::chrono::steady_clock::now();
		cout << "recovering sequences and indexing them (" << output_seq_file << ") ... " ;
		DiskSequenceAccess sequence_access(output_seq_file);
		sequence_access.write_metadata(prokrustean);
		string str;
		sequence_access.update_mode_open();
		// for(uint64_t i=0; i<fm_idx.seq_cnt(); i++){
		// 	fm_idx.recover_text(i, str);
		// 	sequence_access.write_single_sequence(i, str);
		// }
		vector<string> strs;
		fm_idx.recover_all_texts(strs);
		sequence_access.write_strings(strs);
		sequence_access.update_mode_close();
		cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	}
}

