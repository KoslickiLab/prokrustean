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
int num_threads=8;
char TERM = '$';

void help(){
	cout << "prokrustean sequence access [options]" << endl <<
	"Input: ebwt of a collection of sequences. Input2: prokrustean graph" << endl <<
	"Output: A data structure representing supporting sequence access for Prokrustean Graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-o <arg>    output file. default: {input}.prokrustean.sequences" << endl <<
	"-t <arg>    thread count. default: " << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "i:o:h")) != -1){
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
			case 't':
				num_threads = stoi(string(optarg));
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_bwt.size()==0){
		cout << "input1 bwt empty" << endl;
		help();
	}
	if(output_file.size()==0) {
		output_file=input_bwt+".prokrustean.sequences";
	};

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output sequences file: " << output_file << endl;
	cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
    cout << "wavelet tree construction ... " << endl;
	WaveletString str=WaveletString(input_bwt);
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	if(!str.valid){ cout << "wavelet tree is not built. invalid input_bwt path?" << endl; exit(0);}
    
	start = std::chrono::steady_clock::now();
	
	FmIndex fm_idx(str);
	start = std::chrono::steady_clock::now();
	cout << "recover sequences from bwt... ";

    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, num_threads);

	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	start = std::chrono::steady_clock::now();
	cout << "indexing sequences (" << output_file << ") ... " << endl;
	vector<SequenceSize> seq_sizes;
	for(auto &s: seq_texts){
		seq_sizes.push_back(s.size());
	}
	DiskSequenceAccess sequence_access(output_file);
	sequence_access.write_open();
	sequence_access.write_metadata(seq_sizes);
	sequence_access.write_strings(seq_texts);
	sequence_access.write_close();
	// sequence_access.update_open();
	// for(int i=0; i< seq_texts.size(); i++){
	// 	sequence_access.update_single_sequence(i, seq_texts[i]);
	// }
	// sequence_access.update_close();
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

