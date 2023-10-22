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
string input_prokrustean;
string output_file;
int num_threads=12;
char TERM = '$';

void help(){
	cout << "prokrustean sequence access [options]" << endl <<
	"Input: ebwt of a collection of sequences. Input2: prokrustean graph" << endl <<
	"Output: A data structure representing supporting sequence access for Prokrustean Graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-p <arg>    (REQUIRED) input prokrustean graph" << endl <<
	"-o <arg>    output file. default: {input}.prokrustean.sequences" << endl <<
	"-t <arg>    thread count. default: " << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "i:p:o:h")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_bwt = string(optarg);
			break;
			case 'p':
				input_prokrustean = string(optarg);
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
	if(input_bwt.size()==0){
		cout << "input2 prokrustean empty" << endl;
		help();
	}
	if(output_file.size()==0) {
		output_file=input_bwt+".prokrustean.sequences";
	};

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Input prokrustean file: " << input_prokrustean << endl;
	cout << "Output sequences file: " << output_file << endl;
	cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
    cout << "wavelet tree construction ... " << endl;
	WaveletString str=WaveletString(input_bwt);
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	if(!str.valid){ cout << "wavelet tree is not built. invalid input_bwt path?" << endl; exit(0);}
    
	start = std::chrono::steady_clock::now();
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		cout << "loading failed. " << input_prokrustean << endl;
	}
	prokrustean.print_abstract();
	
	start = std::chrono::steady_clock::now();
	cout << "annotate strata example occurrences (so that they can be printed) ... ";
	
	ProkrusteanExtension ext(prokrustean);
	// setup_stratum_example_occ(ext);
    setup_stratum_example_occ_parallel(ext, num_threads);
	
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	
	FmIndex fm_idx(str);
	if(prokrustean.sequence_count != fm_idx.seq_cnt()){
		cout << "discrepancy between prokrustean and bwt sequences: " << "prokrustean sequence count " << prokrustean.sequence_count << " bwt sequence count " << fm_idx.seq_cnt() << endl;
		exit(0);
	}
	start = std::chrono::steady_clock::now();
	cout << "recover sequences from bwt... ";

    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, num_threads);

	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	start = std::chrono::steady_clock::now();
	cout << "indexing sequences (" << output_file << ") ... " << endl;
	
	DiskSequenceAccess sequence_access(output_file);
	sequence_access.write_open();
	sequence_access.write_metadata(prokrustean);
	// sequence_access.write_strings(strs);
	sequence_access.write_close();
	sequence_access.update_open();
	for(int i=0; i< seq_texts.size(); i++){
		sequence_access.update_single_sequence(i, seq_texts[i]);
	}
	sequence_access.update_close();
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

