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
#include "src/application/kmers.hpp"
#include "src/util/string.access.hpp"
#include "src/util/data.store.hpp"

using namespace std;
using namespace sdsl;

string input_sequences;
string input_prokrustean;
string output_file;
int num_threads=12;
int k=-1;
char TERM = '$';

void help(){

	cout << "kmer [options]" << endl <<
	"Input: ebwt of a collection of sequences. Input2: Prokrustean Graph." << endl <<
	"Output: distinct kmers in a text file." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input prokrustean file name" << endl <<
	"-k <arg>    (REQUIRED) k - at least lmin." << endl <<
	"-s <arg>    sequence file name. (generated with prokrustean by -r) Default: {input}.sequences" << endl <<
	"-o <arg>    output file name. Default: {input}.kmer.txt" << endl <<
	"-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "h:i:k:o:g:t:k")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_prokrustean = string(optarg);
			break;
			case 's':
				input_sequences = string(optarg);
			break;
			case 'k':
				k = stoi(string(optarg));
			break;
			case 'p':
				num_threads = stoi(string(optarg));
			break;
			case 't':
				TERM = atoi(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_prokrustean.size()==0) {
		cout << "input file name not set" << endl;
		exit(0);
	};
	if(input_sequences.size()==0) {
		input_sequences=input_prokrustean+".sequences";
	};
	if(output_file.size()==0) {
		output_file=input_prokrustean+".kmer." + "k"+to_string(k) + ".txt";
	};
	if(k==-1){
		cout << "k not set" << endl;
		exit(0);
	}

	cout << "Input prokrustean file: " << input_prokrustean << endl;
	cout << "Input sequences file: " << input_sequences << endl;
	cout << "k: " << k << endl;
	cout << "threads: " << num_threads << endl;

	WaveletString str;

	auto start = std::chrono::steady_clock::now();
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		exit(0);
	}
	prokrustean.print_abstract();

	if(k<prokrustean.lmin){
		cout << "k has to be at least lmin. given k: " << k << ", lmin of prokrustean: " << prokrustean.lmin  << endl;
		exit(0);
	}

	start = std::chrono::steady_clock::now();
	cout << "annotate strata example occurrences (so that they can be printed) ... " << endl;
	
	ProkrusteanExtension ext(prokrustean);
	// setup_stratum_example_occ(ext);
    setup_stratum_example_occ_parallel(ext, num_threads);
	
	DiskSequenceAccess sequence_access(input_sequences);
	sequence_access.read_open();
	// sequence_access.load_all_strings();
	
	start = std::chrono::steady_clock::now();
	cout << "find and store distinct kmers... " << endl;
	vector<string> output;
	DiskStringDataStore string_store(output_file);
	get_distinct_kmers(k, ext, sequence_access, string_store);
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	cout << "stored: " << output_file << endl;
}

