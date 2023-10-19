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

using namespace std;
using namespace sdsl;

string input_bwt;
string input_prokrustean;
string output_file;
int num_threads=12;
int k=-1;
char TERM = '$';

void help(){

	cout << "kmer [options]" << endl <<
	"Input: Prokrustean Graph file." << endl <<
	"Output: distinct kmers count." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-k <arg>    (REQUIRED) k - at least lmin." << endl <<
	"-g <arg>    prokrustean file name. Default: input file name + .prokrustean" << endl <<
	"-o <arg>    output kmer file name. Default: input file name + .kmer.txt" << endl <<
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
				input_bwt = string(optarg);
			break;
			case 'g':
				input_prokrustean = string(optarg);
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

	if(input_bwt.size()==0){
		cout << "input empty" << endl;
		help();
	}
	if(input_prokrustean.size()==0) {
		input_prokrustean=input_bwt+".prokrustean";
	};
	if(output_file.size()==0) {
		output_file=input_bwt+".kmer." + "k"+to_string(k) + ".txt";
	};
	if(k==-1){
		cout << "k not set" << endl;
		exit(0);
	}

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Input prokrustean file: " << input_prokrustean << endl;
	cout << "k: " << k << endl;
	cout << "threads: " << num_threads << endl;

	WaveletString str;

	auto start = std::chrono::steady_clock::now();
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		cout << "loading failed. Prokrustean is built and saved at: " << input_prokrustean << endl;
		start = std::chrono::steady_clock::now();
		cout << "make wavelet tree ... ";
		
		str=WaveletString(input_bwt);
		if(!str.valid){ cout << "wavelet tree is not built. invalid input_bwt path?"; exit(0);}

		cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
		FmIndex fm_idx(str);
		construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, k);
		store_prokrustean(prokrustean, input_prokrustean);
	}
	prokrustean.print_abstract();

	if(k<prokrustean.lmin){
		cout << "k has to be at least lmin. given k: " << k << ", lmin of prokrustean: " << prokrustean.lmin  << endl;
		exit(0);
	}

	start = std::chrono::steady_clock::now();
	cout << "annotate strata example occurrences (so that they can be printed) ... ";
	
	ProkrusteanExtension ext(prokrustean);
	// setup_stratum_example_occ(ext);
    setup_stratum_example_occ_parallel(ext, num_threads);
	
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	if(!str.valid){
		start = std::chrono::steady_clock::now();
		cout << "make wavelet tree ... ";
		
		str=WaveletString(input_bwt);
		if(!str.valid){ cout << "wavelet tree is not built. invalid input_bwt path?"; exit(0);}

		cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	}
	
	FmIndex fm_idx(str);
	
	start = std::chrono::steady_clock::now();
	cout << "recover sequences from bwt... ";

    vector<string> seq_texts;
    recover_sequences_parallel(fm_idx, seq_texts, num_threads);
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	start = std::chrono::steady_clock::now();
	cout << "find distinct kmers... ";

    vector<string> mers;
    get_distinct_kmers_parallel(k, ext, seq_texts, mers);

	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	start = std::chrono::steady_clock::now();
	cout << "kmers total: " << mers.size() << endl;
	cout << "save distinct kmers... ";

	store_kmers(mers, output_file);
	
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

