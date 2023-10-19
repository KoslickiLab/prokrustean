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

using namespace std;
using namespace sdsl;

int lmin=-1;
string input_bwt;
string input_pkg;
string output_kmer;
int num_threads=12;
char TERM = '$';

void help(){

	cout << "kmer [options]" << endl <<
	"Input: ebwt of a collection of sequences. Input2: Prokrustean Graph." << endl <<
	"Output: distinct kmers." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-k <arg>    (REQUIRED) k - has to be at least lmin used in prokrustean construction" << endl <<
	"-o <arg>    output kmer file name. Default: input file name + .kmer.txt" << endl <<
	"-g <arg>    prokrustean file name. Default: input file name + .pkg" << endl <<
	"-p <arg>    thread no. Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "h:i:k:o:g:p")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_bwt = string(optarg);
			break;
			case 'g':
				input_pkg = string(optarg);
			break;
			case 'l':
				lmin = stoi(string(optarg));
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
	if(input_pkg.size()==0) {
		input_pkg=input_bwt+".prokrustean";
	};
	if(output_kmer.size()==0) {
		output_kmer=input_bwt+".kmer.txt";
	};

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Input prokrustean file: " << input_pkg << endl;
	cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
    cout << "constructing sdsl wavelet tree ... " << endl;
    auto str = WaveletString(input_bwt);
    auto fm_idx = FmIndex(str);
    cout << "finished wavelet tree construction " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
	Prokrustean prokrustean;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, lmin);

	start = std::chrono::steady_clock::now();
	cout << "storing prokrustean ... " ;
    store_prokrustean(prokrustean, output_kmer);
	cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

}

