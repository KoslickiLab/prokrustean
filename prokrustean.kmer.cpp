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
string output_file;
int num_threads=12;
char TERM = '$';

void help(){

	cout << "prokrustean kmer [options]" << endl <<
	"Input: ebwt of a collection of sequences. Output: A data structure representing Prokrustean Graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-g <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-o <arg>    (REQUIRED) output file" << endl <<
	"-l <arg>    (REQUIRED) lmin - minimum length of a stratum(maximal repeat)" << endl <<
	"-p <arg>    thread no. Default:" << num_threads << endl <<
	"-t <arg>    ASCII code of the terminator. Default:" << TERM << " ($). Cannot be the code for A,C,G,T,N." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "h:i:o:l:p:t:")) != -1){
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
	if(output_file.size()==0) {
		cout << "output empty" << endl;
		help();
	};

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output prokrustean file: " << output_file << endl;
	cout << "Lmin: " << lmin << endl;
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
    store_prokrustean(prokrustean, output_file);
	cout << "finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

}

