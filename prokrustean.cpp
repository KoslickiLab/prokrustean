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
bool output_txt=false;
int num_threads=12;
char TERM = '$';

void help(){

	cout << "prokrustean [options]" << endl <<
	"Input: ebwt of a collection of sequences." << endl <<
	"Output: A data structure representing Prokrustean Graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-i <arg>    (REQUIRED) input ebwt file name" << endl <<
	"-l <arg>    (REQUIRED) lmin - minimum length of a stratum(maximal repeat)" << endl <<
	"-o <arg>    output file. default: {input}.prokrustean" << endl <<
	"-t <arg>    thread count default: " << num_threads << endl <<
	"-q <arg>    ASCII code of the terminator. default:" << TERM << "($) Cannot be the code for A,C,G,T,N." << endl <<
	"-c          output is a readable txt file. Cannot be reused for applications. default: none" << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "h:i:o:l:q:t:c")) != -1){
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

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output prokrustean file: " << output_file << endl;
	cout << "Lmin: " << lmin << endl;
	cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
    cout << "stage0 wavelet tree ... ";
    WaveletString str(input_bwt);
    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

	FmIndex fm_idx(str);
	Prokrustean prokrustean;
    construct_prokrustean_parallel(fm_idx, prokrustean, num_threads, lmin);

	prokrustean.print_abstract();
	
	if(!output_txt){
		store_prokrustean(prokrustean, output_file);
	} else {
		store_prokrustean_text(prokrustean, output_file);
	}
}

