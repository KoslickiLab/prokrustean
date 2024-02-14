// tbd: license comment
/*
 * prokrustean.cpp
 *  Created on: Jul 29, 2023
 *      Author: Adam Park
 */

#include <iostream>
#include <unistd.h>
#include <set>
#include "src/fm_index/index.hpp"
#include "src/fm_index/string.sdsl.hpp"
#include "src/construction/algorithms.hpp"
#include "src/prokrustean.hpp"
#include "src/prokrustean.support.hpp"
#include "src/application/overlap.hpp"
#include "src/util/string.access.hpp"
#include "src/util/data.store.hpp"

using namespace std;
using namespace sdsl;

// string input_access;
string input_prokrustean;
string output_file;
// int num_threads=8;
int k=-1;
char TERM = '$';

void help(){

	cout << "hierarchical overlap graph construction [options]" << endl <<
	"Input1: -p Prokrustean Graph." << endl <<
	// "Input2: -s Sequence Access. If it is not acquired constructing prokrustean (with -r), then use prokrustean.access to build the index" << endl <<
	"Output: -o Hierarchical overlap graph in a text file." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) input prokrustean file name" << endl <<
	"-k <arg>    k - at least lmin." << endl <<
	// "-s <arg>    input sequence access file name. Default: {input prokrustean}.sequences" << endl <<
	"-o <arg>    output file name. Default: {input prokrustean}.hog.txt" << endl;
	// "-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "hp:k:s:o:t:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			// case 's':
			// 	input_access = string(optarg);
			// break;
			case 'p':
				input_prokrustean = string(optarg);
			break;
			case 'k':
				k = stoi(string(optarg));
			break;
			// case 't':
			// 	num_threads = stoi(string(optarg));
			// break;
			default:
				help();
			return -1;
		}
	}

	if(input_prokrustean.size()==0){
		cout << "input prokrustean empty" << endl;
		help();
	}
	// if(input_access.size()==0) {
	// 	input_access=input_prokrustean+".sequences";
	// };
	// if(output_file.size()==0) {
	// 	output_file=input_prokrustean+".hog." + "k"+to_string(k) + ".txt";
	// };
	// if(k==-1){
	// 	cout << "k not set" << endl;
	// 	exit(0);
	// }

	cout << "Input prokrustean file: " << input_prokrustean << endl;
	// cout << "Input access file: " << input_access << endl;
	cout << "k: " << k << endl;
	// cout << "threads: " << num_threads << endl;

	auto start = std::chrono::steady_clock::now();
	auto start_total = std::chrono::steady_clock::now();
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		cout << "[failed] prokrustean loading failed. " << input_prokrustean << endl;
		exit(0);
	}
	prokrustean.print_abstract();
	
	if(k<prokrustean.lmin){
		k=prokrustean.lmin;
		cout << "k: " << prokrustean.lmin << ", lmin of prokrustean. " << endl;
	}
	if(output_file.size()==0) {
		output_file=input_prokrustean+".hog." + "k"+to_string(k) + ".txt";
	};
	
	start = std::chrono::steady_clock::now();
	cout << "annotate strata example occurrences (so that they can be printed) ... " << endl;
	
	ProkrusteanExtension ext(prokrustean);
    // setup_stratum_example_occ_parallel(ext, num_threads);
	
	// DiskSequenceAccess sequence_access(input_access);
	// sequence_access.load_metadata();
	// sequence_access.load_all_strings();
	
	start = std::chrono::steady_clock::now();
	cout << "compute and save hog... " << endl;
	DiskStringDataStore string_store(output_file);
	extract_hierarchical_overlap_graph(ext, k, string_store);
    
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	cout << "stored: " << output_file << endl;
	
	// cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	// cout << "total " << (std::chrono::steady_clock::now()-start_total).count()/1000000 << "ms" << endl;
}

