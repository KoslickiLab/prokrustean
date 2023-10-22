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
#include "src/application/cdbg.hpp"
#include "src/util/string.access.hpp"
#include "src/util/data.store.hpp"

using namespace std;
using namespace sdsl;

string input_access;
string input_prokrustean;
string output_file;
int num_threads=12;
int k=-1;
char TERM = '$';

void help(){

	cout << "compacted de bruijn graph construction [options]" << endl <<
	"Input1: -p Prokrustean Graph. Input2: Sequence Access." << endl <<
	"Input2: -s Sequence Access. If it is not acquired constructing prokrustean (with -r), then use prokrustean.access to build the index" << endl <<
	"Output: -o Compacted de bruijn graph in a text file." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) input prokrustean file name" << endl <<
	"-k <arg>    (REQUIRED) k - at least lmin." << endl <<
	"-s <arg>    input sequence access file name. Default: {input prokrustean}.sequences" << endl <<
	"-o <arg>    output file name. Default: {input prokrustean}.cdbg.txt" << endl <<
	"-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "hp:k:s:o:t")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 's':
				input_access = string(optarg);
			break;
			case 'p':
				input_prokrustean = string(optarg);
			break;
			case 'k':
				k = stoi(string(optarg));
			break;
			case 't':
				num_threads = stoi(string(optarg));
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_prokrustean.size()==0){
		cout << "input prokrustean empty" << endl;
		help();
	}
	if(input_access.size()==0) {
		input_access=input_prokrustean+".sequences";
	};
	if(output_file.size()==0) {
		output_file=input_prokrustean+".cdbg." + "k"+to_string(k) + ".txt";
	};
	if(k==-1){
		cout << "k not set" << endl;
		exit(0);
	}

	cout << "Input prokrustean file: " << input_prokrustean << endl;
	cout << "Input access file: " << input_access << endl;
	cout << "k: " << k << endl;
	cout << "threads: " << num_threads << endl;

	WaveletString str;

	auto start = std::chrono::steady_clock::now();
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		cout << "prokrustean loading failed. " << input_prokrustean << endl;
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
	
	DiskSequenceAccess sequence_access(input_access);
	sequence_access.load_metadata();
	sequence_access.read_open();
	// sequence_access.load_all_strings();
	
	start = std::chrono::steady_clock::now();
	cout << "compute unitigs... " << endl;
	
    CompactedDBGWorkspace workspace;
    extract_paritial_unitigs(k, ext, workspace);
	cout << "hihi" << endl;
	update_stratum_based_loc_to_seq_based_loc(ext, workspace);
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	start = std::chrono::steady_clock::now();
	cout << "complete and save cdbg... " << endl;
	DiskStringDataStore string_store(output_file);
	construct_cdbg(workspace.unitigs, sequence_access, string_store, k);
	// store_kmers(mers, output_file);
	cout << "stored: " << output_file << endl;
	cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

