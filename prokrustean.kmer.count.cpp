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
#include "src/application/kmers.count.hpp"

using namespace std;
using namespace sdsl;

string input_prokrustean;
string output_file;
int num_threads=8;
int from=-1;
int to=-1;
char TERM = '$';

void help(){

	cout << "kmer counting [options]" << endl <<
	"Input: prokrustean graph file." << endl <<
	"Output: kmer counts of the given range." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-l <arg>    k range left. default: lmin." << endl <<
	"-r <arg>    k range right. default: largest sequence size." << endl <<
	"-o <arg>    output file name. Default: input file name + .kmer.count.txt" << endl <<
	"-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "hi:l:r:o:p:t")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'p':
				input_prokrustean = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'l':
				from = stoi(string(optarg));
			break;
			case 'r':
				to = stoi(string(optarg));
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
		cout << "input empty" << endl;
		help();
	}
	if(output_file.size()==0) {
		output_file=input_prokrustean+".kmer.count.txt";
	};
	auto start = std::chrono::steady_clock::now();
	auto start_total = std::chrono::steady_clock::now();
	// cout << "loading prokrustean ... " << input_prokrustean << endl;
	Prokrustean prokrustean;
	bool success=load_prokrustean(input_prokrustean, prokrustean);
	if(!success){
		cout << "loading failed " << input_prokrustean << endl;
		exit(0);
	}
	// cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	if(from==-1){
		from=prokrustean.lmin;
	} else if(from<prokrustean.lmin){
		cout << "from(l) value should be at least the lmin of prokrustean: " << prokrustean.lmin << endl;
		exit(0);
	}
	if(to==-1){
		for(auto &size: prokrustean.sequences__size){
			if(to<size) to=size;
		}
	} else if(to<from){
		cout << "to(r) value should be at least from(l) value: " << from << endl;
		exit(0);
	}
	
	cout << "k range " << from << ":" << to << endl;
	cout << "threads: " << num_threads << endl;
	
	prokrustean.print_abstract();
	
	start = std::chrono::steady_clock::now();
	cout << "counting kmers .. " << endl;

	vector<uint64_t> output;
	count_distinct_kmers_of_range_parallel(from, to, num_threads, prokrustean, output);
	cout << "counting kmers finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	
	cout << "saving kmer counts... " << endl;

	std::ofstream outputFile(output_file);
	for(int i=0; i<output.size(); i++){
		if(from<=i && i<=to){
//			outputFile << std::left << std::setw(20) << i;
//			outputFile << output[i] << endl;
			outputFile << i << "," << output[i] << endl;
		}
	}
		
	outputFile.close();
	cout << "total " << (std::chrono::steady_clock::now()-start_total).count()/1000000 << "ms" << endl;
}
