// tbd: license comment
/*
 * prokrustean.braycurtis.cpp
 *  Created on: April, 2024
 *      Author: Adam Park
 */

#include <iostream>
#include <unistd.h>
#include <random>
#include "src/construction/algorithms.hpp"
#include "src/application/overlap.hpp"

using namespace std;
using namespace sdsl;

string input_prokrustean;
string output_file;
int num_threads=8;
int min_length=-1;

void help(){

	cout << "overlap computing [options]" << endl <<
	"Input: prokrustean graph file, read sample information " << endl <<
	"Output: overlap graph degrees per each sequence." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-o <arg>    output file name. Default: input file name + .overlap.txt" << endl <<
	"-l <arg>    overlap length threshold. Default: Kmin of the prokrustean graph" << endl <<
	"-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "hi:o:p:t:l:")) != -1){
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
				min_length = stoi(string(optarg));
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
		output_file=input_prokrustean+".overlap.txt";
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
	if(min_length==-1){
		min_length=prokrustean.kmin;
	} else if(min_length<prokrustean.kmin){
		cout << "l (min_length) value is overriden with the kmin of prokrustean: " << prokrustean.kmin << endl;
		min_length=prokrustean.kmin;
	}
	
	cout << "min overlap length: " << min_length  << endl;
	cout << "threads: " << num_threads << endl;
	
	prokrustean.print_abstract();
	
	start = std::chrono::steady_clock::now();
	cout << "counting overlap degrees .. " << endl;

	vector<uint32_t> degrees_in;
	vector<uint32_t> degrees_out;
	auto ext = ProkrusteanExtension(prokrustean);
    compute_incoming_degrees_parallel(ext, num_threads, ext.stratum_incoming_degrees);
	count_overlap_degrees(ext, min_length, ext.stratum_incoming_degrees, degrees_in, degrees_out, false);
	
	cout << "computing overlap degrees finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	
	cout << "writing outputs... " << endl;
	// int i=0;
	// for(auto value: output){
	// 	cout << "k: " << value.k << ", value: " << value.value << endl;
	// }

	std::ofstream outputFile(output_file);
	for(int i=0; i<prokrustean.sequence_count; i++){
		outputFile << i << "," << degrees_in[i] << "," << degrees_out[i] << endl;
	}
		
	outputFile.close();
	cout << "total " << (std::chrono::steady_clock::now()-start_total).count()/1000000 << "ms" << endl;
}
