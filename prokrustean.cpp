// tbd: license comment
/*
 * prokrustean.cpp
 *
 *  Created on: Jul 29, 2023
 *      Author: Adam Park
 */


#include <iostream>
// #include "internal/lcp.hpp"
#include <unistd.h>
// #include "internal/dna_bwt_n.hpp"

using namespace std;

string input_bwt;
string output_file;

// bool out_da = false;
// uint8_t lcp_size = 1;

// bool containsN = false;

char TERM = '#';

void help(){

	cout << "prokrustean [options]" << endl <<
	"Input: ebwt of a collection of sequences. Output: Prokrustean Index." << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-i <arg>    Input ebwt (REQUIRED)" << endl <<
	"-o <arg>    Output file name (REQUIRED)" << endl <<
	"-t          ASCII code of the terminator. Default:" << int('#') << " (#). Cannot be the code for A,C,G,T,N." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();

	int opt;
	while ((opt = getopt(argc, argv, "hi:o:t:")) != -1){
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
			case 't':
				TERM = atoi(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	if(TERM == 'A' or TERM == 'C' or TERM == 'G' or TERM == 'T' or TERM == 'N'){

		cout << "Error: invalid terminator '" << TERM << "'" << endl;
		help();
	}

	if(input_bwt.size()==0) help();
	if(output_file.size()==0) help();

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output prokrustean file: " << output_file << endl;

	cout << "Constructing wavelet tree ... " << endl;
	cout << "Constructing sampled suffix array ... " << endl;
	cout << "Collecting maximal repeats and representative occurrences ... " << endl;
	cout << "Constructing prokrustean index ... " << endl;
	cout << "Writing output ... " << endl;

	// dna_bwt_t BWT = dna_bwt_t(input_bwt, TERM);

	// cout << "Done. Size of BWT: " << BWT.size()<< endl;

	// lcp<dna_bwt_t, uint64_t> M8(&BWT);
	// cout << "Storing output to file ... " << endl;
	// M8.save_to_file(output_file);

	cout << "Done. " << endl;

}

