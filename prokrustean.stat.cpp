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

	cout << "kmer stat [options]" << endl <<
	"Input: prokrustean graph file." << endl <<
	"Output: a statistics file of the graph." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-l <arg>    lmin range left. default: lmin." << endl <<
	"-r <arg>    lmin range right. default: largest sequence size." << endl <<
	"-o <arg>    output file name. Default: input file name + .stat.txt" << endl <<
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
		output_file=input_prokrustean+".stat.txt";
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
	cout << "enumerating graph .. " << endl;

	uint64_t largest_size=0;
	for(auto &size: prokrustean.sequences__size){
		if(largest_size<size) largest_size=size;
	}
	cout << "largest size: " << largest_size << endl;
	
	vector<uint64_t> vertices(largest_size+1);
	vector<uint64_t> edges(largest_size+1);
	vector<uint64_t> cumul_vertices(largest_size+1);
	vector<uint64_t> cumul_edges(largest_size+1);
	
	Vertex vertex;
	for(int i=0; i<prokrustean.sequence_count; i++){
        prokrustean.get_sequence(i, vertex);
		vertices[vertex.size]+=1;
		for(auto& edge: vertex.s_edges){
			edges[edge.size()]+=1;
		}
    }
	for(int i=0; i<prokrustean.stratum_count; i++){
        prokrustean.get_stratum(i, vertex);
		vertices[vertex.size]+=1;
		for(auto& edge: vertex.s_edges){
			edges[edge.size()]+=1;
		}
    }
	cumul_vertices[largest_size]=vertices[largest_size];
	cumul_edges[largest_size]=edges[largest_size];
	cout << "largest size: " << largest_size << " vertex max " << vertices[largest_size] << " edge max " << edges[largest_size] << endl; 
	for(int i=largest_size-1; i>0; i--){
		cumul_vertices[i]=vertices[i]+cumul_vertices[i+1];
		cumul_edges[i]=edges[i]+cumul_edges[i+1];
	}


	cout << "enumration finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	
	cout << "saving graph sizes... " << endl;

	std::ofstream outputFile(output_file);
	for(int i=from; i<=to; i++){
		outputFile << std::left << std::setw(20) << i;
		outputFile << std::left << std::setw(20) << cumul_vertices[i];
		outputFile << cumul_edges[i] << endl;
	}
		
	outputFile.close();
	cout << "total " << (std::chrono::steady_clock::now()-start_total).count()/1000000 << "ms" << endl;
}
