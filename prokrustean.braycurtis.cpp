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
#include "src/application/braycurtis.hpp"

using namespace std;
using namespace sdsl;

string input_prokrustean;
string input_ids;
string output_file;
int num_threads=8;
int from=-1;
int to=-1;
char TERM = '$';

void help(){

	cout << "bray-curtis computing [options]" << endl <<
	"Input: prokrustean graph file, read sample information " << endl <<
	"Output: bray-curtis of the given range of k for all pairs." << endl <<
	"Options:" << endl <<
	"-h          help" << endl <<
	"-p <arg>    (REQUIRED) prokrustean file name" << endl <<
	"-s <arg>    (REQUIRED) A sample ids txt file containing a list of sample ids(0~255) meaning sample ids by sequences. For rows equivalent to each sequence, sample id 0~255 should be written" << endl <<
	"-l <arg>    k range left. default: kmin." << endl <<
	"-r <arg>    k range right. default: largest sequence size." << endl <<
	"-o <arg>    output file name. Default: input file name + .braycurtis.txt" << endl <<
	"-t <arg>    thread count Default:" << num_threads << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2) help();
	int opt;
	while ((opt = getopt(argc, argv, "hi:l:r:o:p:s:t:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'p':
				input_prokrustean = string(optarg);
			break;
			case 's':
				input_ids = string(optarg);
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
		cout << "input prokrustean empty" << endl;
		help();
	}
	if(input_ids.size()==0){
		cout << "input sample ids empty" << endl;
		help();
	}
	if(output_file.size()==0) {
		output_file=input_prokrustean+".braycurtis.txt";
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
		from=prokrustean.kmin;
	} else if(from<prokrustean.kmin){
		cout << "from(l) value should be at least the kmin of prokrustean: " << prokrustean.kmin << endl;
		exit(0);
	}
	if(to==-1){
		// *** the sequence size can be too large.
			// to=prokrustean.sequences__size[0];
			// for(auto &size: prokrustean.sequences__size){
			// 	if(to<size) to=size;
			// }
		to=100;
	} else if(to<from){
		cout << "to(r) value should be at least from(l) value: " << from << endl;
		exit(0);
	}
	
	cout << "k range " << from << ":" << to << endl;
	cout << "threads: " << num_threads << endl;
	
	prokrustean.print_abstract();
	
	cout << "loading sample ids .. " << endl;
	vector<BrayCurtisOutput> output;
    std::ifstream file(input_ids);
    std::vector<uint8_t> ids_by_sequence;
    int number; // Use an int to safely read values that may be out of uint8_t range
    while (file >> number) {
        if (number < 0 || number > 255) {
            std::cerr << "Error: the sequence ids should be in range for uint8_t (0~255): " << number << std::endl;
            exit(0);
        }
        ids_by_sequence.push_back(static_cast<uint8_t>(number));
    }
    file.close();
	if(ids_by_sequence.size()!=prokrustean.sequence_count){
		std::cerr << "Error: the number of sequence ids is different from the number of sequences in the prokrustean graph: " << ids_by_sequence.size() << " vs " << prokrustean.sequence_count << std::endl;
		exit(0);
	}
	// this is too naive because the maximum dataset index governs the space!
	uint8_t max_id = *max_element(ids_by_sequence.begin(), ids_by_sequence.end());
	uint8_t dataset_count = max_id+1;
	// randomly assign
	// std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<> distrib(0, 1); 
	// for(SeqId i=0; i<prokrustean.sequence_count; i++){
	// 	// ids_by_sequence.push_back(distrib(gen));
	// 	if(i<prokrustean.sequence_count/2){
	// 		ids_by_sequence.push_back(0);
	// 	} else {
	// 		ids_by_sequence.push_back(1);
	// 	}
	// }
	start = std::chrono::steady_clock::now();
	cout << "computing bray-curtis with " << (int)dataset_count << " samples " << endl;

	auto ext = ProkrusteanExtension(prokrustean);

    compute_incoming_degrees_parallel(ext, num_threads, ext.stratum_incoming_degrees);
	cout << "computed incoming degrees of the prokrustean graph " << endl;
	compute_braycurtis_k_range_parallel(from, to, ext, num_threads, ids_by_sequence, dataset_count, output);
	cout << "computing bray-curtis finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	assert(output.size()>0);
	cout << "saving results... " << endl;
	cout << "columns: id1, id2, k, value" << endl;
	std::ofstream f(output_file);
	if (!f.is_open()) {
		std::cerr << "Failed to open the file: " << output_file << std::endl;
		exit(0); // Or handle the error appropriately
	}
	for(int i=0; i<output.size(); i++){
		if(from<=output[i].k && output[i].k<=to){
			f << (int)output[i].id1 << "," << (int)output[i].id2 << "," << output[i].k << "," << output[i].value << endl;
		}
		
		// cout << output[i].id1 << "," << (int)output[i].id2 << "," << (int)output[i].k << "," << output[i].value << endl;
	}
	f.close();
}
