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
	"-l <arg>    kmin range left. default: kmin." << endl <<
	"-r <arg>    kmin range right. default: largest sequence size." << endl <<
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
		from=prokrustean.kmin;
	} else if(from<prokrustean.kmin){
		cout << "from(l) value should be at least the kmin of prokrustean: " << prokrustean.kmin << endl;
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
	
	// start = std::chrono::steady_clock::now();
	// cout << "enumerating graph .. " << endl;

	// uint64_t largest_size=0;
	// for(auto &size: prokrustean.sequences__size){
	// 	if(largest_size<size) largest_size=size;
	// }
	// cout << "largest size: " << largest_size << endl;
	
	// vector<uint64_t> vertices(largest_size+1);
	// vector<uint64_t> edges(largest_size+1);
	// vector<uint64_t> cumul_vertices(largest_size+1);
	// vector<uint64_t> cumul_edges(largest_size+1);
	
	Vertex vertex;
	// for(int i=0; i<prokrustean.sequence_count; i++){
    //     prokrustean.get_sequence(i, vertex);
	// 	vertices[vertex.size]+=1;
	// 	for(auto& edge: vertex.s_edges){
	// 		edges[edge.size()]+=1;
	// 	}
    // }
	// for(int i=0; i<prokrustean.stratum_count; i++){
    //     prokrustean.get_stratum(i, vertex);
	// 	vertices[vertex.size]+=1;
	// 	for(auto& edge: vertex.s_edges){
	// 		edges[edge.size()]+=1;
	// 	}
    // }
	// cumul_vertices[largest_size]=vertices[largest_size];
	// cumul_edges[largest_size]=edges[largest_size];
	// cout << "largest size: " << largest_size << " vertex max " << vertices[largest_size] << " edge max " << edges[largest_size] << endl; 
	// for(int i=largest_size-1; i>0; i--){
	// 	cumul_vertices[i]=vertices[i]+cumul_vertices[i+1];
	// 	cumul_edges[i]=edges[i]+cumul_edges[i+1];
	// }


	// cout << "enumration finished " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
	
	if(prokrustean.contains_stratum_extension_count){
		cout << "started calculating trapezoids " << endl;
		uint64_t over_1byte_count, over_2byte_count, zero_strat_count, one_strat_count, two_strat_count, many_stra_count, front_strat_count, end_strat_count;
		uint64_t k_search_two_stra, k_search_max;
		for(int i=0; i<prokrustean.sequence_count; i++){
		    prokrustean.get_sequence(i, vertex);
			if(vertex.s_edges.size()==0){
				zero_strat_count++;
			} else {
				if(vertex.s_edges.size()==1){
					one_strat_count++;
				} else if(vertex.s_edges.size()==2){
					two_strat_count++;
					// find the suffix region of left
					auto left=prokrustean.get_stratum(vertex.s_edges[0].stratum_id);
					uint32_t smallest;
					if(left.s_edges.size()>0 && left.s_edges[left.s_edges.size()-1].to==left.size){
						smallest=left.s_edges[left.s_edges.size()-1].size();
					} else {
						smallest=0;
					}
					auto right=prokrustean.get_stratum(vertex.s_edges[1].stratum_id);
					int search=1;
					if(right.s_edges.size()>0 && right.s_edges[0].from==0){
						auto stra_id = right.s_edges[0].stratum_id;
						while(true){
							right=prokrustean.get_stratum(stra_id);
							if(right.s_edges.size()>0 && right.s_edges[0].from==0 && right.size>smallest){
								search++;
								stra_id = right.s_edges[0].stratum_id;
							} else {
								break;
							}
						}
					}
					
					k_search_two_stra+=search;
					if(k_search_max<search) search;
				} else {
					many_stra_count++;
				}
				if(vertex.s_edges[0].from==0){
					front_strat_count++;
				}
				if(vertex.s_edges[vertex.s_edges.size()-1].to==vertex.size){
					end_strat_count++;
				}
			}
		}
		cout << "seq zero, one, two, many stratifying: " << zero_strat_count << ", " << one_strat_count << ", " << two_strat_count << ", " << many_stra_count << endl;
		cout << "two search avg: " << (float)k_search_two_stra/two_strat_count << endl;
		cout << "seq front stratifying: " << front_strat_count << endl;
		cout << "seq end stratifying: " << end_strat_count << endl;
		zero_strat_count=one_strat_count=two_strat_count=many_stra_count=front_strat_count=end_strat_count=0;
		k_search_two_stra=0;

		for(int i=0; i<prokrustean.stratum_count; i++){
		    prokrustean.get_stratum(i, vertex);
			if(vertex.size>255){
				over_1byte_count++;
				if(vertex.size>65535){
					over_2byte_count++;
				}
			}
			if(vertex.s_edges.size()==0){
				zero_strat_count++;
			} else {
				if(vertex.s_edges.size()==1){
					one_strat_count++;
				} else if(vertex.s_edges.size()==2){
					two_strat_count++;
					// find the suffix region of left
					auto left=prokrustean.get_stratum(vertex.s_edges[0].stratum_id);
					uint32_t smallest;
					if(left.s_edges.size()>0 && left.s_edges[left.s_edges.size()-1].to==left.size){
						smallest=left.s_edges[left.s_edges.size()-1].size();
					} else {
						smallest=0;
					}
					auto right=prokrustean.get_stratum(vertex.s_edges[1].stratum_id);
					int search=1;
					if(right.s_edges.size()>0 && right.s_edges[0].from==0){
						auto stra_id = right.s_edges[0].stratum_id;
						while(true){
							right=prokrustean.get_stratum(stra_id);
							if(right.s_edges.size()>0 && right.s_edges[0].from==0 && right.size>smallest){
								search++;
								stra_id = right.s_edges[0].stratum_id;
							} else {
								break;
							}
						}
					}
					
					k_search_two_stra+=search;
					if(k_search_max<search) search;
				} else {
					many_stra_count++;
				}
				if(vertex.s_edges[0].from==0){
					front_strat_count++;
				}
				if(vertex.s_edges[vertex.s_edges.size()-1].to==vertex.size){
					end_strat_count++;
				}
			}
		}
		cout << "stra zero, one, two, many stratifying: " << zero_strat_count << ", " << one_strat_count << ", " << two_strat_count << ", " << many_stra_count << endl;
		cout << "two search avg: " << (float)k_search_two_stra/two_strat_count << endl;
		cout << "stra front stratifying: " << front_strat_count << endl;
		cout << "stra end stratifying: " << end_strat_count << endl;
		cout << "stra over 1 and 2 bytes length: " << over_1byte_count << ", " << over_2byte_count << endl;
		cout << "finished calculating trapezoids " << endl;
	}

	// cout << "saving graph sizes... " << endl;

	// std::ofstream outputFile(output_file);
	// for(int i=from; i<=to; i++){
	// 	outputFile << std::left << std::setw(20) << i;
	// 	outputFile << std::left << std::setw(20) << cumul_vertices[i];
	// 	outputFile << cumul_edges[i] << endl;
	// }
		
	// outputFile.close();
	// cout << "total " << (std::chrono::steady_clock::now()-start_total).count()/1000000 << "ms" << endl;
}
