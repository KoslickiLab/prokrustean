#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <set>
#include "util.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.support.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

class CompactedDeBruijnGraph {
private:
    struct Node {
        std::set<std::string> incoming;
        std::set<std::string> outgoing;
    };

    std::unordered_map<std::string, Node> graph;

    bool isLinear(const std::string& node) {
        if(graph[node].outgoing.size()!=1){
            return false;
        } else {
            auto next= graph[node].outgoing.begin();
            if(graph[*next].incoming.size() == 1){
                if(*next!=node){
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
    }

public:
    void generatePlainFromStrings(const std::vector<std::string>& strings, int k) {
        graph.clear();
        // Constructing the initial de Bruijn graph
        for (const auto& str : strings) {
            if (str.size() < k) continue;
            for (size_t i = 0; i <= str.size() - k; ++i) {
                std::string kmer = str.substr(i, k);
                if (i > 0) {
                    std::string prefix = str.substr(i - 1, k);
                    graph[prefix].outgoing.insert(kmer);
                    graph[kmer].incoming.insert(prefix);
                }
            }
        }
    }

    
    void generateFromStrings(const std::vector<std::string>& strings, int k) {
        graph.clear();
        // Constructing the initial de Bruijn graph
        for (const auto& str : strings) {
            if (str.size() < k) continue;
            for (size_t i = 0; i <= str.size() - k; ++i) {
                std::string kmer = str.substr(i, k);
                if (i > 0) {
                    std::string prefix = str.substr(i - 1, k);
                    graph[prefix].outgoing.insert(kmer);
                    graph[kmer].incoming.insert(prefix);
                }
            }
        }

        // Compaction step
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& [node, nodeInfo] : graph) {
                if (isLinear(node)) {
                    std::string next = *nodeInfo.outgoing.begin();

                    // Merge `node` into `prev`
                    std::string merged = node + next.substr(k-1);
                    cout << "single-single found: " << node << " -> " << next << " merged " << merged << endl;
                    
                    // graph[node].outgoing.erase(next);
                    // graph[node].outgoing.insert(merged);

                    // Update connections for adjacent nodes
                    for (const auto& adj : graph[node].incoming) {
                        graph[adj].outgoing.erase(node);
                        graph[adj].outgoing.insert(merged);
                    }
                    for (const auto& adj : graph[next].outgoing) {
                        graph[adj].incoming.erase(next);
                        graph[adj].incoming.insert(merged);
                    }
                    graph[merged].incoming = graph[node].incoming;
                    graph[merged].outgoing = graph[next].outgoing;

                    // Remove old nodes
                    graph.erase(node);
                    graph.erase(next);
                    cout << "graph size: " << graph.size() << endl;

                    changed = true;
                    break;  // Break to avoid iterator invalidation
                }
            }
        }
    }
    int maximal_unitig_cnt(){
        return graph.size();
    }
    void print() {
        for (const auto& [node, nodeInfo] : graph) {
            std::cout << node << " -> ";
            for (const auto& neighbor : nodeInfo.outgoing) {
                std::cout << neighbor << " ";
            }
            std::cout << "\n";
        }
    }
};


int count_maximal_unitigs(ProkrusteanSupport& support, int k){
    int cnt=0;
    Prokrustean& prokrustean = support.prokrustean;
    /* single thread version of filling extensions */
    for(int i=0; i<prokrustean.sequence_count(); i++){
        Sequence seq = prokrustean.get_sequence(i);
        if(seq.s_regions.size()==0){
            cnt++;
        } else if(seq.s_regions[0].from==0 && seq.s_regions[0].size()>=k){
            // 0 is not reflectum
        } else {
            cnt++;
        }
    }
    for(int i=0; i<prokrustean.stratum_count(); i++){
        Stratum stra = prokrustean.get_stratum(i);
        // cout << "left ext count: " << support.stratum_left_ext_count[i]<< endl;
        // if(support.stratum_left_ext_count[i]==1){
        //     continue;
        // }
        if(stra.s_regions.size()==0){
            cnt++;
        } else if(stra.s_regions[0].from==0 && stra.s_regions[0].size()>=k){
            // 0 is not reflectum
        } else {
            cnt++;
        }
    }
    return cnt;
}

void test_unitig_counting(){
    int Lmin = 1;
    WaveletString str(PATH5_CDBG_SAMPLE, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanOptional prokrustean_optional(true);
    ProkrusteanSupport support(prokrustean, prokrustean_optional);
    construct_prokrustean(fm_idx, prokrustean, Lmin, &prokrustean_optional);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    for(auto &txt: seq_texts){
        cout << txt << endl;
    }
    int k = 4;

    CompactedDeBruijnGraph cdbg;
    cdbg.generatePlainFromStrings(seq_texts, k);
    cdbg.print();
    cdbg.generateFromStrings(seq_texts, k);
    cdbg.print();
    cout << "maximal unitigs naive: " << cdbg.maximal_unitig_cnt() << endl;

    // int mu_cnt= count_maximal_unitigs(support, k);
    // cout << "maximal unitigs: " << mu_cnt << endl;
    support.fill_stratum_left_bound_single_right_extension(k);
}


void main_application_unitig_count() {
    test_unitig_counting();
}
