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
        set<char> characters;
        graph.clear();
        // Constructing the initial de Bruijn graph
        for (const auto& str : strings) {
            for(char c: str){
                characters.insert(c);
            }
            if (str.size() < k) continue;
            for (size_t i = 0; i <= str.size() - k; ++i) {
                std::string kmer = str.substr(i, k);
                if (i > 0) {
                    std::string prefix = str.substr(i - 1, k);
                    graph[prefix].outgoing.insert(kmer);
                    graph[kmer].incoming.insert(prefix);
                } else {
                    if(graph.count(kmer)==0){
                        graph[kmer]=Node();
                    }
                }
                auto kmer_prefix=kmer.substr(0, k-1);
                auto kmer_suffix=kmer.substr(kmer.size()-(k-1), k-1);
                for(auto c: characters){
                    auto node = c+kmer_prefix;
                    if(graph.count(node)>0){
                        graph[kmer].incoming.insert(node);
                        graph[node].outgoing.insert(kmer);
                    }
                    node = kmer_suffix + c;
                    if(graph.count(node)>0){
                        graph[kmer].outgoing.insert(node);
                        graph[node].incoming.insert(kmer);
                    }
                }
                // for (const auto& [node, nodeInfo] : graph){
                //     // node is predecessor
                //     if(node.substr(node.size()-(k-1), k-1)==kmer.substr(0, k-1)){
                //         graph[kmer].incoming.insert(node);
                //         graph[node].outgoing.insert(kmer);
                //     }
                //     // kmer is predecessor
                //     if(node.substr(0, k-1)==kmer.substr(kmer.size()-(k-1), k-1)){
                //         graph[kmer].outgoing.insert(node);
                //         graph[node].incoming.insert(kmer);
                //     }
                // }
            }
        }
    }

    
    void generateFromStrings(const std::vector<std::string>& strings, int k) {
        graph.clear();
        generatePlainFromStrings(strings, k);
        
        // Compaction step
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& [node, nodeInfo] : graph) {
                if (isLinear(node)) {
                    std::string next = *nodeInfo.outgoing.begin();

                    // Merge `node` into `prev`
                    std::string merged = node + next.substr(k-1);
                    // cout << "single-single found: " << node << " -> " << next << " merged " << merged << endl;
                    
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
                    // cout << "graph size: " << graph.size() << endl;

                    changed = true;
                    break;  // Break to avoid iterator invalidation
                }
            }
        }
    }
    int maximal_unitig_cnt(){
        int cnt=0;
        for (const auto& [node, nodeInfo] : graph) {
            if(graph[node].outgoing.size()==1){
                auto next= graph[node].outgoing.begin();
                if(node==*next){
                    continue;
                }
            }
            cnt++;
        }
        return cnt;
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


void test_unitig_counting(){
    int Lmin = 1;
    WaveletString str(PATH4_SREAD_PARTITIONED, '$');
    auto fm_idx = FmIndex(str);
    
    Prokrustean prokrustean;
    ProkrusteanOptional prokrustean_optional(true);
    ProkrusteanSupport support(prokrustean, prokrustean_optional);
    construct_prokrustean(fm_idx, prokrustean, Lmin, &prokrustean_optional);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    // for(auto &txt: seq_texts){
    //     cout << txt << endl;
    // }
    int k = 15;
    // cout << "maximal unitigs " << support.compute_maximal_unitigs(k) << endl;
    // return;

    CompactedDeBruijnGraph cdbg;
    // cdbg.generatePlainFromStrings(seq_texts, k);
    // cdbg.print();
    cdbg.generateFromStrings(seq_texts, k);
    // cdbg.print();

    auto cnt=support.compute_maximal_unitigs(k);
    if(cdbg.maximal_unitig_cnt()==cnt){
        cout << "maximal unitigs matched: " << cnt << endl;
    } else {
        cout << "maximal unitigs unmatched: " << cnt << ", " << cdbg.maximal_unitig_cnt() << endl;
        assert(false);
    }
}


void main_application_unitig_count() {
    test_unitig_counting();
}
