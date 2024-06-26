#ifndef TEST_NAIVE_CPP
#define TEST_NAIVE_CPP

#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_map>

using namespace std;

class NaiveFmIndex{

public:
	NaiveFmIndex(vector<string> &sequences, int s_factor=1){
        // sort(sequences.begin(), sequences.end());
        this->sequences = sequences;
        this->sa = get_gsa(this->sequences);
        this->ssa = get_ssa(sa, s_factor);
        this->ebwt = get_ebwt(this->sa);
    }

    char TERM = '$';
    vector<string> sequences;
    vector<pair<int, string>> sa;
    vector<pair<int, string>> ssa;
    vector<char> ebwt;

    vector<pair<int, string>> get_gsa(vector<string> sequences){
        vector<pair<int, string>> suffixes;
        int global_idx = 0;
        for (int i=0; i< sequences.size(); i++){
            auto seq = sequences[i];
            for (int j=0; j< seq.size(); j++){
                // string suffix = seq.substr(j, seq.size()) + TERM;
                string suffix = seq.substr(j);
                suffixes.push_back(make_pair(global_idx, suffix));
                global_idx++;
            }
            // string s(1, TERM);
            // suffixes.push_back(make_tuple(global_idx, s));
            // global_idx++;
        }

        // affects the ordering rule
        std::sort(suffixes.begin(), suffixes.end(), 
        [](tuple<int, string> const &t1, tuple<int, string> const &t2) {
            return get<1>(t1) != get<1>(t2)? get<1>(t1) < get<1>(t2) : get<0>(t1) < get<0>(t2); 
        });

        return suffixes;
    }

    vector<char> get_ebwt(vector<pair<int, string>> sa){
        string concatenated_string = "";
        for (int i=0; i< sequences.size(); i++){
            // concatenated_string += (sequences[i]+TERM);
            concatenated_string += sequences[i];
        }

        vector<char> ebwt;
        for (int i=0; i< sa.size(); i++){
            int global_idx = sa[i].first;
            if(global_idx == 0) ebwt.push_back(TERM);
            else ebwt.push_back(concatenated_string[global_idx-1]);
        }
        return ebwt;
    }

     vector<pair<int, string>> get_ssa(vector<pair<int, string>> sa, int s_factor){
        vector<pair<int, string>> ssa;
        for(int i=0; i<sa.size(); i++){
            if(sa[i].first%s_factor == 0){
                ssa.push_back(sa[i]);
            }
        }

        return ssa;
    }

    void print_sa(){
        cout << "---- suffix array ----" << endl;
        for (int i=0; i< sa.size(); i++){
            cout << sa[i].first << " : " << sa[i].second << endl;
        }
    }

    void print_ssa(){
        cout << "---- sampled suffix array ----" << endl;
        for (int i=0; i< ssa.size(); i++){
            cout << ssa[i].second << endl;
        }
    }

    void print_ebwt(){
        cout << "---- ebwt ----" << endl;
        for (int i=0; i< ebwt.size(); i++){
            cout << ebwt[i] << endl;
        }
    }

    uint64_t get_suffix(uint64_t i){
        return sa[i].first;
    }
};

vector<string> get_sequences_naive(string path, char TERM='$'){
    ifstream ifs(path);
    string seq;
    vector<string> sequences;
    while(ifs.peek()!=EOF){
        char c;
        ifs.read((char*)&c, sizeof(char));
        seq += c;
        if(c==TERM){
            sequences.push_back(seq);
            seq.clear();
        } 

        // if(c==TERM){
        //     sequences.push_back(seq);
        //     seq.clear();
        // } 
        // else {
        //     seq += c;
        // }
    }
    return sequences;
}

struct NaiveNode {
    std::set<std::string> incoming;
    std::set<std::string> outgoing;
};

class NaiveCompactedDeBruijnGraph {
public:
    std::unordered_map<std::string, NaiveNode> graph;
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
    void construct_plain(const std::vector<std::string>& strings, int k) {
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
                        graph[kmer]=NaiveNode();
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
            }
        }
    }

    void construct_compacted(const std::vector<std::string>& strings, int k) {
        graph.clear();
        construct_plain(strings, k);
        
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
    int maximal_unitig_cnt(vector<string> *output=nullptr){
        // int cnt=0;
        // for (const auto& [node, nodeInfo] : graph) {
        //     if(graph[node].outgoing.size()==1 && graph[node].incoming.size()==1){
        //         auto next= graph[node].outgoing.begin();
        //         if(node==*next){
        //             cout << node << endl;
        //             continue;
        //         }
        //     }
        //     cnt++;
        //     if(output!=nullptr){
        //         (*output).push_back(node);
        //     }
        // }
        // return cnt;
        return graph.size();
    }

    void get_right_navigation(unordered_map<string, vector<string>> &unitig_outgoing){
        for(auto& [s, node]: graph){
            if(node.outgoing.size()>0){
                for(auto &o: node.outgoing){
                    unitig_outgoing[s].push_back(o);
                }
            } else {
                unitig_outgoing[s]=vector<string>();
            }
        }
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

vector<string> get_distinct_kmers_naive(vector<string> sequences, unsigned int k){
    set<string> mers;
    for(auto seq: sequences){
        if(seq.size()<k) continue;

        for(int i=0; i<seq.size()-(k-1); i++){
            string mer = seq.substr(i, k);
            mers.insert(mer);
        }
    }
    vector<string> output(mers.begin(), mers.end());
    return output;
}

uint64_t get_distinct_kmer_count_naive(const std::vector<std::string>& sequences, int k, int least_frequency) {
    std::unordered_map<std::string, uint64_t> kmer_count;

    // Count all k-mers across all sequences
    for (const std::string& seq : sequences) {
        if (seq.size() >= k) {  // Check if the sequence is long enough to contain any k-mers
            for (size_t i = 0; i <= seq.size() - k; i++) {
                // Extract the k-mer starting at position i
                std::string kmer = seq.substr(i, k);
                kmer_count[kmer]++;
            }
        }
    }

    // Sum up all frequencies
    uint64_t total_sum = 0;
    for (const auto& pair : kmer_count) {
        if(least_frequency<=pair.second){
            total_sum += 1;
        }
    }
    return total_sum;
}

vector<int> get_pattern_frequencies_naive(vector<string> sequences, vector<string> patterns){
    vector<int> frequencies;
    for(auto pattern: patterns){
        int count = 0;  // Initialize count of occurrences

        // Iterate over each string in the vector
        for (const std::string& seq : sequences) {
            // Search for the pattern in the string starting from position 0
            std::size_t pos = seq.find(pattern, 0);

            // Continue finding the pattern until the end of the string
            while (pos != std::string::npos) {
                count++;  // Increment count for each occurrence found
                // Move to the position after the last found occurrence to find new occurrences
                pos = seq.find(pattern, pos + 1);
            }
        }

        frequencies.push_back(count);
    }
    return frequencies;
}

uint64_t get_kmer_frequency_sum_naive(const std::vector<std::string>& sequences, int k, int least_frequency) {
    std::unordered_map<std::string, uint64_t> kmer_count;

    // Count all k-mers across all sequences
    for (const std::string& seq : sequences) {
        if (seq.size() >= k) {  // Check if the sequence is long enough to contain any k-mers
            for (size_t i = 0; i <= seq.size() - k; i++) {
                // Extract the k-mer starting at position i
                std::string kmer = seq.substr(i, k);
                kmer_count[kmer]++;
            }
        }
    }

    // Sum up all frequencies
    uint64_t total_sum = 0;
    for (const auto& pair : kmer_count) {
        if(least_frequency<=pair.second){
            total_sum += pair.second;
        }
    }
    return total_sum;
}

uint64_t get_braycurtis_nominator_naive(const std::vector<std::string>& sequences, vector<uint8_t> dataset_ids, int k, int least_frequency, uint8_t dt_id1, uint8_t dt_id2) {
    std::unordered_map<std::string, uint64_t> kmer_count1;
    std::unordered_map<std::string, uint64_t> kmer_count2;
    std::unordered_map<std::string, uint64_t> kmer_count_merged;

    for(int id=0; id<sequences.size(); id++){
        // assert(dataset_ids[id]==dt_id1 || dataset_ids[id]==dt_id2);
        auto &seq=sequences[id];
        if (seq.size() >= k) {  // Check if the sequence is long enough to contain any k-mers
            for (size_t i = 0; i <= seq.size() - k; i++) {
                // Extract the k-mer starting at position i
                std::string kmer = seq.substr(i, k);
                kmer_count_merged[kmer]+=1;
                if(dataset_ids[id]==dt_id1){
                    kmer_count1[kmer]+=1;
                    kmer_count2[kmer]+=0; // for key registeration
                }
                if(dataset_ids[id]==dt_id2){
                    kmer_count1[kmer]+=0; // for key registeration
                    kmer_count2[kmer]+=1;
                }
            }
        }
    }

    // sum up min frequencies
    uint64_t total_sum = 0;
    for (const auto& pair : kmer_count_merged) {
        if(least_frequency<=pair.second){
            total_sum += min(kmer_count1[pair.first], kmer_count2[pair.first]);
        }
    }
    return total_sum;
}


// Function to calculate overlap degree of each sequence in a vector
// std::vector<uint32_t> count_overlap_degrees_naive(const std::vector<std::string> &sequences, size_t threshold) {
//     std::vector<uint32_t> overlapDegrees(sequences.size(), 0);

//     for (size_t i = 0; i < sequences.size(); ++i) {
//         for (size_t j = 0; j < sequences.size(); ++j) {
//             size_t overlap = 0;
//             for (size_t p = 0; p < std::min(sequences[i].size(),sequences[j].size()); ++p) {
//                 if (sequences[i][sequences[i].size() - p - 1] == sequences[j][p]) {
//                     overlap++;
//                     if(overlap>=threshold){
//                         overlapDegrees[i]++;
//                         overlapDegrees[j]++;
//                         break;
//                     }
//                 } else {
//                     break;
//                 }
//             }
//         }
//     }

//     return overlapDegrees;
// }


void count_overlap_degrees_naive(const std::vector<std::string> &sequences, uint32_t threshold, std::vector<uint32_t> &degrees_in, std::vector<uint32_t> &degrees_out, bool allow_multi_edge) {
    degrees_in.resize(sequences.size(), 0);
    degrees_out.resize(sequences.size(), 0);

    for (int i = 0; i < sequences.size(); ++i) {
        for (int j = 0; j < sequences.size(); ++j) {
            uint32_t maxOverlap = std::min(sequences[i].size(), sequences[j].size());
            if(i==j){
                maxOverlap--;
            }
            // threshold <= len < maxOverlap (not the whole)
            for (uint32_t len = threshold; len <= maxOverlap; ++len) {
                // Check if the suffix of sequences[i] of length 'len' matches the prefix of sequences[j] of the same length
                auto suffix = sequences[i].substr(sequences[i].size() - len);
                auto prefix = sequences[j].substr(0, len);
                if (suffix==prefix) {
                    degrees_out[i]++;
                    degrees_in[j]++;
                    if(allow_multi_edge){
                        // multiple edges can be defined
                    } else {
                        break; // Overlaps in a pair mean at most 1 edge  
                    }
                }
            }
        }
    }
}
#endif