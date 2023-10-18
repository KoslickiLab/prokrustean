#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <set>
#include "const.cpp"	
#include "../tests/naive_impl.cpp"
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.support.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.hpp"
#include "../src/application/kmers.hpp"
#include "../src/application/cdbg.count.hpp"

using namespace std;
using namespace sdsl;


void _fully_verify_compacted_de_bruijn_graph(vector<string> &sequences, unsigned int k, unordered_map<string, NaiveNode> &graph){
    set<char> characters;
    // get distinct kmers
    set<string> mers;
    for(auto seq: sequences){
        if(seq.size()<k) continue;

        for(int i=0; i<seq.size()-(k-1); i++){
            string mer = seq.substr(i, k);
            mers.insert(mer);
        }
        
        for(auto c: seq){
            characters.insert(c);
        }
    }
    auto func__get_prev_cnt = [](string &mer, set<string> &mers, int k, set<char> &characters) {
        int cnt=0;
        auto mer_prefix=mer.substr(0, k-1);
        for(auto c: characters){
            if(mers.count(c+mer_prefix)>0) 
            cnt++;
        }
        return cnt;
    };
    auto func__get_next_cnt = [](string &mer, set<string> &mers, int k, set<char> &characters) {
        int cnt=0;
        auto mer_suffix=mer.substr(1, k-1);
        for(auto c: characters){
            if(mers.count(mer_suffix+c)>0) 
            cnt++;
        }
        return cnt;
    };
    
    int kmer_cnt_in_unitigs=0;
    for (const auto& [unitig, node]: graph) {
        assert(unitig.size()>=k);

        int last_idx=unitig.size()-k;
        for(int i=0; i<=last_idx; i++){
            auto mer=unitig.substr(i, k);
            assert(mers.count(mer)>0);
            if(i>0){
                // every kmer inside of unitig except first has a single prev
                assert(func__get_prev_cnt(mer, mers, k, characters)==1);
            } else {
                // first mer has the same number of suffix-prefix rel as the incoming of the node
                auto first_mer=unitig.substr(0, k);        
                assert(func__get_prev_cnt(first_mer, mers, k, characters)==node.incoming.size());
            }
            if(i<last_idx){
                // every kmer inside of unitig except last has a single next
                assert(func__get_next_cnt(mer, mers, k, characters)==1);    
            } else {
                // last mer has the same number of suffix-prefix rel as the outgoing of the node
                auto last_mer=unitig.substr(last_idx, k);        
                assert(func__get_next_cnt(last_mer, mers, k, characters)==node.outgoing.size());
            }
            kmer_cnt_in_unitigs++;
        }
        
        // suffix-prefix navigation validation
        auto unitig_prefix=unitig.substr(0, k-1);
        for(auto &i: node.incoming){
            // cout << unitig << ", " << o << endl;
            assert(unitig_prefix==i.substr(i.size()-(k-1), k-1));
        }
        auto unitig_suffix=unitig.substr(unitig.size()-(k-1), k-1);
        for(auto &o: node.outgoing){
            // cout << unitig << ", " << o << endl;
            assert(unitig_suffix==o.substr(0, k-1));
        }
    }

    // check distinct kmers included in unitigs
    cout << "Naive CDBG validation finished: mers: " << kmer_cnt_in_unitigs << " in unitigs " << graph.size() << ", given mers " << mers.size() << endl;
    assert(kmer_cnt_in_unitigs==mers.size());
    
}



void test_naive_cdbg(){
    int Lmin = 1;
    WaveletString str(PATH_SREAD_001001_GRLBWT_BWT, '$');
    auto fm_idx = FmIndex(str);

    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    // int k = 30;
    for(int k=2; k<100; k++){
        if(k%4!=0){
            continue;
        }
        cout << " ------- " << "k: " << k << " ------ " << endl;
        NaiveCompactedDeBruijnGraph cdbg;
        cdbg.construct_compacted(seq_texts, k);

        _fully_verify_compacted_de_bruijn_graph(seq_texts, k, cdbg.graph);
    }
}


void main_application_naive_cdbg() {
    test_naive_cdbg();
}
