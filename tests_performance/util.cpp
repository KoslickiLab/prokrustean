#include <iostream>
#include <vector>
#include <set>
#include <tuple>

#ifndef TEST_UTIL_CPP
#define TEST_UTIL_CPP

using namespace std;

/*relative path from tests/build/*/
string PATH1_SEQ = "../data/1_sequences.txt";
string PATH1_BWT = "../data/1_ebwt.txt";
string PATH2_SEQ = "../data/2_sequences_unsorted.txt";
string PATH2_BWT = "../data/2_ebwt.txt";
string PATH3_SEQ = "../data/3_sequences_unsorted_tied.txt";
string PATH3_BWT = "../data/3_ebwt.txt";
string PERFORMANCE_DATA_FOLDER = "../../../prokrustean_data";
string PATH1_PERFORMANCE_SREAD_SEQ = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.txt";
// string PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.bwt";
string PATH_SREAD_001001_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.001001.bwt";
string PATH_SREAD_002_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.002.bwt";
string PATH_SREAD_FULL_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR20044276.bwt";
string PATH_SREAD_GUT_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/ERR3450203_1.bwt";
string PATH_GENOME_GCA_019_ROPEBWT2_BWT = PERFORMANCE_DATA_FOLDER+"/GCA_019.bwt";
string PATH_SREAD_GRLBWT_BWT = PERFORMANCE_DATA_FOLDER+"/SRR19392960.bwt";

// If parameter is not true, test fails
// This check function would be provided by the test framework
#define IS_TRUE(x) { if (!(x)) std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; }

class NaiveFmIndex{

public:
	NaiveFmIndex(vector<string> sequences, int s_factor=1){
        // sort(sequences.begin(), sequences.end());
        this->sequences = sequences;
        this->sa = get_gsa(this->sequences);
        this->ssa = get_ssa(sa, s_factor);
        this->ebwt = get_ebwt(this->sa);
    }

    char TERM = '#';
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
                string suffix = seq.substr(j, seq.size());
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


vector<string> get_sequences(string path){
    ifstream ifs(path);
    string seq;
    vector<string> sequences;
    while(ifs.peek()!=EOF){
        char c;
        ifs.read((char*)&c, sizeof(char));
        seq += c;
        if(c=='#'){
            sequences.push_back(seq);
            seq.clear();
        } 
        // else {
        //     seq += c;
        // }
    }
    return sequences;
}

#endif