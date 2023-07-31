#include <iostream>
#include <vector>
#include <set>

#ifndef TEST_UTIL_CPP
#define TEST_UTIL_CPP

using namespace std;
// If parameter is not true, test fails
// This check function would be provided by the test framework
#define IS_TRUE(x) { if (!(x)) std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; }


class NaiveFmIndex{

public:
	NaiveFmIndex(vector<string> sequences){
        sort(sequences.begin(), sequences.end());
        this->sequences = sequences;
        this->sa = get_gsa(this->sequences);
        this->ebwt = get_ebwt(this->sa);
    }

    vector<string> sequences;
    vector<pair<int, string>> sa;
    vector<char> ebwt;

    vector<pair<int, string>> get_gsa(vector<string> sequences){
        vector<pair<int, string>> suffixes;
        int global_idx = 0;
        for (int i=0; i< sequences.size(); i++){
            auto seq = sequences[i];
            for (int j=0; j< seq.size(); j++){
                string suffix = seq.substr(j, seq.size()) + "$";
                suffixes.push_back(make_tuple(global_idx, suffix));
                global_idx++;
            }
            suffixes.push_back(make_tuple(global_idx, "$"));
            global_idx++;
        }

        std::sort(suffixes.begin(), suffixes.end(), 
        [](tuple<int, string> const &t1, tuple<int, string> const &t2) {
            return get<1>(t1) < get<1>(t2); // or use a custom compare function
        });

        return suffixes;
    }

    vector<char> get_ebwt(vector<pair<int, string>> sa){
        string concatenated_string = "";
        for (int i=0; i< sequences.size(); i++){
            concatenated_string += sequences[i];
        }

        vector<char> ebwt;
        for (int i=0; i< sa.size(); i++){
            int global_idx = sa[i].first;
            if(global_idx == 0) ebwt.push_back('$');
            else ebwt.push_back(concatenated_string[global_idx-1]);
        }
        return ebwt;
    }

    void print_sa(){
        for (int i=0; i< sa.size(); i++){
            cout << sa[i].second << endl;
        }
    }
};

#endif