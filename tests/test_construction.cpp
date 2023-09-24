#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	
#include "../src/construction/algorithms.hpp"
#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/application/kmers.hpp"

using namespace std;
using namespace sdsl;

vector<string> collect_distinct_kmers_naive(vector<string> sequences, unsigned int k){
    set<string> uniques;
    for(auto seq: sequences){
        for(int i=0; i<seq.size()-k; i++){
            string mer = seq.substr(i, k);
            uniques.insert(mer);
        }
    }
    vector<string> mers = vector(uniques.begin(), uniques.end());
    return mers;
}

void test_basic_construction(){
    int Lmin = 2;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    assert(pk.stratums.size()>0);
    assert(pk.seqs.size()>0);
    // for(int i=0; i< pk.stratums.size(); i++){
    //     assert(pk.stratums[i].id == i);
    // }
    // for(int i=0; i< pk.seqs.size(); i++){
    //     assert(pk.seqs[i].id == i);
    // }
}

void test_distinct_kmers(){
    int Lmin = 1;
    auto str = WaveletString(PATH1_BWT);
    auto fm_idx = FmIndex(str);
    Prokrustean pk = build_prokrustean(fm_idx, Lmin, true);
    auto sequences = recover_text(fm_idx);
    print_prokrustean(pk);

    for(int k=2; k<7; k++){
        cout << "k: " << k << endl;
        vector<string> mers = collect_distinct_kmers(pk, k);
        vector<string> mers_naive = collect_distinct_kmers_naive(sequences, k);
        sort(mers.begin(), mers.end());
        sort(mers_naive.begin(), mers_naive.end());
        cout << "-- naive --" << endl;
        for(auto m: mers_naive){
            cout << m << endl;
        }
        cout << "-- pk --" << endl;
        for(auto m: mers){
            cout << m << endl;
        }
        vector<string>::iterator mers_itr = mers.begin();
        vector<string>::iterator mers_naive_itr = mers_naive.begin();
        while(mers_itr<mers.end()){
            if(*mers_itr != *mers_naive_itr){
                cout << "not matched: " << *mers_itr << ", " << *mers_naive_itr << endl;
            }
            assert(*mers_itr == *mers_naive_itr);
            mers_itr++;
            mers_naive_itr++;
        }
        assert(mers_naive.size()==mers.size());
    }
    
    // for(auto s: mers_naive){
    //     cout << s << endl;
    // }
    // for(auto mer: collect_distinct_kmer(pk, 2)){
    //     cout << mer << endl;
    // }
}

void test_min_cover_algo(){
    /* I recorded an error case */
    uint64_t s = 151;
    vector<uint64_t> p = {0, 3, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 43, 45, 47, 48, 53, 54, 55, 59, 62, 71, 75, 76, 79, 88, 99, 121, 122, 123, 128};
    vector<uint64_t> sa = {94131, 93683, 17560, 66398, 302443, 300246, 290773, 249229, 68842, 12606, 42678, 177676, 77776, 35722, 156667, 10381, 162138, 28551, 124515, 201465, 82990, 54168, 243317, 26557, 114284, 164341, 169825, 251220, 31493, 133800, 219397, 229670, 272721, 280427, 52794, 235050, 192379, 140763, 293542, 67257, 184591, 13523, 44335, 183885, 90817};
    vector<vector<uint64_t>> r = {{88, 89, 90}, {85, 86}, {72, 73, 74}, {71, 72, 73}, {70, 71, 72}, {69, 70}, {68, 69}, {67, 68, 69}, {66, 67, 68, 70}, {65, 66, 67}, {64, 65, 66}, {63, 64}, {62, 63, 64}, {61, 62, 63, 64}, {60, 61}, {59, 60}, {57, 58, 59, 60, 61, 63}, {56, 57, 58}, {55, 56}, {54, 55}, {52, 53}, {51, 52, 53, 54}, {50, 51}, {49, 50, 51}, {48, 49, 50, 51}, {47, 48, 49}, {45, 46, 47, 48}, {43, 44, 45}, {41, 42}, {40, 41, 42}, {35, 36, 37, 38}, {34, 35}, {33, 34, 35}, {29, 30}, {26, 27, 28}, {17, 18}, {13, 14, 15, 16, 17, 18}, {12, 13, 14}, {10}, {63}, {10}, {10}, {10}, {10}, {11}};
    vector<vector<uint64_t>> f = {{94124, 94131, 94131}, {93676, 93683}, {17552, 17560, 17560}, {66383, 66398, 66398}, {302402, 302443, 302443}, {300204, 300246}, {290730, 290773}, {249186, 249229, 249229}, {68797, 68842, 68842, 68842}, {12558, 12606, 12606}, {42629, 42678, 42678}, {177626, 177676}, {77726, 77776, 77776}, {35669, 35722, 35722, 35722}, {156612, 156667}, {10326, 10381}, {162082, 162138, 162138, 162138, 162138, 162138}, {28492, 28551, 28551}, {124455, 124515}, {201403, 201465}, {82925, 82990}, {54099, 54168, 54168, 54168}, {243247, 243317}, {26487, 26557, 26557}, {114213, 114284, 114284, 114284}, {164268, 164341, 164341}, {169744, 169825, 169825, 169825}, {251138, 251220, 251220}, {31407, 31493}, {133712, 133800, 133800}, {219305, 219397, 219397, 219397}, {229578, 229670}, {272629, 272721, 272721}, {280335, 280427}, {52700, 52794, 52794}, {234953, 235050}, {192167, 192379, 192379, 192379, 192379, 192379}, {140550, 140763, 140763}, {293542}, {67257}, {184591}, {13521}, {44333}, {183883}, {90816}};
    
    vector<PositionAnnotation> a;
    for(int i=0; i< p.size(); i++){
        vector<MaximalRepeatAnnotation> m;
        for(int j=0; j<r[i].size(); j++){
            m.push_back({r[i][j], {}, f[i][j]});
        }
        a.push_back({p[i], sa[i], m, r[i]});
    }
    SequenceAnnotation annot = {0, s, a};
    auto mc = get_min_covers(annot);
}


void main_construction_mc() {
    // test_basic_construction();
    test_distinct_kmers();
    // test_min_cover_algo();
}
