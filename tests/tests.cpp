#include "test_fm_index.cpp"
#include "test_naive_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.repr.annot.cpp"
#include "test_construction.maxrep.cpp"
#include "test_construction.navigation.cpp"
#include "test_construction.mincover.cpp"
#include "test_construction.repr.rule.cpp"

int main(void){
    cout << "--- main_fm_index ---" << endl;
    main_fm_index();
    cout << "--- main_naive_fm_index ---" << endl;
    main_naive_fm_index();
    cout << "--- main_construction_navigation ---" << endl;
    main_construction_navigation();
    cout << "--- main_construction_max_rep ---" << endl;
    main_construction_max_rep();
    cout << "--- main_construction_reprrank ---" << endl;
    main_construction_repr_annot();
    cout << "--- main_bitvector ---" << endl;
    main_bitvector();   
    cout << "--- main_mincover ---" << endl;
    main_construction_mc();
    cout << "--- main_construction_repr ---" << endl;
    main_construction_repr();
}