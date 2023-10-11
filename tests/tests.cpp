#include "test_naive_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_application.kmer.cpp"
#include "test_application.cdbg.count.cpp"
#include "test_application.cdbg.cpp"

int main(void){
    // cout << "--- main_naive_fm_index ---" << endl;
    // main_naive_fm_index();
    // cout << "--- main_bitvector ---" << endl;
    // main_bitvector();   
    // cout << "--- main_construction ---" << endl;
    // main_construction();

    cout << "--- application kmer ---" << endl;
    main_application_kmer();

    // cout << "--- application unitig counting ---" << endl;
    // main_application_unitig_count();
    
    cout << "--- application cdbg ---" << endl;
    main_application_cdbg();
}