#include "test_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_application.kmer.cpp"
#include "test_application.kmer.count.cpp"
#include "test_application.cdbg.count.cpp"
#include "test_application.cdbg.cpp"

int main(void){
    // cout << "--- test bitvector ---" << endl;
    // main_bitvector();   

    // cout << "--- test fm_index ---" << endl;
    // main_fm_index();

    cout << "--- test construction ---" << endl;
    main_construction();

    // cout << "--- test kmer ---" << endl;
    // main_application_kmer();

    // cout << "--- test kmer counting ---" << endl;
    // main_application_kmer_count();
    
    // cout << "--- test cdbg ---" << endl;
    // main_application_cdbg();

    // cout << "--- test unitig counting ---" << endl;
    // main_application_unitig_count();
}