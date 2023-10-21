#include "test_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_application.kmer.cpp"
#include "test_application.kmer.count.cpp"
#include "test_application.dbg.count.cpp"
#include "test_application.cdbg.cpp"
#include "test_storage.cpp"
#include "test_construction.support.cpp"
#include "test_text_indexing.cpp"

int main(void){
    // cout << "--- test bitvector ---" << endl;
    // main_bitvector();   

    // cout << "--- test fm_index ---" << endl;
    // main_fm_index();

    // cout << "--- test construction ---" << endl;
    // main_construction();

    // cout << "--- test kmer ---" << endl;
    // main_application_kmer();

    cout << "--- test kmer counting ---" << endl;
    main_application_kmer_count();
    
    // cout << "--- test cdbg ---" << endl;
    // main_application_cdbg();

    // cout << "--- test unitig counting ---" << endl;
    // main_application_unitig_count();

    // cout << "--- test store ---" << endl;
    // main_test_storage();

    // cout << "--- test support ---" << endl;
    // main_prokrustean_support();

    // cout << "--- test text indexing ---" << endl;
    // main_text_indexing();
}