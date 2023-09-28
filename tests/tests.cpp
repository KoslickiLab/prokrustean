#include "test_naive_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.cpp"

int main(void){
    cout << "--- main_naive_fm_index ---" << endl;
    main_naive_fm_index();
    cout << "--- main_bitvector ---" << endl;
    main_bitvector();   
    cout << "--- main_construction ---" << endl;
    main_construction();
}