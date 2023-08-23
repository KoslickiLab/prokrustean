#include "test_performance.bitvector.cpp"
#include "test_performance.construction.cpp"

int main(void){
    cout << "--- bitvector ---" << endl;
    main_bitvector();
    cout << "--- construction ---" << endl;
    main_construction_mc();
}