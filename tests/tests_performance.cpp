#include "test_performance.bitvector.cpp"
#include "test_performance.construction.cpp"

int main(void){
    cout << "--- bitvector ---" << endl;
    main_performance_bitvector();
    cout << "--- construction ---" << endl;
    main_performance_construction();
}