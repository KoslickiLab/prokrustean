#include "test_performance.bitvector.cpp"
#include "test_performance.construction.cpp"
#include "test_performance.tree.cpp"

int main(void){
    // cout << "--- bitvector ---" << endl;
    // main_performance_bitvector();
    // cout << "--- construction ---" << endl;
    // main_performance_construction();
    cout << "--- tree ---" << endl;
    main_performance_tree();
}