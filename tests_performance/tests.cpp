#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_tree.cpp"
#include "test_new_tree.cpp"

int main(void){
    // cout << "--- bitvector ---" << endl;
    // main_performance_bitvector();
    // cout << "--- construction ---" << endl;
    // main_performance_construction();
    cout << "--- tree ---" << endl;
    main_performance_new_tree();
}