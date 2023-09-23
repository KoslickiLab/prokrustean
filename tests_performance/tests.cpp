#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_tree.cpp"
#include "test_new_tree.cpp"
#include "test_new_repr_blocks.cpp"

int main(void){
    // cout << "--- bitvector ---" << endl;
    // main_performance_bitvector();
    // cout << "--- construction ---" << endl;
    // main_performance_construction();
    // cout << "--- tree(step1) ---" << endl;
    // main_performance_new_tree();
    cout << "--- block(step2) ---" << endl;
    main_performance_new_repr_block();
}