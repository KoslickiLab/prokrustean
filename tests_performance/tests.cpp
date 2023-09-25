#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_tree.cpp"
#include "test_step1_tree_explore.cpp"
#include "test_step2_repr_blocks.cpp"
#include "test_step3_build_prokrustean.cpp"
#include "test_full_process_analysis.cpp"

int main(void){
    // cout << "--- bitvector ---" << endl;
    // main_performance_bitvector();
    // cout << "--- construction ---" << endl;
    // main_performance_construction();
    cout << "--- tree(step1) ---" << endl;
    main_performance_new_tree();
    cout << "--- block(step2) ---" << endl;
    main_performance_new_repr_block();
    cout << "--- full ---" << endl;
    
}