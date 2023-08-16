#include "test_fm_index.cpp"
#include "test_naive_fm_index.cpp"
#include "test_bitvector.cpp"
// #include "test_construction.repr.cpp"
#include "test_construction.maxrep.cpp"
#include "test_construction.navigation.cpp"

int main(void){
    main_fm_index();
    main_naive_fm_index();
    main_bitvector();
    // main_construction_repr();
    main_construction_navigation();
    main_construction_max_rep();
}