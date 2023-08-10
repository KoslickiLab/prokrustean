#include "test_fm_index.cpp"
#include "test_naive_fm_index.cpp"
#include "test_bitvector.cpp"
#include "test_construction.cpp"
#include "test_construction.navigation.cpp"

int main(void){
    main_fm_index();
    main_naive_fm_index();
    main_bitvector();
    main_construction();
    main_construction_navigation();
}