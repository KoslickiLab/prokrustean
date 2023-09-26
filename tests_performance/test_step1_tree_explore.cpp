#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include "util.cpp"	

#include "../src/construction/models.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/locate.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.procedures_new.hpp"
#include "../src/fm_index/tree_new.hpp"

using namespace std;
using namespace sdsl;

void test_comparison(){
    int Lmin = 10;
    auto start = std::chrono::steady_clock::now();
    // auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');

    // auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto fm_idx = FmIndex(str);
    cout << "bwt $ cout: " << fm_idx.seq_cnt() << " wv tree took " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    
    SuffixArrayNode root = get_root(fm_idx);
    vector<MaximalRepeatAnnotation> repeats;
    SuffixArrayNode_NEW root_new = get_root_new(fm_idx);
    vector<MaximalRepeatAnnotation> repeats_new;
    
    start = std::chrono::steady_clock::now();
    // navigate_tree_new<MaximalRepeatAnnotation, get_repeat_annotations_new>(root_new, Lmin, fm_idx, repeats_new);
    navigate_tree_new<MaximalRepeatAnnotation, report_repr_sa>(root_new, Lmin, fm_idx, repeats_new);
    cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    int repr_cnt=0;
    int repr_distinct_cnt=0;
    auto repr_distinct = vector<bool>(str.size());
    for(int i=0; i<repeats_new.size(); i++){
        for(auto idx: repeats_new[i].repr_indexes){
            repr_cnt++;
            if(!repr_distinct[idx]){
                repr_distinct_cnt++;
                repr_distinct[idx]=true;
            }
        }
    }
    cout << "repeats reprs cnt: "<< repr_cnt <<", distinct count: " << repr_distinct_cnt  << endl;

    start = std::chrono::steady_clock::now();
    navigate_tree<MaximalRepeatAnnotation, get_repeat_annotations>(root, Lmin, fm_idx, repeats);
    cout << "original algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    assert(repeats.size()==repeats_new.size());
    cout << "repeats are same: " << repeats.size() << endl;
    
    for(int i=0; i<repeats.size(); i++){
        assert(repeats[i].repr_indexes.size()==repeats_new[i].repr_indexes.size());
    }
    cout << "repeats reprs are same" << endl;
}

void test_memory_reserved_stack(){
    std::vector<int> myVector;
    myVector.reserve(100);
    
    std::stack<int, std::vector<int>> myStack(myVector);
    cout << "stack size: "<< myStack.size() << endl;
    // Now push and pop without incurring dynamic allocations up to 100 elements
    myStack.push(1);
    cout << "stack size: "<< myStack.size() << endl;
    myStack.push(2);
    cout << "stack size: "<< myStack.size() << endl;
    myStack.pop();
    cout << "stack size: "<< myStack.size() << endl;
}

struct RepMcDummy{
    int size;
};

class AtomicStratumIdGenerator {
    int curr_rep_id;
    std::vector<RepMcDummy>* rep_mcs;

    std::atomic_flag lock = ATOMIC_FLAG_INIT; // Used as a spinlock for vector push_back

public:
    AtomicStratumIdGenerator(std::vector<RepMcDummy> &rep_mcs) : curr_rep_id(0), rep_mcs(&rep_mcs) {}

    int generateNumber() {
        int number = 0;
        // Spinlock to ensure that push_back is thread-safe
        while (lock.test_and_set(std::memory_order_acquire));
        number = curr_rep_id;
        curr_rep_id++;
        (*rep_mcs).push_back(RepMcDummy());
        lock.clear(std::memory_order_release);  // release lock

        return number;
    }
};

void worker(AtomicStratumIdGenerator& generator) {
    for (int i = 0; i < 10000; ++i) {
        int number = generator.generateNumber();
        // std::cout << "Generated number: " << number << " by thread: " << std::this_thread::get_id() << "\n";
        // std::this_thread::sleep_for(std::chrono::milliseconds(10));  // Just to space out the output
    }
}


void test_atomic_rep_id() {
    vector<RepMcDummy> rep_mcs;
    AtomicStratumIdGenerator generator(rep_mcs);

    auto start = std::chrono::steady_clock::now();

    std::vector<std::thread> threads;
    for (int i = 0; i < 12; ++i) {
        threads.emplace_back(worker, std::ref(generator));
    }

    for (auto& thread : threads) {
        thread.join();
    }
    
    std::cout << "size of rep mc:" << rep_mcs.size() << endl;

    cout << "elapsed: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

void test_infrastructural_components() {
    std::tuple<uint8_t, uint16_t> tuple;
    assert(sizeof(tuple)==4);
    // std::cout << "Size of tuple: " << sizeof(tuple) << " bytes" << std::endl;
}

void main_performance_new_tree() {
    // test_comparison();
    // test_memory_reserved_stack();
    // test_atomic_rep_id();
    test_infrastructural_components();
}
