#ifndef BASE_HPP_
#define BASE_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

/* These types influence the size of outputs A LOT */

typedef uint32_t StratumId; // 4,294,967,295
typedef uint32_t SeqId; // 4,294,967,295
typedef uint32_t Pos; // 4,294,967,295 == Sequence max length
typedef uint16_t StratumSize; //65,535



class SpinLock {
    std::atomic_flag flag = ATOMIC_FLAG_INIT;

public:
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire)) {}
    }

    void unlock() {
        flag.clear(std::memory_order_release);
    }
};

#endif