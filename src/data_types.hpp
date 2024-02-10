#ifndef BASE_HPP_
#define BASE_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

/* These types influence the output size */

/* Has to be under no. of sequences */
typedef uint32_t SeqId; //
/* Has to be under maximum no. of stratums(maximal repeats) */
typedef uint32_t StratumId; // 
/* Has to be at least (SeqId or StratumId) */
typedef uint32_t StratumOrSeqId; //
/* Has to be under possible sequence max length. Has to be synced together */
typedef uint32_t Pos; 
typedef uint32_t SequenceSize; 
/* Has to be under possible stratum length */
typedef uint32_t StratumSize; 
/* Can covering region be over 255? meaning stratified regions are at least 122? */ 
/*I doubt about it even if the dataset is very largein normal cases */
typedef uint32_t CoveringRegionIdx; 
typedef uint32_t CoveringRegionCount; 
//Alphabet
typedef uint8_t CharId; 
typedef uint8_t CharCount; 
// frequency
typedef uint32_t FrequencyCount; 

/* For prokrustean construction */
typedef uint16_t SuffixArrayIdx_InBlock;
typedef uint64_t SuffixArrayIdx;
typedef std::vector<uint64_t> RankArray; 
/* For application */
typedef uint32_t UnitigId; // 
typedef uint32_t UnitigOrMaxUnitigId; // 


class SpinLock {
    std::atomic_flag flag = ATOMIC_FLAG_INIT;

public:
    SpinLock(){}
    void lock() {
        while (flag.test_and_set(std::memory_order_acquire)) {}
    }

    void unlock() {
        flag.clear(std::memory_order_release);
    }
};
#endif