#ifndef BASE_HPP_
#define BASE_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

/* These types influence the size of outputs A LOT */
typedef uint32_t StratumId; // 
typedef uint32_t SeqId; //
typedef uint16_t Pos; // Sequence max length
typedef uint16_t StratumSize; 
typedef uint32_t SequenceSize; 

/* during constructions */
typedef uint16_t SuffixArrayIdx_InBlock;
//character ids -> lexicographical order of characters 
typedef uint8_t CharId;
typedef uint64_t SuffixArrayIdx;
//C array of fm index. Each index is CharId 
typedef vector<uint64_t> RankArray;
#endif