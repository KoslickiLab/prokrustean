#ifndef BASE_HPP_
#define BASE_HPP_

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

/* These types influence the output size */
typedef uint32_t StratumId; // 
typedef uint32_t SeqId; //
typedef uint32_t StratumOrSeqId; //
typedef uint16_t Pos; // Sequence max length
typedef uint16_t StratumSize; 
typedef uint32_t SequenceSize; 
typedef uint8_t BandIdx; // 

/* For prokrustean construction */
typedef uint16_t SuffixArrayIdx_InBlock;
typedef uint8_t CharId; //character ids -> lexicographical order of characters 
typedef uint64_t SuffixArrayIdx;
typedef vector<uint64_t> RankArray; 
/* For application */
typedef uint32_t UnitigId; // 
typedef uint32_t UnitigOrMaxUnitigId; // 
#endif