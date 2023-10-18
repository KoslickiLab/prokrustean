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
typedef uint16_t Pos; 
typedef uint16_t SequenceSize; 
/* Has to be under possible stratum length */
typedef uint16_t StratumSize; 
/* Can covering region be over 255? meaning stratified regions are at least 122? */ 
/*I doubt about it even if the dataset is very largein normal cases */
typedef uint8_t CoveringRegionIdx; 

/* For prokrustean construction */
typedef uint16_t SuffixArrayIdx_InBlock;
typedef uint8_t CharId; //character ids -> lexicographical order of characters 
typedef uint64_t SuffixArrayIdx;
typedef vector<uint64_t> RankArray; 
/* For application */
typedef uint32_t UnitigId; // 
typedef uint32_t UnitigOrMaxUnitigId; // 
#endif