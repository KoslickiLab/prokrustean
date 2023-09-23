#ifndef APPLICATION_UTIL_HPP_
#define APPLICATION_UTIL_HPP_
#include <algorithm>
#include <tuple>
#include "../prokrustean.hpp"

struct Interval{
    Pos from;
    Pos to;

    uint64_t size(){
        return to - from;
    }
    bool include(Interval other){
        return from <= other.from && other.to <= to;
    }

    vector<Interval> chop(unsigned int k){
        assert(size()>=k);
        vector<Interval> mers;
        for(Pos l=from; l<=to-k; l++){
            mers.push_back({l, l+k});
        }
        return mers;
    }

    bool is_suffix(Interval other){
        return from <= other.from && other.to == to;
    }

    bool is_prefix(Interval other){
        return from == other.from && other.to <= to;
    }

    bool operator==(const Interval& other)
    {
        return from == other.from && to == other.to;
    }
};

struct Occurrence: Interval{
    SeqId seq_id;

    Occurrence(SeqId seq_id, Pos from, Pos to){
        this->seq_id=seq_id;
        this->from=from;
        this->to=to;
    }
};

vector<Interval> get_gaps(uint64_t string_size, int degree, vector<Interval> reps){
    vector<Interval> gaps;
    if(reps.size()==0){
        gaps.push_back({0, string_size});
        return gaps;
    }
    
    for(int i=0; i<reps.size(); i++){
        assert(reps[i].size()>=degree);
        //prefix rep
        if(i==0 && reps[i].from>0){
            gaps.push_back({0, reps[i].from+degree-1});
        }
        //two reps -> need to have 
        if(i<reps.size()-1){
            Interval g = {reps[i].to-(degree-1), reps[i+1].from+(degree-1)};
            if(g.size()>=degree){
                gaps.push_back(g);
            }
        }
        //suffix rep
        if(i==reps.size()-1 && reps[i].to < string_size)
        {
            gaps.push_back({reps[i].to-(degree-1), string_size});
        }
    }
    return gaps;
}

vector<Interval> get_gaps_seq(Prokrustean pk, SeqId id, int degree){
    vector<Interval> rIntervals;
    for(auto r: pk.seqs[id].regions){
        Pos rpos = get<0>(r);
        StratumId rid = get<1>(r);
        uint64_t size = pk.stratums[rid].size;
        if(size>=degree){
            rIntervals.push_back({rpos, rpos + pk.stratums[rid].size});
        }
    }
    return get_gaps(pk.seqs[id].size, degree, rIntervals);
}

vector<Interval> get_gaps_rep(Prokrustean pk, StratumId id, int degree){
    vector<Interval> rIntervals;
    for(auto r: pk.stratums[id].regions){
        Pos rpos = get<0>(r);
        StratumId rid = get<1>(r);
        uint64_t size = pk.stratums[rid].size;
        if(size>=degree){
            rIntervals.push_back({rpos, rpos + pk.stratums[rid].size});
        }
    }
    return get_gaps(pk.stratums[id].size, degree, rIntervals);
}

vector<Occurrence> collect_rep_example_occurrences(Prokrustean pk){
    vector<optional<Occurrence>> occs(pk.stratums.size());
    stack<tuple<Stratum, Occurrence>> mcs_w_occ;
    for(auto mc: pk.seqs){
        mcs_w_occ.push(make_tuple(mc, Occurrence(mc.id, 0, mc.size)));
    }
    while(!mcs_w_occ.empty()){
        auto pair = mcs_w_occ.top();
        mcs_w_occ.pop();
        Stratum mc = get<0>(pair);
        Occurrence occ = get<1>(pair);
        for(auto r: mc.regions){
            Pos rpos = get<0>(r);
            StratumId rid = get<1>(r);
            if(occs[rid].has_value()){
                continue;
            }
            uint64_t size = pk.stratums[rid].size;
            auto r_occ = Occurrence(occ.seq_id, occ.from + rpos, occ.from + rpos + size);
            occs[rid] = r_occ;
            mcs_w_occ.push(make_tuple(pk.stratums[rid], r_occ));
        }
    }
    vector<Occurrence> non_optional;
    for(auto occ: occs){
        assert(occ.has_value()); //check exhaustive collection
        non_optional.push_back(occ.value());
    }
    return non_optional;
}
#endif