
#ifndef CONSTRUCTION_ALGO_HPP_
#define CONSTRUCTION_ALGO_HPP_
#include <algorithm>
#include <numeric>
#include <thread>
#include <future>
#include <list>
#include "algorithms.stage1_stratification.hpp"
#include "algorithms.stage2_prokrustean.hpp"

using namespace std;

struct Configuration{
    // 
    bool collect_ext_count=false;
};

void construct_prokrustean(FmIndex &fm_idx, Prokrustean &prokrustean, uint64_t Lmin=1, ProkrusteanEnhancement* opt=nullptr){
    auto start = std::chrono::steady_clock::now();
    cout << "step1 start" << endl;
    SuffixArrayNode root = get_root(fm_idx);
    StratumProjectionWorkspace workspace_step1(prokrustean, fm_idx, opt);
    navigate_maximals<StratumProjectionWorkspace, report_representative_locations>(root, Lmin, fm_idx, workspace_step1);
    cout << "step1 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    start = std::chrono::steady_clock::now();
    cout << "step2 start" << endl;
    workspace_step1.setup_prokrustean();
    StratificationWorkSpace workspace_step2;
    uint64_t cnt=prokrustean.sequence_count();
    for(uint64_t i=0; i<cnt; i++){
        workspace_step2.update_contexts_for_seq(i, fm_idx, workspace_step1, prokrustean.stratums__size);
        build_prokrustean(workspace_step2, prokrustean);
    }
    cout << "step2 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    cout << "stratum: " << prokrustean.stratum_count() << ", stratified regions new: "<< prokrustean.get_cardinality() << endl;
}

#endif