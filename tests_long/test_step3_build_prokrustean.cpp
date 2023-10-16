#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "const.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.stage1_stratum.hpp"
#include "../src/construction/algorithms.stage2_stratified_rgn.hpp"
#include "../src/construction/algorithms.stage1_tree.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"

using namespace std;
using namespace sdsl;


vector<StratumSize> get_random_stratum_sizes(int stra_cnt){
    auto stratum_length=30;
    auto sequence_length=10000;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> stratum_id_dist(1, sequence_length);
    std::uniform_int_distribution<uint32_t> stratum_len_dist(1, stratum_length);
    

    vector<StratumSize> sizes;
    for(int i=0; i<stra_cnt; i++){
        sizes.push_back(stratum_len_dist(gen));
    }
    return sizes;
}

vector<ProjectedStratifiedRegion> get_random_projected_regions(int regions_cnt, vector<StratumSize>* stratum_sizes){
    auto sequence_length=100;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> stratum_id_dist(0, stratum_sizes->size());
    std::uniform_int_distribution<uint32_t> seq_len_dist(0, sequence_length);

    vector<ProjectedStratifiedRegion> regions;
    for(int i=0; i<regions_cnt; i++){
        regions.push_back(ProjectedStratifiedRegion(stratum_id_dist(gen), seq_len_dist(gen), seq_len_dist(gen)%2));
    }
    return regions;
}


// void test_basic_step3_operation(){
//     auto pos_cnt=10;
//     auto stra_cnt=3;
    
//     StratificationWorkSpace workspace; 
//     Prokrustean prokrustean;
    
//     auto stratum_sizes=get_random_stratum_sizes(stra_cnt);
//     auto regions=get_random_projected_regions(pos_cnt*stra_cnt, &stratum_sizes);
//     //reset function
//     workspace.update_contexts_by_seq(0, regions, stratum_sizes);
//     assert(workspace.consume_refs.size()==workspace.seq_annot.position_annots.size());
    

//     workspace.seq_annot.position_annots[5].regions[2].is_primary=true;
//     auto primary=workspace.first_primary();
//     assert(primary.has_value());
//     assert(primary.value().pidx==5 && primary.value().sidx==2);
//     // earlier 
//     workspace.seq_annot.position_annots[6].regions[2].is_primary=true;
//     primary=workspace.first_primary();
//     assert(primary.has_value());
//     assert(primary.value().pidx==6 && primary.value().sidx==2);
//     // earlier 
//     workspace.seq_annot.position_annots[6].regions[1].is_primary=true;
//     primary=workspace.first_primary();
//     assert(primary.has_value());
//     assert(primary.value().pidx==6 && primary.value().sidx==1);
    
//     // next 
//     primary=workspace.first_primary();
//     primary=workspace.next_primary(primary.value());
//     assert(primary.value().pidx==6 && primary.value().sidx==2);
//     primary=workspace.next_primary(primary.value());
//     assert(primary.value().pidx==5 && primary.value().sidx==2);
//     primary=workspace.next_primary(primary.value());
//     assert(!primary.has_value());
// }

struct region {
    int pos;
    int stratum_id;
    int size;   
    bool is_primary;   
    region(int pos, int stratum_id, int size, bool is_primary):
    pos(pos), stratum_id(stratum_id), size(size), is_primary(is_primary)
    {}
};

void _generate_step3_scenario(vector<region> rgns, StratificationWorkSpace &workspace, Prokrustean &prokrustean){
    vector<ProjectedStratifiedRegion> projected_regions;
    prokrustean.seqs.clear();
    prokrustean.seqs.push_back(Sequence());

    int max=0;
    for(auto &r: rgns){
        max=max<r.stratum_id? r.stratum_id: max; 
        projected_regions.push_back(ProjectedStratifiedRegion(r.stratum_id, r.pos, r.is_primary));
    }
    prokrustean.stratums__size.clear();
    prokrustean.stratums__size.resize(max+1);
    prokrustean.stratums__region.clear();
    prokrustean.stratums__region.resize(max+1);
    prokrustean.stratums__region_cnt.clear();
    prokrustean.stratums__region_cnt.resize(max+1);
    
    for(auto &r: rgns){
        prokrustean.stratums__size[r.stratum_id]=r.size;
    }
    
    workspace.update_contexts_by_seq(0, projected_regions, prokrustean.stratums__size);

    for(auto &pos: workspace.seq_annot.position_annots){
        cout << "pos: " << pos.pos << endl;
        for(auto &region: pos.regions){
            cout << "region: " << "stratum id:"<< region.stratum_id << " stratum size:" << region.stratum_size << " is primary: "<< region.is_primary << endl;
        }
    }
}

void test_step3_scenarios(){
    StratificationWorkSpace workspace;
    Prokrustean prokrustean;
    
    // scenario: empty. one region
    auto regions= vector<region>({
        region(1, 1, 3, true) // int pos, int stratum_id, int size, bool is_primary
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.seqs[0].region_cnt==1);
    
    // leftmost relationship
    prokrustean=Prokrustean();
    regions= vector<region>({
        region(1, 1, 3, false), // int pos, int stratum_id, int size, bool is_primary
        region(1, 2, 5, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==0);
    assert(prokrustean.stratums__region_cnt[2]==1);
    assert(prokrustean.stratums__size[2]==5);
    

    // rightmost relationship
    regions= vector<region>({
        region(3, 1, 3, true), // int pos, int stratum_id, int size, bool is_primary
        region(1, 2, 5, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==0);
    assert(prokrustean.stratums__region_cnt[2]==1);
    assert(prokrustean.stratums__size[2]==5);


    // middle relationship
    regions= vector<region>({
        region(3, 1, 3, true), // int pos, int stratum_id, int size, bool is_primary
        region(1, 2, 7, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[2]==1);


    // middle relationship2 layered
    regions= vector<region>({
        region(3, 1, 3, false), // int pos, int stratum_id, int size, bool is_primary
        region(3, 2, 4, false),
        region(1, 3, 7, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==0);
    assert(prokrustean.stratums__region_cnt[2]==0);
    assert(prokrustean.stratums__region_cnt[3]==1);


    // double inclusion relationship
    regions= vector<region>({
        region(3, 1, 3, true), // int pos, int stratum_id, int size, bool is_primary
        region(1, 2, 5, true),
        region(3, 3, 7, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==0);
    assert(prokrustean.stratums__region_cnt[2]==1);
    assert(prokrustean.stratums__region_cnt[3]==1);
    assert(prokrustean.seqs[0].region_cnt==2);

    // double inclusion relationship2 - largest one includes all
    regions= vector<region>({
        region(3, 1, 3, true), // int pos, int stratum_id, int size, bool is_primary
        region(2, 2, 4, true),
        region(3, 3, 7, true),
        region(1, 4, 10, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==0);
    assert(prokrustean.stratums__region_cnt[2]==1);
    assert(prokrustean.stratums__region_cnt[3]==1);
    assert(prokrustean.stratums__region_cnt[4]==2);
    assert(prokrustean.seqs[0].region_cnt==1);

    // one includes over two
    regions= vector<region>({
        region(1, 1, 10, true), // int pos, int stratum_id, int size, bool is_primary
        region(2, 2, 2, true),
        region(2, 3, 4, true),
        region(3, 4, 5, true)
    });
    _generate_step3_scenario(regions, workspace, prokrustean);
    build_prokrustean(workspace, prokrustean);
    assert(prokrustean.stratums__region_cnt[1]==2);
    assert(prokrustean.stratums__region_cnt[2]==0);
    assert(prokrustean.stratums__region_cnt[3]==1);
    assert(prokrustean.stratums__region_cnt[4]==0);
    assert(prokrustean.seqs[0].region_cnt==1);
}

void main_performance_build_prokrustean() {
    test_step3_scenarios();
}
