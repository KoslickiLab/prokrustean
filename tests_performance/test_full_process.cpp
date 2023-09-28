// #include <fstream>
// #include <vector>
// #include <cassert>
// #include <iostream>
// #include <random>
// #include <future>
// #include "util.cpp"	
// #include "../src/construction/algorithms.step1_project_stratums.hpp"
// #include "../src/construction/algorithms.step2_build_prokrustean.hpp"
// #include "../src/prokrustean.hpp"
// #include "../src/fm_index/index.hpp"
// #include "../src/fm_index/ssa.hpp"
// #include "../src/fm_index/string.sdsl.hpp"
// #include "../src/fm_index/tree_new.hpp"
// #include "../src/sdsl/int_vector.hpp"
// #include "../src/sdsl/rank_support_v.hpp"
// #include "../src/sdsl/rrr_vector.hpp"

// using namespace std;
// using namespace sdsl;

////////////////////////////////////////////////////////////////////////////////
/* Use ssa. for archiving                                                     */
////////////////////////////////////////////////////////////////////////////////

// void test_step1_push(){
//     int Lmin=30;
//     auto num_threads=12;
//     auto sampling_factor=8;

//     // auto str=WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
//     WaveletString str(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
//     // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
//     FmIndex fm_idx(str);
//     SampledSuffixArray ssa(fm_idx, sampling_factor);
//     Prokrustean prokrustean;
    
//     auto start = std::chrono::steady_clock::now();
//     vector<future<void>> futures;
//     // collect blocks
//     auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
//     for (auto &f : futures) {f.wait();}
//     futures.clear();
//     // consume blocks
//     auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
//     for (auto &f : futures) {f.wait();}
//     futures.clear();
//     cout << "finished sampling: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
//     SuffixArrayNode_NEW root = get_root_new(fm_idx);
//     StratumProjectionOutput output(prokrustean, fm_idx.seq_cnt());
    
//     fm_idx.set_sampled_suffix_array(ssa);
//     atomic<int> idx_gen;
//     start = std::chrono::steady_clock::now();
//     vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
//     auto func__navigate = [](vector<SuffixArrayNode_NEW> &roots, FmIndex &fm_idx, int Lmin, StratumProjectionOutput &output, atomic<int> &idx_gen) {
//         while(true){
//             auto idx = idx_gen.fetch_add(1);
//             if(idx>=roots.size()){
//                 break;
//             } else {
//                 navigate_maximals<StratumProjectionOutput, report_representative_locations>(roots[idx], Lmin, fm_idx, output);
//             }
//         }
//     };
    
//     // futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));}
//     for (auto &f : futures) {f.wait();}
    
//     cout << "new algorithm: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
// }

// void test_full_process_push(){
//     int Lmin=20;
//     auto num_threads=12;
//     auto sampling_factor=16;
//     int sleep1=0;
//     int sleep2=0;
//     auto start = std::chrono::steady_clock::now();

//     // WaveletString str(PATH1_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
//     // auto str=WaveletString(PATH2_PERFORMANCE_SREAD_FULL_GRLBWT_BWT, '$');
//     // auto str=WaveletString(PATH5_PERFORMANCE_SREAD_GRLBWT_BWT, '$');
//     auto str=WaveletString(PATH3_PERFORMANCE_SREAD_GUT_GRLBWT_BWT, '$');
//     cout << "wavelet string: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

//     if(sleep1>0) std::cout << "1. Wavelete string. sleeping... " << sleep1 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep1));
//     if(sleep1>0) std::cout << "2. start " << sleep1 << std::endl;
//     FmIndex fm_idx(str);
//     SampledSuffixArray ssa(fm_idx, sampling_factor);
//     Prokrustean prokrustean;
    
//     start = std::chrono::steady_clock::now();
//     vector<future<void>> futures;
//     // collect blocks
//     auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
//     for (auto &f : futures) {f.wait();}
//     futures.clear();
    
//     if(sleep1>0) std::cout << "2. sample raw. sleeping... " << sleep1 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep1));
//     if(sleep1>0) std::cout << "3. start " << sleep1 << std::endl;
    
//     // consume blocks
//     auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
//     for (auto &f : futures) {f.wait();}
//     futures.clear();
//     fm_idx.set_sampled_suffix_array(ssa);
//     cout << "finished sampling: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    
//     if(sleep1>0) std::cout << "3. sample consumed. sleeping... " << sleep1 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep1));
//     if(sleep1>0) std::cout << "4. start " << sleep1 << std::endl;

//     SuffixArrayNode_NEW root = get_root_new(fm_idx);
//     StratumProjectionOutput workspace(prokrustean, fm_idx.seq_cnt());
    
//     atomic<int> idx_gen;
//     start = std::chrono::steady_clock::now();
//     vector<SuffixArrayNode_NEW> roots = collect_nodes(root, fm_idx, 3);
    
//     auto func__navigate = [](vector<SuffixArrayNode_NEW> &roots, FmIndex &fm_idx, int Lmin, StratumProjectionOutput &output, atomic<int> &idx_gen) {
//         while(true){
//             auto idx = idx_gen.fetch_add(1);
//             if(idx>=roots.size()){
//                 break;
//             } else {
//                 navigate_maximals<StratumProjectionOutput, report_representative_locations>(roots[idx], Lmin, fm_idx, output);
//             }
//         }
//     };
    
//     // futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(output), ref(idx_gen)));
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__navigate, ref(roots), ref(fm_idx), Lmin, ref(workspace), ref(idx_gen)));}
//     for (auto &f : futures) {f.wait();}
    
//     cout << "step1 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
//     if(sleep1>0) std::cout << "4. navigated sleeping... " << sleep2 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep1));
    
//     fm_idx.dispose();
    
//     if(sleep2>0) std::cout << "5. fm indx dropped sleeping... " << sleep2 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep2));
    
//     ssa.set_sequence_lengths(prokrustean);
//     ssa.dispose();
    
//     if(sleep2>0) std::cout << "6. set_and_drop_sequences sleeping... " << sleep2 << std::endl;
//     std::this_thread::sleep_for(std::chrono::seconds(sleep2));
    
//     workspace.prokrustean->stratums__region.resize(workspace.prokrustean->stratums__size.size());
//     workspace.prokrustean->stratums__region_cnt.resize(workspace.prokrustean->stratums__size.size(), 0);
    

//     auto func__build = [](Prokrustean &prokrustean, StratumProjectionOutput &output, atomic<int> &idx_gen) {
//         StratificationWorkSpace workspace;
//         while(true){
//             auto idx = idx_gen.fetch_add(1);
//             if(idx>=prokrustean.seqs.size()){
//                 break;
//             }
//             workspace.update_contexts_by_seq(idx, output.sequence_regions[idx], prokrustean.stratums__size);
//             build_prokrustean(workspace, prokrustean);
//         }
//     };

//     start = std::chrono::steady_clock::now();
//     atomic<int> seq_id_iter;
//     for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__build, ref(prokrustean), ref(workspace), ref(seq_id_iter)));}
//     for (auto &f : futures) {f.wait();}
//     cout << "step2 finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
// }

// void main_performance_full() {
//     // test_step1_push();
//     test_full_process_push();
// }
