#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <future>
#include <random>
#include "const.cpp"	
#include "../src/prokrustean.hpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/ssa.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/construction/algorithms.stage1_tree_navigation.hpp"
#include "../src/sdsl/int_vector.hpp"
#include "../src/sdsl/rank_support_v.hpp"
#include "../src/sdsl/rrr_vector.hpp"
#include "../src/util/data.store.hpp"

using namespace std;
using namespace sdsl;

void test_sampling_push(){
    auto num_threads=6;
    auto sampling_factor=8;
    // auto str = WaveletString(PATH3_PERFORMANCE_SREAD_GUT_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH2_PERFORMANCE_SREAD_FULL_ROPEBWT2_BWT, '$');
    // auto str = WaveletString(PATH1_PERFORMANCE_SREAD_ROPEBWT2_BWT, '$');
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx=FmIndex(str);
    auto ssa=SampledSuffixArray(fm_idx, sampling_factor);
    
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_basic: " << "thread:" << num_threads << ", sampling:" << sampling_factor << endl;
    vector<future<void>> futures;
    // collect blocks
    auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    // consume blocks
    auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    ssa.validate();
}


void test_sampling_still_works_if_factor_exceeds_seq_length(){
    auto num_threads=6;
    auto sampling_factor=500;// only first positions of sequences may be sampled
    auto str=WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_idx=FmIndex(str);
    auto ssa=SampledSuffixArray(fm_idx, sampling_factor);
    
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_basic: " << "thread:" << num_threads << ", sampling:" << sampling_factor << endl;
    vector<future<void>> futures;
    // collect blocks
    auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    // consume blocks
    auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(ssa)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();
    cout << "finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;

    ssa.validate();
}

void test_sampling_works_the_same_for_sampling_factors(){
    int num_threads=6;
    auto str = WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_index = FmIndex(str);
    atomic<uint64_t> idx_generator;
    sdsl::bit_vector is_sampled_bv(str.size());
    auto ssa1=SampledSuffixArray(fm_index, 3);
    auto ssa2=SampledSuffixArray(fm_index, 7);
    auto ssa3=SampledSuffixArray(fm_index, 10);
    vector<SampledSuffixArray*> ssa_list;
    ssa_list.push_back(&ssa1);
    ssa_list.push_back(&ssa2);
    ssa_list.push_back(&ssa3);
    auto start = std::chrono::steady_clock::now();
    cout << "test_sampling_works_the_same_for_sampling_factors: " << "thread:" << num_threads << endl;

    vector<future<void>> futures;
    for(auto &ssa: ssa_list){
        start = std::chrono::steady_clock::now();
        // collect blocks
        auto func__sample_to_parallel = [](SampledSuffixArray &ssa) {while(ssa.sample_one_sequence_and_store()){}};
        for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__sample_to_parallel, ref(*ssa)));}
        for (auto &f : futures) {f.wait();}
        futures.clear();
        // consume blocks
        auto func__consume_block = [](SampledSuffixArray &ssa) {while(ssa.consume_one_block_and_release()){}};
        for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__consume_block, ref(*ssa)));}
        for (auto &f : futures) {f.wait();}
        futures.clear();
        cout << "one ssa finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    }

    auto func__compare_loc = [](vector<SampledSuffixArray*> &ssa_list, atomic<uint64_t> &idx_generator, uint64_t idx_limit) {
        uint64_t idx=0;
        optional<tuple<SeqId, Pos>> loc=nullopt;
        optional<tuple<SeqId, Pos>> other_loc=nullopt;
        while(true){
            idx=idx_generator.fetch_add(1);
            if(idx>=idx_limit){
                break;
            }
            loc=nullopt;
            other_loc=nullopt;
            for(auto &ssa: ssa_list){
                loc=ssa->get_location(idx);
                if(other_loc.has_value()){
                    assert(loc==other_loc);
                } else{
                    other_loc=loc;
                }
            }
        }
    };

    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__compare_loc, ref(ssa_list), ref(idx_generator), fm_index.size()));}
    for (auto &f : futures) {f.wait();}
    cout << "location comparison finished: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    futures.clear();
}

void test_temp(){
    int num_threads=6;
    auto str = WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_index = FmIndex(str);
    vector<int> lengths(fm_index.seq_cnt());
    auto start = std::chrono::steady_clock::now();
    string file_name="temp";
    cout << "start " << "thread:" << num_threads << endl;

    vector<future<void>> futures;
    auto func__compare_loc = [](FmIndex &fm_index, string &file_name, int thread_idx, int thread_cnt, vector<int> &lengths) {
        std::ofstream out_file(file_name+to_string(thread_idx)+".txt", std::ios::binary | std::ios::in | std::ios::out);
        uint64_t seq_cnt=fm_index.seq_cnt();
        for(int i=thread_idx; i<seq_cnt; i+=thread_cnt){
            SuffixArrayIdx L = i;
            SuffixArrayIdx F = fm_index.LF(L);
            uint64_t idx;
            Pos reverse_pos=0;
            vector<Pos> reverse_positions;
            vector<char> string;
            auto characters = fm_index.STRING->get_characters();
            while(F >= seq_cnt){
                L = F;
                F = fm_index.LF(L);
                out_file.write(reinterpret_cast<const char*>(&i), sizeof(i));
                out_file.write(reinterpret_cast<const char*>(&reverse_pos), sizeof(reverse_pos));
                reverse_pos++;
            }
            lengths[i]=reverse_pos;
        }
    };

    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__compare_loc, ref(fm_index), ref(file_name), i, num_threads, ref(lengths)));}
    for (auto &f : futures) {f.wait();}
    cout << "finished write: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    futures.clear();

    auto func__read = [](FmIndex &fm_index, string &file_name, int thread_idx, int thread_cnt, vector<int> &lengths) {
        std::ifstream in_file(file_name+to_string(thread_idx)+".txt", std::ios::binary);
        uint64_t seq_cnt=fm_index.seq_cnt();
        uint64_t idx;
        Pos pos=0;
        for(int i=thread_idx; i<seq_cnt; i+=thread_cnt){
            for(int p=0; p<lengths[i]; p++){
                in_file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
                in_file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
            }
        }
    };

    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__read, ref(fm_index), ref(file_name), i, num_threads, ref(lengths)));}
    for (auto &f : futures) {f.wait();}
    cout << "finished read: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}

struct Element {
    int x, y;
};

class FileStorage {
private:
    std::string baseFilename; // base name for all files
    size_t elementSize; // size in bytes for each element
    vector<std::ofstream> outs;
    vector<std::ifstream> ins;
    uint64_t sequence_length_per_file;
    uint64_t sequence_total_length;
    int numFiles; // total number of files
    vector<SpinLock> locks;

public:
    FileStorage(const std::string& base, int numFiles, uint64_t sequence_total_length)
        : baseFilename(base), numFiles(numFiles), sequence_total_length(sequence_total_length), elementSize(sizeof(Element)) {
            this->locks=vector<SpinLock>(numFiles);
            this->sequence_length_per_file=sequence_total_length/numFiles+1;
            cout << "this->sequence_length_per_file " << this->sequence_length_per_file<<endl;
        }
    
    void write_open() {
        for(int i=0; i<numFiles; i++){
            std::string filename = baseFilename + std::to_string(i) + ".dat";
            this->outs.push_back(std::ofstream(filename, std::ios::binary | std::ios::out));
        }
    }
    void write_close() {
        for(auto &out: this->outs){
            out.close();
        }
        this->outs.clear();
    }

    void read_open() {
        for(int i=0; i<numFiles; i++){
            std::string filename = baseFilename + std::to_string(i) + ".dat";
            this->ins.push_back(std::ifstream(filename, std::ios::binary));
        }
    }
    void read_close() {
        for(auto &in: this->ins){
            in.close();
        }
        this->ins.clear();
    }

    void write(uint64_t index, Element &e) {
        int file_idx = index/this->sequence_length_per_file;
        this->locks[file_idx].lock();
        // outs[file_idx].seekp(index%this->sequence_length_per_file * elementSize);
        outs[file_idx].write(reinterpret_cast<char*>(&e), elementSize);
        this->locks[file_idx].unlock();
        // cout << "store at " << file_idx << " , " << index%this->sequence_length_per_file << endl;
    }

    Element read(uint64_t index, Element &e) {
        int file_idx = index/this->sequence_length_per_file;
        ins[file_idx].seekg(index%this->sequence_length_per_file * elementSize);
        ins[file_idx].read(reinterpret_cast<char*>(&e), elementSize);
        // cout << "read at " << file_idx << " , " << index%this->sequence_length_per_file << endl;
        return e;
    }
};

void test_temp2() {
    int num_threads=5;
    auto str = WaveletString(PATH_SREAD_FULL_GRLBWT_BWT, '$');
    auto fm_index = FmIndex(str);
    auto start = std::chrono::steady_clock::now();
    vector<int> lengths(fm_index.seq_cnt());
    string file_name="temp";
    cout << "start " << "thread:" << num_threads << endl;
    FileStorage storage("element", 5*num_threads, fm_index.size());
    storage.write_open();

    vector<future<void>> futures;
    auto func__compare_loc = [](FmIndex &fm_index, FileStorage &storage, int thread_idx, int thread_cnt, vector<int> &lengths) {
        uint64_t seq_cnt=fm_index.seq_cnt();
        cout << "seq _cnt "<< seq_cnt << endl;
        vector<Pos> reverse_positions;
        Element e;
        for(int i=thread_idx; i<seq_cnt; i+=thread_cnt){
            reverse_positions.clear();
            SuffixArrayIdx L = i;
            SuffixArrayIdx F = fm_index.LF(L);
            uint64_t idx;
            Pos reverse_pos=0;
            while(F >= seq_cnt){
                L = F;
                F = fm_index.LF(L);
                e.x=i;
                e.y=reverse_pos;
                storage.write(L, e);
                reverse_pos++;
            }
            lengths[i]=reverse_pos;
        }
    };
    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__compare_loc, ref(fm_index), ref(storage), i, num_threads, ref(lengths)));}
    for (auto &f : futures) {f.wait();}
    futures.clear();

    cout << "finished write: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // vector<SuffixArrayIdx> suffix_indices;
    // vector<Pos> reverse_positions;
    // for(int i=0; i<fm_index.seq_cnt(); i++){
    //     suffix_indices.clear();
    //     reverse_positions.clear();

    //     SuffixArrayIdx L = i;
    //     SuffixArrayIdx F = fm_index.LF(L);
    //     uint64_t idx;
    //     Pos reverse_pos=0;
    //     while(F >= fm_index.seq_cnt()){
    //         L = F;
    //         F = fm_index.LF(L);
    //         suffix_indices.push_back(L);
    //         reverse_positions.push_back(reverse_pos);
    //         reverse_pos++;
    //     }
    //     for(int p=0; p<suffix_indices.size();p++){
    //         storage.write(suffix_indices[p], p, reverse_positions[p]);
    //     }
    //     seq_sizes.push_back(reverse_pos);
    // }
    
    storage.write_close();


    // storage.read_open();
    // auto func__read = [](FmIndex &fm_index, FileStorage &storage, int thread_idx, int thread_cnt, vector<int> &lengths) {
    //     uint64_t seq_cnt=fm_index.seq_cnt();
    //     cout << "seq _cnt "<< seq_cnt << endl;
    //     vector<Pos> reverse_positions;
    //     Element e;
    //     for(int i=thread_idx; i<seq_cnt; i+=thread_cnt){
    //         reverse_positions.clear();
    //         SuffixArrayIdx L = i;
    //         SuffixArrayIdx F = fm_index.LF(L);
    //         uint64_t idx;
    //         Pos reverse_pos=0;
    //         while(F >= seq_cnt){
    //             L = F;
    //             F = fm_index.LF(L);
    //             storage.read(L,e);
    //             reverse_pos++;
    //         }
    //         lengths[i]=reverse_pos;
    //     }
    // };
    // start = std::chrono::steady_clock::now();
    // for(int i=0; i<num_threads; i++){futures.push_back(std::async(std::launch::async, func__read, ref(fm_index), ref(storage), i, num_threads, ref(lengths)));}
    // for (auto &f : futures) {f.wait();}
    // futures.clear();
    // storage.read_close();
    // cout << "finished read: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // std::cout << "Read: (" << e.x << ", " << e.y << ")" << std::endl;
}


void main_performance_ssa() {
    // test_sampling_push();
    // test_sampling_still_works_if_factor_exceeds_seq_length();
    // test_sampling_works_the_same_for_sampling_factors();
    // test_temp();
    test_temp2();
}
