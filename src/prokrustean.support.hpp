#ifndef PROKRUSTEAN_SUPPORT_HPP_
#define PROKRUSTEAN_SUPPORT_HPP_
#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "prokrustean.hpp"

using namespace std;

struct ProkrusteanExtension {
    Prokrustean& prokrustean;

    vector<SeqId> stratum_sample_occ_seq_id;
    vector<Pos> stratum_sample_occ_pos;

    vector<uint8_t> stratum_left_ext_count;
    vector<uint8_t> stratum_right_ext_count;

    // construction context.
    bool collect_left_right_extensions=false;
    
    // CDBG
    // for counting maximal unitigs of single k
    vector<optional<bool>> stratum__left_bound_stratum_is_branching;
    vector<bool> stratum__is_start_of_maximal_unitig;
    vector<bool> sequence__is_start_of_maximal_unitig;

    vector<SpinLock> stratum_locks;
    int stratum_lock_scale;
    
    ProkrusteanExtension(Prokrustean& prokrustean):prokrustean(prokrustean) {}

    void set_stratum_locks(int thread_cnt){
        this->stratum_lock_scale=thread_cnt*1000;
        this->stratum_locks=vector<SpinLock>(this->stratum_lock_scale);
    }
    void lock_stratum(StratumId id){
        this->stratum_locks[id%this->stratum_lock_scale].lock();
    }
    void unlock_stratum(StratumId id){
        this->stratum_locks[id%this->stratum_lock_scale].unlock();
    }
};

bool no_stratified_region_in_stra(StratumId id, int k, ProkrusteanExtension &ext){
    if(ext.prokrustean.stratums__region_cnt[id]==0)
    return true;
    for(int i=0; i< ext.prokrustean.stratums__region_cnt[id]; i++){
        if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][i].stratum_id]>=k){
            return false;
        }
    }
    return true;
}

// bool no_stratified_region_in_seq(SeqId id, int k, ProkrusteanExtension &ext){
//     if(ext.prokrustean.sequences__region_cnt[id]==0)
//     return true;
//     for(int i=0; i< ext.prokrustean.sequences__region_cnt[id]; i++){
//         if(ext.prokrustean.stratums__size[ext.prokrustean.sequences__region[id][i].stratum_id]>=k-1){
//             return false;
//         }
//     }
//     return true;
// }

uint8_t first_stratified_region_idx_in_stra(StratumId id, int k, ProkrusteanExtension &ext){
    for(int i=0; i< ext.prokrustean.stratums__region_cnt[id]; i++){
        if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][i].stratum_id]>=k){
            return i;
        }
    }
    // no stratified region over k-1 found
    assert(false);
}

uint8_t last_stratified_region_idx_in_stra(StratumId id, int k, ProkrusteanExtension &ext){
    for(int i=ext.prokrustean.stratums__region_cnt[id]-1; i>=0; i--){
        if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][i].stratum_id]>=k){
            return i;
        }
    }
    // no stratified region over k-1 found
    assert(false);
}

uint8_t next_stratified_region_idx_in_stra(StratumId id, uint8_t idx, int k, ProkrusteanExtension &ext){
    idx++;
    while(idx<ext.prokrustean.stratums__region_cnt[id]){
        if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][idx].stratum_id]>=k){
            return idx;
        }
        idx++;
    }
    // no stratified region over k-1 found
    assert(false);
}

uint8_t prev_stratified_region_idx_in_stra(StratumId id, uint8_t idx, int k, ProkrusteanExtension &ext){
    idx--;
    while(idx>=0){
        if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][idx].stratum_id]>=k){
            return idx;
        }
        idx--;
    }
    // no stratified region over k-1 found
    assert(false);
}
// uint8_t first_stratified_region_idx_in_seq(SeqId id, int k, ProkrusteanExtension &ext){
//     for(int i=0; i< ext.prokrustean.sequences__region_cnt[id]; i++){
//         if(ext.prokrustean.stratums__size[ext.prokrustean.stratums__region[id][i].stratum_id]>=k-1){
//             return i;
//         }
//     }
//     // no stratified region over k-1 found
//     assert(false);
// }

void setup_stratum_example_occ(ProkrusteanExtension &ext){
    ext.stratum_sample_occ_seq_id.resize(ext.prokrustean.stratum_count);
    ext.stratum_sample_occ_pos.resize(ext.prokrustean.stratum_count);
    vector<bool> visits(ext.prokrustean.stratum_count);
    std::stack<StratumId> stratum_stack;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        for(auto &rgn: ext.prokrustean.get_sequence(i).s_edges){
            ext.stratum_sample_occ_seq_id[rgn.stratum_id]=i;
            ext.stratum_sample_occ_pos[rgn.stratum_id]=rgn.from;
            visits[rgn.stratum_id]=true;

            stratum_stack.push(rgn.stratum_id);
            while(!stratum_stack.empty()){
                auto stratum_id=stratum_stack.top();
                stratum_stack.pop();
                for(auto &c_rgn: ext.prokrustean.get_stratum(stratum_id).s_edges){
                    if(visits[c_rgn.stratum_id]) 
                    continue;
                    SeqId seq_id=ext.stratum_sample_occ_seq_id[stratum_id];
                    Pos rel_pos=ext.stratum_sample_occ_pos[stratum_id]+c_rgn.from;
                    ext.stratum_sample_occ_seq_id[c_rgn.stratum_id]=seq_id;
                    ext.stratum_sample_occ_pos[c_rgn.stratum_id]=rel_pos;
                    visits[c_rgn.stratum_id]=true;
                    stratum_stack.push(c_rgn.stratum_id);
                }
            }
        }
    }

}

void _set_stratum_example_occ_for_seq(SeqId id, ProkrusteanExtension &ext, vector<uint8_t> &visits){
    std::stack<StratumId> stratum_stack;
    for(auto &rgn: ext.prokrustean.get_sequence(id).s_edges){
        ext.stratum_sample_occ_seq_id[rgn.stratum_id]=id;
        ext.stratum_sample_occ_pos[rgn.stratum_id]=rgn.from;
        ext.lock_stratum(rgn.stratum_id);
        if(visits[rgn.stratum_id]==1){
            ext.unlock_stratum(rgn.stratum_id);
            continue;
        } else {
            visits[rgn.stratum_id]=1;
            ext.unlock_stratum(rgn.stratum_id);
        }

        stratum_stack.push(rgn.stratum_id);
        while(!stratum_stack.empty()){
            auto stratum_id=stratum_stack.top();
            stratum_stack.pop();
            for(auto &c_rgn: ext.prokrustean.get_stratum(stratum_id).s_edges){
                ext.lock_stratum(c_rgn.stratum_id);
                if(visits[c_rgn.stratum_id]==1){
                    ext.unlock_stratum(c_rgn.stratum_id);
                    continue;
                } else {
                    visits[c_rgn.stratum_id]=1;
                    ext.unlock_stratum(rgn.stratum_id);
                }

                SeqId seq_id=ext.stratum_sample_occ_seq_id[stratum_id];
                Pos rel_pos=ext.stratum_sample_occ_pos[stratum_id]+c_rgn.from;
                ext.stratum_sample_occ_seq_id[c_rgn.stratum_id]=seq_id;
                ext.stratum_sample_occ_pos[c_rgn.stratum_id]=rel_pos;
                stratum_stack.push(c_rgn.stratum_id);
            }
        }
    }
}

void setup_stratum_example_occ_parallel(ProkrusteanExtension &ext, int thread_cnt){
    ext.set_stratum_locks(thread_cnt);
    ext.stratum_sample_occ_seq_id.resize(ext.prokrustean.stratum_count);
    ext.stratum_sample_occ_pos.resize(ext.prokrustean.stratum_count);
    vector<uint8_t> visits(ext.prokrustean.stratum_count, 0);
    vector<future<void>> futures;
    atomic<int> seq_idx_gen;

    auto func_ = [](vector<uint8_t> &visits, ProkrusteanExtension &ext, atomic<int> &seq_idx_gen) {
        while(true){
            auto idx = seq_idx_gen.fetch_add(1);
            if(idx>=ext.prokrustean.sequence_count)
            break;
            _set_stratum_example_occ_for_seq(idx, ext, visits);
            cout << idx << endl;
        }
    };
    // for(int i=0; i<thread_cnt; i++){futures.push_back(std::async(std::launch::async, func_, ref(visits), ref(ext), ref(seq_idx_gen)));}
    for (auto &f : futures) {f.wait();}
}

// void fill_stratum_left_bound_single_right_extension(int k){
//     // for(int i=0; i<prokrustean.sequence_count; i++){
//     //     Sequence seq = prokrustean.get_sequence(i);
//     //     seq.print();
//     // }
//     // for(int i=0; i<prokrustean.stratum_count; i++){
//     //     Stratum stra = prokrustean.get_stratum(i);
//     //     stra.print();
//     // }
//     this->stratum__left_bound_stratum_is_branching=vector<optional<bool>>(prokrustean.stratum_count);
//     vector<uint8_t> indices;
//     /* single thread version of filling extensions */
//     for(int i=0; i<prokrustean.sequence_count; i++){
//         Sequence seq = prokrustean.get_sequence(i);
//         // seq.print();
//         seq.get_valid_indices(k, indices);
//         // cout<< "seq " << i << " indices " << indices.size() << endl;
//         // has to pair two
//         for(int r=0; r+1<indices.size(); r++){
//             fill_for_two_(seq.s_regions[indices[r]], seq.s_regions[indices[r+1]], k);
//         }
//     }
//     for(int i=0; i<prokrustean.stratum_count; i++){
//         Stratum stra = prokrustean.get_stratum(i);
//         stra.print();
//         stra.get_valid_indices(k, indices);
//         // has to pair two
//         for(int r=0; r+1<indices.size(); r++){
//             fill_for_two_(stra.s_regions[indices[r]], stra.s_regions[indices[r+1]], k);
//         }
//     }

//     this->stratum__is_start_of_maximal_unitig=vector<bool>(prokrustean.stratum_count, false);
//     this->sequence__is_start_of_maximal_unitig=vector<bool>(prokrustean.sequence_count, false);
//     // count
//     int cnt=0;
//     for(int i=0; i<prokrustean.sequence_count; i++){
//         Sequence seq = prokrustean.get_sequence(i);
//         if(seq.size<k){
//             continue;
//         }
        
//         if(seq.s_regions.size()>0 && seq.s_regions[0].size()>=k){
//             this->sequence__is_start_of_maximal_unitig[i]=false;
//         } else {
//             // tip location: k-reflectum will have single occurrence (in sequence) 
//             cout << "seq tip " << i << " found" << endl;
//             this->sequence__is_start_of_maximal_unitig[i]=true;
//             cnt++;
//         }
//     }

//     for(int i=0; i<prokrustean.stratum_count; i++){
//         Stratum stra = prokrustean.get_stratum(i);
//         if(stra.size<k){
//             continue;
//         }
//         if(stra.s_regions.size()>0 && stra.s_regions[0].size()>=k){
//             this->stratum__is_start_of_maximal_unitig[i]=false;
//         } else if(prokrustean_optional.stratum_left_ext_count[i]!=1) {
//             // if count is 0, then tip location.
//             // if count is >1, then branch convergence.
//             this->stratum__is_start_of_maximal_unitig[i]=true;
//             cnt++;
//         } else if(!stratum__left_bound_stratum_is_branching[i].has_value()){
//             // 1. maybe stratified region is on pos 0
//             // 2. maybe ancestors of every occurrence is leftmost but not pos 0 -> always reflectum defines the maximal
//             this->stratum__is_start_of_maximal_unitig[i]=false;
//         } else if(stratum__left_bound_stratum_is_branching[i].value()){
//             this->stratum__is_start_of_maximal_unitig[i]=true;
//             cnt++;
//         } else {
//             this->stratum__is_start_of_maximal_unitig[i]=false;
//         }
//     }

//     vector<Region> spectrum;
//     for(int i=0; i<prokrustean.sequence_count; i++){
//         Sequence seq = prokrustean.get_sequence(i);
//         prokrustean.get_spectrum(seq, k, spectrum);
//         for(int r=0; r< spectrum.size(); r++){
//             if(r>0 && spectrum[r].is_reflected){
//                 // must be stratified
//                 if(1<prokrustean_optional.stratum_right_ext_count[spectrum[r-1].stratum_id]){
//                     cnt++;
//                     cout << "found reflectum after branching stratum id: " << spectrum[r-1].stratum_id << " rgn: " << r << endl;
//                 }
//             }
//         }
//     }

//     for(int i=0; i<prokrustean.stratum_count; i++){
//         Stratum stra = prokrustean.get_stratum(i);
//         prokrustean.get_spectrum(stra, k, spectrum);
//         for(int r=0; r< spectrum.size(); r++){
//             if(r>0 && spectrum[r].is_reflected){
//                 // must be stratified
//                 if(1<prokrustean_optional.stratum_right_ext_count[spectrum[r-1].stratum_id]){
//                     cnt++;
//                     cout << "found reflectum after branching stratum id: " << spectrum[r-1].stratum_id << " rgn: " << r << endl;
//                 }
//             }
//         }
//     }

//     cout << "maximal unitig counted: " << cnt << endl;
//     ran__fill_stratum_left_bound_single_right_extension=true;
// }

// void fill_for_two_(StratifiedRegion& l_rgn, StratifiedRegion& r_rgn, int k){
//     // go deepest to rightmost of left and leftmost of right.
//     // if k is fixed, there is no way but go into the deepest stratum such that the rightmost/leftmost is not (stratified and >=k)
//     StratumId l_id= get_deepest_stratum(l_rgn, k, true);
//     // note: since it is left side region of a pair, right ext is always 1 or more
//     bool left_desendant_is_branching=prokrustean_optional.stratum_right_ext_count[l_id]>1;
//     cout << "left descendant: " << left_desendant_is_branching << endl;
//     // StratumId r_id= get_deepest_stratum(r_rgn, k, true);
//     // bool left_descendant_is_diverging=prokrustean_optional.stratum_right_ext_count[l_id]!=1;
//     auto r_id=r_rgn.stratum_id;
//     while(true){
//         //already computed 
//         if(stratum__left_bound_stratum_is_branching[r_id].has_value()){
//             break;
//         }
//         auto child_stratum=prokrustean.get_stratum(r_id);
//         //leftmost shares edge
//         if(child_stratum.s_regions.size()>0){
//             auto& leftmost_rgn = child_stratum.s_regions[0]; 
//             if(leftmost_rgn.from==0 && leftmost_rgn.size()>=k){
//                 r_id=leftmost_rgn.stratum_id;
//                 continue;
//             }
//         }    
//         stratum__left_bound_stratum_is_branching[r_id]=left_desendant_is_branching;
//         break;
//     }

// }

// StratumId get_deepest_stratum(StratifiedRegion& rgn, int k, bool is_direction_left){
//     auto stratum_id=rgn.stratum_id;
//     if(is_direction_left){
//         while(true){
//             auto child_stratum=prokrustean.get_stratum(stratum_id);
//             //rightmost of the left shares edge
//             if(child_stratum.s_regions.size()>0){
//                 auto& rightmost_rgn = child_stratum.s_regions[child_stratum.s_regions.size()-1]; 
//                 if(rightmost_rgn.to==child_stratum.size && rightmost_rgn.size()>=k){
//                     stratum_id=rightmost_rgn.stratum_id;
//                     continue;
//                 }
//             }
//             // stop
//             break;
//         }
//     } else {
//         while(true){
//             auto child_stratum=prokrustean.get_stratum(stratum_id);
//             //leftmost shares edge
//             if(child_stratum.s_regions.size()>0){
//                 auto& leftmost_rgn = child_stratum.s_regions[0]; 
//                 if(leftmost_rgn.from==0 && leftmost_rgn.size()>=k){
//                     stratum_id=leftmost_rgn.stratum_id;
//                     continue;
//                 }
//             }
//             // stop
//             break;
//         }
//     }
//     return stratum_id;
// }

// void setup_stratum_example_occ(){
//     stratum_occ_samples.resize(prokrustean.stratum_count);
//     vector<bool> visits(prokrustean.stratum_count);
//     std::stack<Stratum> stratum_stack;
//     for(int i=0; i<prokrustean.sequence_count; i++){
//         for(auto &rgn: prokrustean.get_sequence(i).s_regions){
//             auto stratum = prokrustean.get_stratum(rgn.stratum_id);
//             stratum.set_occ(i, rgn.from);
//             stratum_occ_samples[rgn.stratum_id]=stratum.example_occ.value();
//             visits[rgn.stratum_id]=true;

//             stratum_stack.push(stratum);
//             while(!stratum_stack.empty()){
//                 auto stratum=stratum_stack.top();
//                 stratum_stack.pop();
//                 for(auto &c_rgn: stratum.s_regions){
//                     if(visits[c_rgn.stratum_id]) continue;
//                     SeqId seq_id = get<0>(stratum.example_occ.value());
//                     Pos rel_pos = get<1>(stratum.example_occ.value())+c_rgn.from;
//                     auto c_stratum = prokrustean.get_stratum(c_rgn.stratum_id);
//                     c_stratum.set_occ(seq_id, rel_pos);
//                     stratum_occ_samples[c_rgn.stratum_id]=c_stratum.example_occ.value();
//                     visits[c_rgn.stratum_id]=true;
//                     stratum_stack.push(c_stratum);
//                 }
//             }
//         }
//     }
// }
/********************************************************************************************************/
/*                              save/load                                                                */
/********************************************************************************************************/

// Define serialization function for StratifiedData
void _serializeStratifiedData(std::ofstream& file, StratifiedData* data) {
    file.write(reinterpret_cast<const char*>(&data->stratum_id), sizeof(StratumId));
    file.write(reinterpret_cast<const char*>(&data->pos), sizeof(Pos));
}

// Define deserialization function for StratifiedData
void _deserializeStratifiedData(std::ifstream& file, StratifiedData* data) {
    // StratifiedData* data = new StratifiedData();
    file.read(reinterpret_cast<char*>(&data->stratum_id), sizeof(StratumId));
    file.read(reinterpret_cast<char*>(&data->pos), sizeof(Pos));
    // return data;
}

// Serialize the Prokrustean structure
void store_prokrustean(const Prokrustean& data, const std::string& filename) {
    auto start = std::chrono::steady_clock::now();
	cout << "storing prokrustean (" << filename << ") ... " ;
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        // Serialize meta
        file.write(reinterpret_cast<const char*>(&data.version), sizeof(data.version));
        file.write(reinterpret_cast<const char*>(&data.lmin), sizeof(data.lmin));
        file.write(reinterpret_cast<const char*>(&data.sequence_count), sizeof(data.sequence_count));
        file.write(reinterpret_cast<const char*>(&data.total_sequence_region_count), sizeof(data.total_sequence_region_count));
        file.write(reinterpret_cast<const char*>(&data.stratum_count), sizeof(data.stratum_count));
        file.write(reinterpret_cast<const char*>(&data.total_strata_region_count), sizeof(data.total_strata_region_count));

        // Serialize the SequenceSize vector
        for (const auto& size : data.sequences__size) {
            file.write(reinterpret_cast<const char*>(&size), sizeof(SequenceSize));
        }

        // Serialize sequences__region_cnt and StratifiedData vectors
        for (int i = 0; i < data.sequence_count; ++i) {
            uint8_t count = data.sequences__region_cnt[i];
            file.write(reinterpret_cast<const char*>(&count), sizeof(uint8_t));
            
            for (int j = 0; j < count; ++j) {
                _serializeStratifiedData(file, &data.sequences__region[i][j]);
            }
        }

        // Serialize the StratumSize vector
        for (const auto& size : data.stratums__size) {
            file.write(reinterpret_cast<const char*>(&size), sizeof(StratumSize));
        }

        // Serialize stratums__region_cnt and StratifiedData vectors
        for (int i = 0; i < data.stratum_count; ++i) {
            uint8_t count = data.stratums__region_cnt[i];
            file.write(reinterpret_cast<const char*>(&count), sizeof(uint8_t));
            
            for (int j = 0; j < count; ++j) {
                _serializeStratifiedData(file, &data.stratums__region[i][j]);
            }
        }
        
        file.close();
        cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    } else {
        std::cerr << "Unable to open the file for writing." << std::endl;
    }
}

// Deserialize the Prokrustean structure
bool load_prokrustean(const std::string& filename, Prokrustean& data) {
    auto start = std::chrono::steady_clock::now();
	cout << "loading prokrustean (" << filename << ") ... " ;

    std::ifstream file(filename, std::ios::binary);
    
    if (file.is_open()) {
        size_t sequence_count, stratum_count;
        
        // Deserialize meta
        file.read(reinterpret_cast<char*>(&data.version), sizeof(data.version));
        file.read(reinterpret_cast<char*>(&data.lmin), sizeof(data.lmin));
        file.read(reinterpret_cast<char*>(&sequence_count), sizeof(sequence_count));
        file.read(reinterpret_cast<char*>(&data.total_sequence_region_count), sizeof(data.total_sequence_region_count));
        file.read(reinterpret_cast<char*>(&stratum_count), sizeof(stratum_count));
        file.read(reinterpret_cast<char*>(&data.total_strata_region_count), sizeof(data.total_strata_region_count));

        data.set_seq_count(sequence_count);
        data.set_stratum_count(stratum_count);
        // Deserialize the SequenceSize vector
        for (int i = 0; i < sequence_count; ++i) {
            file.read(reinterpret_cast<char*>(&data.sequences__size[i]), sizeof(SequenceSize));
        }
        
        // Deserialize sequences__region_cnt and StratifiedData vectors
        for (int i = 0; i < data.sequence_count; ++i) {
            uint8_t count;
            file.read(reinterpret_cast<char*>(&count), sizeof(uint8_t));
            data.sequences__region_cnt[i]=count;
            data.sequences__region[i]=new StratifiedData[count];
            
            for (int j = 0; j < count; ++j) {
                _deserializeStratifiedData(file, &data.sequences__region[i][j]);
            }
        }
        
        // Deserialize the StratumSize vector
        for (int i = 0; i < stratum_count; ++i) {
            file.read(reinterpret_cast<char*>(&data.stratums__size[i]), sizeof(StratumSize));
        }

        // Deserialize sequences__region_cnt and StratifiedData vectors
        for (int i = 0; i < data.stratum_count; ++i) {
            uint8_t count;
            file.read(reinterpret_cast<char*>(&count), sizeof(uint8_t));
            data.stratums__region_cnt[i]=count;
            data.stratums__region[i]=new StratifiedData[count];
            for (int j = 0; j < count; ++j) {
                _deserializeStratifiedData(file, &data.stratums__region[i][j]);
            }
        }
        
        file.close();
        cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
        return true;
    } else {
        cout << "failed" << endl;
        return false;
    }
}


void debug(ProkrusteanExtension &ext, StratumId target_stratum_id){
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        for(auto& edge: ext.prokrustean.get_sequence(i).s_edges){
            if(edge.stratum_id==target_stratum_id){
                cout << "parent seq: " << i << endl;
                ext.prokrustean.get_sequence(i).print();
            }
        }
    }
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        for(auto& edge: ext.prokrustean.get_stratum(i).s_edges){
            if(edge.stratum_id==target_stratum_id){
                cout << "parent stratum: " << i << endl;
                ext.prokrustean.get_stratum(i).print();
            }
        }
    }
}

void store_prokrustean_text(const Prokrustean& prokrustean, const std::string& filename) {
    auto start = std::chrono::steady_clock::now();
	cout << "storing prokrustean (" << filename << ") ... " ;

    std::ofstream outputFile(filename);

    // Set the width for each column and specify left alignment
    const int columnWidth = 20; 

    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    outputFile << "sequences: " << prokrustean.sequence_count << endl; 
    outputFile << "sequence regions: " << prokrustean.total_sequence_region_count << endl; 
    outputFile << "strata: " << prokrustean.stratum_count << endl; 
    outputFile << "strata regions: " << prokrustean.total_strata_region_count << endl; 
    outputFile << "Lmin: " << prokrustean.lmin << endl;
    outputFile << "** strata are listed after all "<< prokrustean.sequence_count << " sequences are listed " << endl;
    outputFile << std::left << std::setw(columnWidth) << "sequence"
               << std::left << std::setw(columnWidth) << "length"
               << std::left << std::setw(columnWidth) << "stratified regions  [from:to (stratum id)]" << std::endl;
    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    

    for (size_t i = 0; i < prokrustean.sequences__size.size(); ++i) {
        outputFile << std::left << std::setw(columnWidth) << i;
        outputFile << std::left << std::setw(columnWidth) << prokrustean.sequences__size[i];
        for (uint8_t j = 0; j < prokrustean.sequences__region_cnt[i]; ++j) {
            string expr;
            expr+=to_string(prokrustean.sequences__region[i][j].pos);
            expr+=":";
            expr+=to_string(prokrustean.sequences__region[i][j].pos+prokrustean.stratums__size[prokrustean.sequences__region[i][j].stratum_id]);
            expr+=" (";
            expr+=to_string(prokrustean.sequences__region[i][j].stratum_id);
            expr+=")";
            outputFile << std::left << std::setw(columnWidth) <<  expr;
        }
        outputFile<< std::endl;
    }

    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    outputFile << std::left << std::setw(columnWidth) << "stratum"
               << std::left << std::setw(columnWidth) << "length"
               << std::left << std::setw(columnWidth) << "stratified regions  [from:to (stratum id)]" << std::endl;
    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    

    for (size_t i = 0; i < prokrustean.stratums__size.size(); ++i) {
        outputFile << std::left << std::setw(columnWidth) << i;
        outputFile << std::left << std::setw(columnWidth) << prokrustean.stratums__size[i];
        for (uint8_t j = 0; j < prokrustean.stratums__region_cnt[i]; ++j) {
            string expr;
            expr+=to_string(prokrustean.stratums__region[i][j].pos);
            expr+=":";
            expr+=to_string(prokrustean.stratums__region[i][j].pos+prokrustean.stratums__size[prokrustean.stratums__region[i][j].stratum_id]);
            expr+=" (";
            expr+=to_string(prokrustean.stratums__region[i][j].stratum_id);
            expr+=")";
            outputFile << std::left << std::setw(columnWidth) <<  expr;
        }
        outputFile<< std::endl;
    }

    outputFile.close();
    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}



#endif