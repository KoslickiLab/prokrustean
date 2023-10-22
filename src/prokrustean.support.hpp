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

    bool is_sample_occ_available;
    vector<SeqId> stratum_sample_occ_seq_id;
    vector<Pos> stratum_sample_occ_pos;

    vector<CharCount> stratum_left_ext_count;
    vector<CharCount> stratum_right_ext_count;

    vector<SpinLock> stratum_locks;
    int stratum_lock_scale;
    bool stratum_lock_active=false;
    
    ProkrusteanExtension(Prokrustean& prokrustean):prokrustean(prokrustean) {}

    void set_stratum_locks(int thread_cnt){
        this->stratum_lock_scale=thread_cnt*10000;
        this->stratum_locks=vector<SpinLock>(this->stratum_lock_scale);
    }
    void lock_stratum(StratumId id){
        this->stratum_locks[id%this->stratum_lock_scale].lock();
    }
    void unlock_stratum(StratumId id){
        this->stratum_locks[id%this->stratum_lock_scale].unlock();
    }

    bool refracted_at_front_of_cover(int k, uint64_t size, StratifiedData* region, CoveringRegionCount rgn_cnt){
        assert(size>=k);
        if(rgn_cnt==0){
            return true;
        } else if(region[0].pos>0 || prokrustean.get_stratum_size(region[0].stratum_id)<k){
            return true;
        }
        return false;
    }

    bool seq__refracted_at_front_of_cover(SeqId id, int k){
        return refracted_at_front_of_cover(k, prokrustean.sequences__size[id], prokrustean.sequences__region[id], prokrustean.sequences__region_cnt[id]);
    }

    bool stratum__refracted_at_front_of_cover(StratumId id, int k){
        return refracted_at_front_of_cover(k, prokrustean.stratums__size[id], prokrustean.stratums__region[id], prokrustean.stratums__region_cnt[id]);
    }

    bool refracted_at_back_of_cover(int k, uint64_t size, StratifiedData* region, CoveringRegionCount rgn_cnt){
        assert(size>=k);
        if(rgn_cnt==0){
            return true;
        } else{
            auto last_rgn_size = prokrustean.get_stratum_size(region[rgn_cnt-1].stratum_id);
            if(last_rgn_size<k || region[rgn_cnt-1].pos+last_rgn_size<size){
                return true;
            }
        }
        return false;
    }

    bool seq__refracted_at_back_of_cover(SeqId id, int k){
        return refracted_at_back_of_cover(k, prokrustean.sequences__size[id], prokrustean.sequences__region[id], prokrustean.sequences__region_cnt[id]);
    }
    
    bool stratum__refracted_at_back_of_cover(StratumId id, int k){
        return refracted_at_back_of_cover(k, prokrustean.stratums__size[id], prokrustean.stratums__region[id], prokrustean.stratums__region_cnt[id]);
    }

    // bool get_seq_first_stratified(SeqId id, StratifiedEdge &edge){
    //     if(prokrustean.sequences__region_cnt[id]==0){
    //         return false;
    //     } else {
    //         auto &region=prokrustean.sequences__region[id];
    //         edge.update(region[0].pos, region[0].pos + prokrustean.get_stratum_size(region[0].stratum_id), true, region[0].stratum_id);
    //         return true;
    //     }
    // }
    bool get_stratum_first_stratified(StratumId id, StratifiedEdge &edge){
        if(prokrustean.stratums__region_cnt[id]==0){
            return false;
        } else {
            auto &region=prokrustean.stratums__region[id][0];
            edge.update(region.pos, region.pos + prokrustean.get_stratum_size(region.stratum_id), true, region.stratum_id);
            return true;
        }
    }
    bool get_stratum_last_stratified(StratumId id, StratifiedEdge &edge){
        if(prokrustean.stratums__region_cnt[id]==0){
            return false;
        } else {
            auto &region=prokrustean.stratums__region[id][prokrustean.stratums__region_cnt[id]-1];
            edge.update(region.pos, region.pos + prokrustean.get_stratum_size(region.stratum_id), true, region.stratum_id);
            return true;
        }
    }
    optional<StratifiedEdge> first_stratified(StratifiedData* region, CoveringRegionCount rgn_cnt){
        if(rgn_cnt==0){
            return nullopt;
        } else {
            return StratifiedEdge(region[0].pos, region[0].pos + prokrustean.get_stratum_size(region[0].stratum_id), region[0].stratum_id);
        }
    }
    optional<StratifiedEdge> seq__first_stratified(SeqId id){
        return first_stratified(prokrustean.sequences__region[id], prokrustean.sequences__region_cnt[id]);
    }
    optional<StratifiedEdge> stratum__first_stratified(StratumId id){
        return first_stratified(prokrustean.stratums__region[id], prokrustean.stratums__region_cnt[id]);
    }

    optional<StratifiedEdge> last_stratified(StratifiedData* region, CoveringRegionCount rgn_cnt){
        if(rgn_cnt==0){
            return nullopt;
        } else {
            return StratifiedEdge(region[rgn_cnt-1].pos, region[rgn_cnt-1].pos + prokrustean.get_stratum_size(region[rgn_cnt-1].stratum_id), region[rgn_cnt-1].stratum_id);
        }
    }
    optional<StratifiedEdge> seq__last_stratified(SeqId id){
        return last_stratified(prokrustean.sequences__region[id], prokrustean.sequences__region_cnt[id]);
    }
    optional<StratifiedEdge> stratum__last_stratified(StratumId id){
        return last_stratified(prokrustean.stratums__region[id], prokrustean.stratums__region_cnt[id]);
    }

    void increment_left_ext(StratumId id){
        if(stratum_lock_active){
            lock_stratum(id);
            stratum_left_ext_count[id]++;
            unlock_stratum(id);
        } else {
            stratum_left_ext_count[id]++;
        }
    }
    void increment_right_ext(StratumId id){
        if(stratum_lock_active){
            lock_stratum(id);
            stratum_right_ext_count[id]++;
            unlock_stratum(id);
        } else {
            stratum_right_ext_count[id]++;
        }
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

/********************************************************************************************************/
/*                              Setting occurrence information of strata                                */
/********************************************************************************************************/
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
    ext.is_sample_occ_available=true;
}

void _set_stratum_example_occ_for_seq(SeqId id, ProkrusteanExtension &ext, vector<uint8_t> &visits){
    std::stack<StratumId> stratum_stack;
    for(auto &rgn: ext.prokrustean.get_sequence(id).s_edges){
        ext.lock_stratum(rgn.stratum_id);
        if(visits[rgn.stratum_id]==1){
            ext.unlock_stratum(rgn.stratum_id);
            continue;
        } else {
            visits[rgn.stratum_id]=1;
            ext.unlock_stratum(rgn.stratum_id);
            ext.stratum_sample_occ_seq_id[rgn.stratum_id]=id;
            ext.stratum_sample_occ_pos[rgn.stratum_id]=rgn.from;
            stratum_stack.push(rgn.stratum_id);
        }

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
                    ext.unlock_stratum(c_rgn.stratum_id);
                    SeqId seq_id=ext.stratum_sample_occ_seq_id[stratum_id];
                    Pos rel_pos=ext.stratum_sample_occ_pos[stratum_id]+c_rgn.from;
                    ext.stratum_sample_occ_seq_id[c_rgn.stratum_id]=seq_id;
                    ext.stratum_sample_occ_pos[c_rgn.stratum_id]=rel_pos;
                    stratum_stack.push(c_rgn.stratum_id);
                }
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
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(std::async(std::launch::async, func_, ref(visits), ref(ext), ref(seq_idx_gen)));}
    for (auto &f : futures) {f.wait();}
    ext.is_sample_occ_available=true;
}

/********************************************************************************************************/
/*                              Calculating left/right extension characters 
                                 cannot work because actual characters are required - expensive          */
/********************************************************************************************************/
void _increment_leftright_extension_for_non_intersected_descendents(Vertex &vertex, ProkrusteanExtension &ext, stack<StratumId> &stratum_stack){
    for(int e=0; e<vertex.s_edges.size(); e++){
        auto &edge=vertex.s_edges[e];
        if(edge.from>0){
            auto limited_length=0;
            if(0<e && edge.from<vertex.s_edges[e-1].to){
                // intersection exists
                limited_length=vertex.s_edges[e-1].to-edge.from;
            }
            stratum_stack.push(edge.stratum_id);
            while(!stratum_stack.empty()){
                auto stratum_id=stratum_stack.top();
                stratum_stack.pop();
                // ext.prokrustean.get_left_cnt(stratum_id)++;
                ext.increment_left_ext(stratum_id);
                // ext.stratum_left_ext_count[stratum_id]++;
                if(0<ext.prokrustean.stratums__region_cnt[stratum_id]
                && ext.prokrustean.stratums__region[stratum_id][0].pos==0) {
                    auto c_stratum_id=ext.prokrustean.stratums__region[stratum_id][0].stratum_id;
                    if(limited_length<ext.prokrustean.get_stratum_size(c_stratum_id)){
                        stratum_stack.push(c_stratum_id);
                    }
                }
            }
        }
        if(edge.to<vertex.size){
            auto limited_length=0;
            if(e+1<vertex.s_edges.size() && vertex.s_edges[e+1].from<edge.to){
                // intersection exists
                limited_length=edge.to-vertex.s_edges[e+1].from;
            }
            stratum_stack.push(edge.stratum_id);
            while(!stratum_stack.empty()){
                auto stratum_id=stratum_stack.top();
                stratum_stack.pop();
                // ext.stratum_right_ext_count[stratum_id]++;
                ext.increment_right_ext(stratum_id);
                // optional<StratifiedEdge> last_edge= ext.stratum__last_stratified(stratum_id);
                // if(last_edge.has_value()
                // &&last_edge.value().to==ext.prokrustean.get_stratum_size(stratum_id)
                // &&limited_length<last_edge.value().size()
                // ){
                //     stratum_stack.push(last_edge.value().stratum_id);
                // }
                if(0<ext.prokrustean.stratums__region_cnt[stratum_id]){
                    auto &rgn = ext.prokrustean.stratums__region[stratum_id][ext.prokrustean.stratums__region_cnt[stratum_id]-1];
                    if(rgn.pos + ext.prokrustean.get_stratum_size(rgn.stratum_id)==ext.prokrustean.get_stratum_size(stratum_id)
                    && limited_length < ext.prokrustean.get_stratum_size(rgn.stratum_id)){
                        stratum_stack.push(rgn.stratum_id);
                    }
                }
            }
        }       
    }
}
void count_left_right_character_extensions(ProkrusteanExtension &ext){
    ext.stratum_left_ext_count.resize(ext.prokrustean.stratum_count,0);
    ext.stratum_right_ext_count.resize(ext.prokrustean.stratum_count,0);
    std::stack<StratumId> stratum_stack;
    Vertex vertex;
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        ext.prokrustean.get_sequence(i, vertex);
        _increment_leftright_extension_for_non_intersected_descendents(vertex, ext, stratum_stack);
    }
    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        ext.prokrustean.get_stratum(i, vertex);
        _increment_leftright_extension_for_non_intersected_descendents(vertex, ext, stratum_stack);
    }
}
void count_left_right_character_extensions_parallel(ProkrusteanExtension &ext, int thread_cnt){
    ext.set_stratum_locks(thread_cnt);
    ext.stratum_lock_active=true;
    ext.stratum_left_ext_count.resize(ext.prokrustean.stratum_count,0);
    ext.stratum_right_ext_count.resize(ext.prokrustean.stratum_count,0);
    vector<future<void>> futures;
    auto func_ = [](ProkrusteanExtension &ext, int thread_idx, int thread_cnt) {
        stack<StratumId> stratum_stack;
        Vertex vertex;
        for(int i=thread_idx; i<ext.prokrustean.sequence_count; i+=thread_cnt){
            ext.prokrustean.get_sequence(i, vertex);
            _increment_leftright_extension_for_non_intersected_descendents(vertex, ext, stratum_stack);
        }
        for(int i=thread_idx; i<ext.prokrustean.stratum_count; i+=thread_cnt){
            ext.prokrustean.get_stratum(i, vertex);
            _increment_leftright_extension_for_non_intersected_descendents(vertex, ext, stratum_stack);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(std::async(std::launch::async, func_, ref(ext), i, thread_cnt));}
    for (auto &f : futures) {f.wait();}
}

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
        file.write(reinterpret_cast<const char*>(&data.contains_stratum_extension_count), sizeof(data.contains_stratum_extension_count));
        file.write(reinterpret_cast<const char*>(&data.contains_stratum_frequency), sizeof(data.contains_stratum_frequency));

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
            CoveringRegionCount count = data.stratums__region_cnt[i];
            file.write(reinterpret_cast<const char*>(&count), sizeof(CoveringRegionCount));
            
            for (int j = 0; j < count; ++j) {
                _serializeStratifiedData(file, &data.stratums__region[i][j]);
            }
        }
        if(data.contains_stratum_extension_count){
            for (int i = 0; i < data.stratum_count; ++i) {
                file.write(reinterpret_cast<const char*>(&data.stratums__character_extensions_cnt[i]), sizeof(CharCount));
            }
        }

        if(data.contains_stratum_frequency){
            for (int i = 0; i < data.stratum_count; ++i) {
                file.write(reinterpret_cast<const char*>(&data.stratums__frequency_cnt[i]), sizeof(FrequencyCount));
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
	cout << "loading prokrustean (" << filename << ") ... " << endl;

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
        file.read(reinterpret_cast<char*>(&data.contains_stratum_extension_count), sizeof(data.contains_stratum_extension_count));
        file.read(reinterpret_cast<char*>(&data.contains_stratum_frequency), sizeof(data.contains_stratum_frequency));

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

        if(data.contains_stratum_extension_count){
            for (int i = 0; i < stratum_count; ++i) {
                file.read(reinterpret_cast<char*>(&data.stratums__character_extensions_cnt[i]), sizeof(CharCount));
            }
        }

        if(data.contains_stratum_frequency){
            for (int i = 0; i < stratum_count; ++i) {
                file.read(reinterpret_cast<char*>(&data.stratums__frequency_cnt[i]), sizeof(FrequencyCount));
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

void store_prokrustean_text(Prokrustean& prokrustean, const std::string& filename) {
    auto start = std::chrono::steady_clock::now();
	cout << "storing prokrustean (" << filename << ") ... " ;

    std::ofstream outputFile(filename);

    // Set the width for each column and specify left alignment
    const int columnWidth = 15; 
    const int columnWidthMedium = 20; 
    const int columnWidthLarge = 25; 

    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    outputFile << "sequences: " << prokrustean.sequence_count << endl; 
    outputFile << "sequence regions: " << prokrustean.total_sequence_region_count << endl; 
    outputFile << "strata: " << prokrustean.stratum_count << endl; 
    outputFile << "strata regions: " << prokrustean.total_strata_region_count << endl; 
    outputFile << "Lmin: " << prokrustean.lmin << endl;
    outputFile << "** sequences are listed after all "<< prokrustean.stratum_count << " strata are listed " << endl;
    outputFile << std::left << std::setw(columnWidth) << "strata"
               << std::left << std::setw(columnWidth) << "length";
    if(prokrustean.contains_stratum_extension_count){
        outputFile << std::left << std::setw(columnWidthLarge) << "char exts [left|right]";
    }
    if(prokrustean.contains_stratum_frequency){
        outputFile << std::left << std::setw(columnWidth) << "frequency";
    }
    outputFile << std::left << std::setw(columnWidthMedium) << "stratified regions  [from:to (stratum id)]" << std::endl;
    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    
    for (size_t i = 0; i < prokrustean.stratums__size.size(); ++i) {
        outputFile << std::left << std::setw(columnWidth) << i;
        outputFile << std::left << std::setw(columnWidth) << prokrustean.stratums__size[i];
        if(prokrustean.contains_stratum_extension_count){
            string expr;
            expr+=to_string(prokrustean.get_left_cnt(i));
            expr+=" | ";
            expr+=to_string(prokrustean.get_right_cnt(i));
            outputFile << std::left << std::setw(columnWidthLarge) << expr;
        }
        if(prokrustean.contains_stratum_frequency){
            outputFile << std::left << std::setw(columnWidth) << prokrustean.stratums__frequency_cnt[i];
        }
        for (uint8_t j = 0; j < prokrustean.stratums__region_cnt[i]; ++j) {
            string expr;
            expr+=to_string(prokrustean.stratums__region[i][j].pos);
            expr+=":";
            expr+=to_string(prokrustean.stratums__region[i][j].pos+prokrustean.stratums__size[prokrustean.stratums__region[i][j].stratum_id]);
            expr+=" (";
            expr+=to_string(prokrustean.stratums__region[i][j].stratum_id);
            expr+=")";
            outputFile << std::left << std::setw(columnWidthMedium) <<  expr;
        }
        outputFile<< std::endl;
    }

    outputFile << "------------------------------------------------------------------------------------------------" << endl;
    outputFile << std::left << std::setw(columnWidth) << "sequences"
               << std::left << std::setw(columnWidth) << "length";
    outputFile << std::left << std::setw(columnWidthMedium) << "stratified regions  [from:to (stratum id)]" << std::endl;
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
            outputFile << std::left << std::setw(columnWidthMedium) <<  expr;
        }
        outputFile<< std::endl;
    }

    outputFile.close();
    cout << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
}



#endif