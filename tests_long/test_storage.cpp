#ifndef TEST_STORAGE_HPP_
#define TEST_STORAGE_HPP_
#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "../src/prokrustean.hpp"

using namespace std;


// Define serialization function for StratifiedData
void serializeStratifiedData(std::ofstream& file, StratifiedData* data) {
    file.write(reinterpret_cast<const char*>(&data->stratum_id), sizeof(StratumId));
    file.write(reinterpret_cast<const char*>(&data->pos), sizeof(Pos));
}

// Define deserialization function for StratifiedData
void deserializeStratifiedData(std::ifstream& file, StratifiedData* data) {
    // StratifiedData* data = new StratifiedData();
    file.read(reinterpret_cast<char*>(&data->stratum_id), sizeof(StratumId));
    file.read(reinterpret_cast<char*>(&data->pos), sizeof(Pos));
    // return data;
}

// Serialize the Prokrustean structure
void store_prokrustean(const Prokrustean& data, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        // Serialize the sizes
        file.write(reinterpret_cast<const char*>(&data.sequence_count), sizeof(data.sequence_count));
        file.write(reinterpret_cast<const char*>(&data.stratum_count), sizeof(data.stratum_count));

        // Serialize the SequenceSize vector
        for (const auto& size : data.sequences__size) {
            file.write(reinterpret_cast<const char*>(&size), sizeof(SequenceSize));
        }

        // Serialize sequences__region_cnt and StratifiedData vectors
        for (int i = 0; i < data.sequence_count; ++i) {
            uint8_t count = data.sequences__region_cnt[i];
            file.write(reinterpret_cast<const char*>(&count), sizeof(uint8_t));
            
            if(i==0){
                cout << "i count save " << (int)count << endl;
            }
            for (int j = 0; j < count; ++j) {
                serializeStratifiedData(file, &data.sequences__region[i][j]);
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
                serializeStratifiedData(file, &data.stratums__region[i][j]);
            }
        }
        
        file.close();
    } else {
        std::cerr << "Unable to open the file for writing." << std::endl;
    }
}

// Deserialize the Prokrustean structure
Prokrustean load_prokrustean(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    Prokrustean data;
    
    if (file.is_open()) {
        size_t sequence_count, stratum_count;
        
        // Deserialize the sizes
        file.read(reinterpret_cast<char*>(&sequence_count), sizeof(sequence_count));
        file.read(reinterpret_cast<char*>(&stratum_count), sizeof(stratum_count));
        
        data.set_seq_count(sequence_count);
        data.set_stratum_count(stratum_count);
        cout << "data.sequence_count " << data.sequence_count << " data.stratum_count " << data.stratum_count << endl;
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
                deserializeStratifiedData(file, &data.sequences__region[i][j]);
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
                deserializeStratifiedData(file, &data.stratums__region[i][j]);
            }
        }
        
        
        
        file.close();
    } else {
        std::cerr << "Unable to open the file for reading." << std::endl;
    }
    
    return data;
}

void test_store_and_retrieve() {
    // Create an instance of the Prokrustean structure
    int region_per_stratum= 2;
    int sequence_cnt=400000;
    int stratum_cnt=2262035;
    Prokrustean originalData;
    // Populate the data as needed
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    originalData.sequences__size[5]=3;
    for(int i=0; i<sequence_cnt; i++){
        StratifiedData* arr= new StratifiedData[region_per_stratum];
        for(int r=0; r<region_per_stratum; r++){
            arr[r].stratum_id=i;
            arr[r].pos=6;    
        }
        originalData.set_seq_regions(i, 200, arr, region_per_stratum);
    }
    for(int i=0; i<stratum_cnt; i++){
        StratifiedData* arr= new StratifiedData[region_per_stratum];
        for(int r=0; r<region_per_stratum; r++){
            arr[r].stratum_id=i;
            arr[r].pos=5;    
        }
        originalData.set_stratum_regions(i, arr, region_per_stratum);
    }
    auto start = std::chrono::steady_clock::now();

    // Serialize the data to a file
    store_prokrustean(originalData, "data.bin");

    cout << "save: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData = load_prokrustean("data.bin");

    cout << "load: "  << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // Now, loadedData contains the deserialized data
    for(int i=0; i<sequence_cnt; i++){
        assert(originalData.sequences__size[i]==loadedData.sequences__size[i]);
        cout << "i: " << i <<" originalData.sequences__region_cnt[i]: " << (int)originalData.sequences__region_cnt[i] << " loadedData.sequences__region_cnt[i] " << (int)loadedData.sequences__region_cnt[i] <<endl; 
        assert(originalData.sequences__region_cnt[i]==loadedData.sequences__region_cnt[i]);
        for(int r=0; r<loadedData.sequences__region_cnt[i]; r++){
            assert(originalData.sequences__region[i]->stratum_id==loadedData.sequences__region[i]->stratum_id);
            assert(originalData.sequences__region[i]->pos==loadedData.sequences__region[i]->pos);
        }
    }
    for(int i=0; i<stratum_cnt; i++){
        assert(originalData.stratums__size[i]==loadedData.stratums__size[i]);
        assert(originalData.stratums__region_cnt[i]==loadedData.stratums__region_cnt[i]);
        for(int r=0; r<loadedData.stratums__region_cnt[i]; r++){
            assert(originalData.stratums__region[i]->stratum_id==loadedData.stratums__region[i]->stratum_id);
            assert(originalData.stratums__region[i]->pos==loadedData.stratums__region[i]->pos);
        }
    }
    cout << "done " << endl;
}

void test_store_and_retrieve_simple() {
    // Create an instance of the Prokrustean structure
    int region_per_stratum= 2;
    int sequence_cnt=40000;
    int stratum_cnt=226203;
    Prokrustean originalData;
    // Populate the data as needed
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    originalData.sequences__size[5]=3;
    for(int i=0; i<sequence_cnt; i++){
        StratifiedData* arr= new StratifiedData[region_per_stratum];
        for(int r=0; r<region_per_stratum; r++){
            arr[r].stratum_id=i;
            arr[r].pos=6;    
        }
        originalData.set_seq_regions(i, 200, arr, region_per_stratum);
    }
    for(int i=0; i<stratum_cnt; i++){
        StratifiedData* arr= new StratifiedData[region_per_stratum];
        for(int r=0; r<region_per_stratum; r++){
            arr[r].stratum_id=i;
            arr[r].pos=5;    
        }
        originalData.set_stratum_regions(i, arr, region_per_stratum);
    }
    auto start = std::chrono::steady_clock::now();

    // Serialize the data to a file
    store_prokrustean(originalData, "data.bin");

    cout << "save: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData = load_prokrustean("data.bin");

    cout << "load: "  << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // Now, loadedData contains the deserialized data
    for(int i=0; i<sequence_cnt; i++){
        assert(originalData.sequences__size[i]==loadedData.sequences__size[i]);
        assert(originalData.sequences__region_cnt[i]==loadedData.sequences__region_cnt[i]);
        for(int r=0; r<loadedData.sequences__region_cnt[i]; r++){
            assert(originalData.sequences__region[i]->stratum_id==loadedData.sequences__region[i]->stratum_id);
            assert(originalData.sequences__region[i]->pos==loadedData.sequences__region[i]->pos);
        }
    }
    for(int i=0; i<stratum_cnt; i++){
        assert(originalData.stratums__size[i]==loadedData.stratums__size[i]);
        assert(originalData.stratums__region_cnt[i]==loadedData.stratums__region_cnt[i]);
        for(int r=0; r<loadedData.stratums__region_cnt[i]; r++){
            assert(originalData.stratums__region[i]->stratum_id==loadedData.stratums__region[i]->stratum_id);
            assert(originalData.stratums__region[i]->pos==loadedData.stratums__region[i]->pos);
        }
    }
    cout << "done " << endl;
}


template <typename T>
T generate_random_val() {
    // Ensure T is an integer type
    static_assert(std::is_integral<T>::value, "Type must be an integral type");

    std::random_device rd; // Obtain a seed from the system entropy device
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_int_distribution<T> dist(0, std::numeric_limits<T>::max());

    return dist(gen);
}


void test_store_and_retrieve_random_data() {
    // Create an instance of the Prokrustean structure
    int sequence_cnt=40000;
    int stratum_cnt=800000;
    Prokrustean originalData;
    // Populate the data as needed
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    for(int i=0; i<sequence_cnt; i++){
        auto rgn_count = generate_random_val<StratifyingRegionIdx>();
        StratifiedData* arr= new StratifiedData[rgn_count];
        for(int r=0; r<rgn_count; r++){
            arr[r].stratum_id=generate_random_val<StratumId>();
            arr[r].pos=generate_random_val<Pos>();
        }
        originalData.set_seq_regions(i, generate_random_val<SequenceSize>(), arr, rgn_count);
    }
    for(int i=0; i<stratum_cnt; i++){
        auto rgn_count = generate_random_val<StratifyingRegionIdx>();
        StratifiedData* arr= new StratifiedData[rgn_count];
        for(int r=0; r<rgn_count; r++){
            arr[r].stratum_id=generate_random_val<StratumId>();
            arr[r].pos=generate_random_val<Pos>();
        }
        originalData.set_stratum_regions(i, arr, rgn_count);
        originalData.stratums__size[i]=generate_random_val<StratumSize>();
    }
    auto start = std::chrono::steady_clock::now();

    // Serialize the data to a file
    store_prokrustean(originalData, "data.bin");

    cout << "save: " << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData = load_prokrustean("data.bin");

    cout << "load: "  << (std::chrono::steady_clock::now()-start).count()/1000000 << "ms" << endl;
    // Now, loadedData contains the deserialized data
    for(int i=0; i<sequence_cnt; i++){
        assert(originalData.sequences__size[i]==loadedData.sequences__size[i]);
        assert(originalData.sequences__region_cnt[i]==loadedData.sequences__region_cnt[i]);
        for(int r=0; r<loadedData.sequences__region_cnt[i]; r++){
            assert(originalData.sequences__region[i]->stratum_id==loadedData.sequences__region[i]->stratum_id);
            assert(originalData.sequences__region[i]->pos==loadedData.sequences__region[i]->pos);
        }
    }
    for(int i=0; i<stratum_cnt; i++){
        assert(originalData.stratums__size[i]==loadedData.stratums__size[i]);
        assert(originalData.stratums__region_cnt[i]==loadedData.stratums__region_cnt[i]);
        for(int r=0; r<loadedData.stratums__region_cnt[i]; r++){
            assert(originalData.stratums__region[i]->stratum_id==loadedData.stratums__region[i]->stratum_id);
            assert(originalData.stratums__region[i]->pos==loadedData.stratums__region[i]->pos);
        }
    }
    cout << "done " << endl;
}


void main_test_storage(){
    test_store_and_retrieve_simple();
    test_store_and_retrieve_random_data();
}
#endif