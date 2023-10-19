#ifndef TEST_STORAGE_HPP_
#define TEST_STORAGE_HPP_
#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <future>
#include "../src/prokrustean.hpp"
#include "../src/prokrustean.support.hpp"

using namespace std;


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
        originalData.set_seq_size(i, 200);
        originalData.set_seq_regions(i, arr, region_per_stratum);
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

    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData;
    load_prokrustean("data.bin", loadedData);

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
}

void test_store_and_retrieve_simple() {
    // Create an instance of the Prokrustean structure
    int region_per_stratum= 2;
    int sequence_cnt=400;
    int stratum_cnt=2262;
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
        originalData.set_seq_size(i, 200);
        originalData.set_seq_regions(i, arr, region_per_stratum);
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

    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData;
    load_prokrustean("data.bin", loadedData);

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
    int sequence_cnt=50;
    int stratum_cnt=100;
    Prokrustean originalData;
    // Populate the data as needed
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    for(int i=0; i<sequence_cnt; i++){
        auto rgn_count = generate_random_val<CoveringRegionIdx>();
        StratifiedData* arr= new StratifiedData[rgn_count];
        for(int r=0; r<rgn_count; r++){
            arr[r].stratum_id=generate_random_val<StratumId>();
            arr[r].pos=generate_random_val<Pos>();
        }
        originalData.set_seq_size(i, generate_random_val<SequenceSize>());
        originalData.set_seq_regions(i, arr, rgn_count);
    }
    for(int i=0; i<stratum_cnt; i++){
        auto rgn_count = generate_random_val<CoveringRegionIdx>();
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

    start = std::chrono::steady_clock::now();

    // Deserialize the data from the file
    Prokrustean loadedData;  
    load_prokrustean("data.bin", loadedData);

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
}

void test_store_and_retrieve_extensions() {
    // Create an instance of the Prokrustean structure
    int region_per_stratum= 2;
    int sequence_cnt=5;
    int stratum_cnt=10;
    // extensions
    Prokrustean originalData;
    originalData.contains_stratum_extension_count=true;
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    for(int i=0; i<stratum_cnt; i++){
        originalData.set_stratum_left_right_extensions(i, 1, 3);
    }
    store_prokrustean(originalData, "data.bin");
    Prokrustean loadedData;
    load_prokrustean("data.bin", loadedData);
    for(int i=0; i<stratum_cnt; i++){
        assert(originalData.get_left_cnt(i)==loadedData.get_left_cnt(i));
        assert(originalData.get_right_cnt(i)==loadedData.get_right_cnt(i));
    }
    // frequencies
    originalData=Prokrustean();
    // Populate the data as needed
    originalData.contains_stratum_frequency=true;
    originalData.set_seq_count(sequence_cnt);
    originalData.set_stratum_count(stratum_cnt);
    for(int i=0; i<stratum_cnt; i++){
        originalData.set_stratum_frequency(i, 5);
    }
    store_prokrustean(originalData, "data.bin");
    loadedData=Prokrustean();
    load_prokrustean("data.bin", loadedData);
    for(int i=0; i<stratum_cnt; i++){
        assert(originalData.stratums__frequency_cnt[i]==loadedData.stratums__frequency_cnt[i]);
    }
}


void main_test_storage(){
    test_store_and_retrieve_simple();
    test_store_and_retrieve_random_data();
    test_store_and_retrieve_extensions();
}
#endif