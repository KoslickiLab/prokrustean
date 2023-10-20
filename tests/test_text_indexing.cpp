#ifndef TEST_TEXT_INDEXING_HPP_
#define TEST_TEXT_INDEXING_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <iostream>	
#include "const.cpp"
#include "../src/fm_index/index.hpp"
#include "../src/fm_index/string.sdsl.hpp"
#include "../src/util/string.access.hpp"

using namespace std;

void test_indexing_simple(){
    std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
    vector<int> seq_sizes;
    for(auto &str: strings){
        seq_sizes.push_back(str.size());
    }
    DiskSequenceAccess sequnce_access("data.dat");
    sequnce_access.write_open();
    sequnce_access.write_meta(seq_sizes);
    sequnce_access.write_strings(strings);
    sequnce_access.write_close();
    // sequnce_access.load_all_strings();
    sequnce_access.read_open();

    cout << sequnce_access.get_string(1) << endl;
    assert("orange"==sequnce_access.get_string(1));
    assert("rang"==sequnce_access.get_substring(1, 1, 4));
    sequnce_access.read_close();
}

void test_bwt_indexing(){
    WaveletString str(PATH6_CDBG_SAMPLE2, '$');
    auto fm_idx = FmIndex(str);
    
    vector<string> seq_texts;
    fm_idx.recover_all_texts(seq_texts);
    vector<tuple<int, int, int>> substring_locs={
        make_tuple(1, 3, 20),
        make_tuple(3, 8, 37),
        make_tuple(9, 20, 42)
    };
    vector<int> seq_sizes;
    for(auto &str: seq_texts){
        seq_sizes.push_back(str.size());
    }
    DiskSequenceAccess sequnce_access("data1.dat");
    sequnce_access.write_open();
    sequnce_access.write_meta(seq_sizes);
    sequnce_access.write_strings(seq_texts);
    for(auto &str: seq_texts){
        sequnce_access.write_single(str);
    }
    sequnce_access.write_close(); 
    sequnce_access.read_open();
    for(auto [idx, from, to]: substring_locs){
        assert(seq_texts[idx]==sequnce_access.get_string(idx));
        assert(seq_texts[idx].substr(from, to-from+1)==sequnce_access.get_substring(idx, from, to-from+1));
    }
    sequnce_access.read_close(); 
    //in memory
    sequnce_access.load_all_strings();
    for(auto [idx, from, to]: substring_locs){
        assert(seq_texts[idx]==sequnce_access.get_string(idx));
        assert(seq_texts[idx].substr(from, to-from+1)==sequnce_access.get_substring(idx, from, to-from+1));
    }
}

void test_verification(){
    std::vector<std::string> strings = {"apple", "orange", "banana", "grape"};
    DiskSequenceAccess sequnce_access("data.dat");
    vector<int> seq_sizes;
    for(auto &str: strings){
        seq_sizes.push_back(str.size());
    }
    sequnce_access.write_open();
    sequnce_access.write_meta(seq_sizes);
    sequnce_access.write_strings(strings);
    sequnce_access.write_close();
    DiskSequenceAccess sequnce_access1("data.dat");
    sequnce_access1.verify();
}

// class Access{
// public:
//     std::string filename;
//     std::ofstream outfile;
//     std::ifstream infile;
//     int METADATA_SIZE=256;
//     string EVIDENCE="aaa";
//     void save_meta(const std::vector<int>& string_sizes) {
//         // Ensure the metadata string is the expected length
//         this->outfile.write(EVIDENCE.c_str(), METADATA_SIZE);

//         // Compute and write the index table
//         int position = METADATA_SIZE + sizeof(int) * string_sizes.size(); // meta + index table size + all strings size
//         for (int size : string_sizes) {
//             this->outfile.write(reinterpret_cast<char*>(&position), sizeof(position));
//             position += size + sizeof(int); // move by string size + size of the integer that stores the length
//         }
//     }
//     void save_strings(const std::vector<std::string>& strings) {
//         // Write the strings
//         for (const std::string& str : strings) {
//             int len = str.size();
//             this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));  // this could be optional now, as we have sizes saved
//             this->outfile.write(str.c_str(), len);
//         }
//         this->outfile.close();
//     }


//     std::string get_string(int index) {
//         // Open the file for reading
//         std::ifstream infile(this->filename, std::ios::binary);

//         // Skip metadata
//         infile.seekg(METADATA_SIZE);

//         // Calculate the position in the index table for the given index
//         int index_position = METADATA_SIZE + index * sizeof(int);

//         // Seek to that position in the index table
//         infile.seekg(index_position);

//         // Read the position of the string
//         int string_position;
//         infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

//         // Seek to the position of the string
//         infile.seekg(string_position);

//         // Read the length of the string (assuming we stored lengths before the strings, as in the earlier code)
//         int len;
//         infile.read(reinterpret_cast<char*>(&len), sizeof(len));

//         // Read the string itself
//         char* buffer = new char[len + 1];
//         infile.read(buffer, len);
//         buffer[len] = '\0';

//         std::string result(buffer);
//         delete[] buffer;

//         return result;
//     }

//     std::vector<std::string> get_strings(int from, int to) {
//         std::vector<std::string> results;

//         // Open the file for reading
//         std::ifstream infile(this->filename, std::ios::binary);
        
//         // Skip metadata
//         infile.seekg(METADATA_SIZE);

//         // Calculate the position in the index table for the 'from' index
//         int from_index_position = METADATA_SIZE + from * sizeof(int);

//         // Seek to that position in the index table
//         infile.seekg(from_index_position);

//         // Read the position of the 'from' string
//         int from_string_position;
//         infile.read(reinterpret_cast<char*>(&from_string_position), sizeof(from_string_position));

//         // Calculate the position in the index table for the 'to' index
//         int to_index_position = METADATA_SIZE + to * sizeof(int);

//         // Seek to that position in the index table
//         infile.seekg(to_index_position);

//         // Read the position of the 'to' string
//         int to_string_position;
//         infile.read(reinterpret_cast<char*>(&to_string_position), sizeof(to_string_position));

//         // Now, calculate the total length of the chunk of strings we need to read
//         int chunk_length = to_string_position - from_string_position;

//         // Read the entire chunk of strings
//         char* buffer = new char[chunk_length];
//         infile.seekg(from_string_position);
//         infile.read(buffer, chunk_length);

//         // Split the chunk into individual strings based on their lengths and store them in the results vector
//         char* ptr = buffer;
//         char* end = buffer + chunk_length;
//         while (ptr < end) {
//             int len;
//             memcpy(&len, ptr, sizeof(len));
//             ptr += sizeof(len);

//             std::string str(ptr, ptr + len);
//             results.push_back(str);
//             ptr += len;
//         }

//         delete[] buffer;
//         return results;
//     }

//     std::string get_substring(int index, int pos, int size) {
//         // Open the file for reading
//         std::ifstream infile(this->filename, std::ios::binary);

//         // Skip metadata
//         infile.seekg(METADATA_SIZE);

//         // Calculate the position in the index table for the given index
//         int index_position = METADATA_SIZE + index * sizeof(int);

//         // Seek to that position in the index table
//         infile.seekg(index_position);

//         // Read the position of the string
//         int string_position;
//         infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

//         // Seek to the position of the string's length
//         infile.seekg(string_position);

//         // Read the length of the string
//         int len;
//         infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        
//         // Check if the requested substring is valid
//         if (pos < 0 || pos + size > len) {
//             throw std::out_of_range("Requested substring is out of range");
//         }

//         // Calculate the position of the substring within the string and seek to it
//         int substring_position = string_position + sizeof(int) + pos; // +sizeof(int) to account for the string's length
//         infile.seekg(substring_position);

//         // Read the substring
//         char* buffer = new char[size + 1];
//         infile.read(buffer, size);
//         buffer[size] = '\0';

//         std::string result(buffer);
//         delete[] buffer;

//         return result;
//     }


// };

// void test_bwt_indexing2(){
//     WaveletString str(PATH6_CDBG_SAMPLE2, '$');
//     auto fm_idx = FmIndex(str);
    
//     vector<string> seq_texts;
//     fm_idx.recover_all_texts(seq_texts);
//     // seq_texts={ "orange", "apple"};
//     vector<tuple<int, int, int>> substring_locs={
//         make_tuple(1, 3, 20),
//         make_tuple(3, 8, 37),
//         make_tuple(9, 20, 42)
//     };
//     // vector<tuple<int, int, int>> substring_locs={
//     //     make_tuple(1, 1, 3),
//     //     make_tuple(0, 1, 3)
//     // };
//     Access sequnce_access;
//     sequnce_access.filename="test_data.dat";
//     sequnce_access.outfile=std::ofstream(sequnce_access.filename, std::ios::binary);
//     vector<int> sizes;
//     for(auto &str: seq_texts){
//         sizes.push_back(str.size());
//     }
//     sequnce_access.save_meta(sizes);
//     sequnce_access.save_strings(seq_texts);
//     for(auto [idx, from, size]: substring_locs){
//         cout << sequnce_access.get_string(idx) << endl;
//         assert(seq_texts[idx]==sequnce_access.get_string(idx));
//         assert(seq_texts[idx].substr(from, size)==sequnce_access.get_substring(idx, from, size));
//     }
// }


void main_text_indexing() {
    // Call all tests. Using a test framework would simplify this.
    test_indexing_simple();
    test_bwt_indexing();
    test_verification();
}

#endif