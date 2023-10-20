#ifndef UTIL_STRING_ACCESS_HPP_
#define UTIL_STRING_ACCESS_HPP_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


/* for both disk-based, memory-based */
class AbstractSequenceAccess{
public:
	virtual std::string get_string(int index)=0;
    virtual std::string get_substring(int index, uint64_t i, uint64_t size)=0;
};

std::string EVIDENCE = "I am Prokrustes. Do not confuse me with the greek creature - Procrustes. I hate him.";

class DiskSequenceAccess: public AbstractSequenceAccess{
public:
    std::string filename;
    std::ofstream outfile;
    std::ifstream infile;
    bool loaded=false;
    vector<string> strings;
    DiskSequenceAccess(std::string filename): filename(filename){
    }

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    std::string get_metadata() {
        std::ifstream tempfile(filename, std::ios::binary);
        char metadata[METADATA_SIZE + 1];
        tempfile.read(metadata, METADATA_SIZE);
        metadata[METADATA_SIZE] = '\0';
        return std::string(metadata);
    }

    void verify(){
        if (EVIDENCE != get_metadata()){
            throw std::invalid_argument("the file does not exist or seems like not the sequence index generated with the prokrustean. "+filename);
        }
    }

    void write_meta(const std::vector<int>& string_sizes) {
        // Ensure the metadata string is the expected length
        this->outfile.write(EVIDENCE.c_str(), METADATA_SIZE);

        // Compute and write the index table
        int position = METADATA_SIZE + sizeof(int) * string_sizes.size(); // meta + index table size + all strings size
        for (int size : string_sizes) {
            this->outfile.write(reinterpret_cast<char*>(&position), sizeof(position));
            position += size + sizeof(int); // move by string size + size of the integer that stores the length
        }
        
        // Pre-allocate space for all the strings
        outfile.seekp(position - 1);  // Move to the end of the pre-allocated space
        char dummy_byte = 0;
        outfile.write(&dummy_byte, 1);  // Write a single byte to make the file system allocate space up to this position
    }

    void write_strings(const std::vector<std::string>& strings) {
        // Write the strings
        for (const std::string& str : strings) {
            int len = str.size();
            this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));  // this could be optional now, as we have sizes saved
            this->outfile.write(str.c_str(), len);
        }
    }

    void write_open() {
        this->outfile=std::ofstream(filename, std::ios::binary);
    }

    void write_single(const std::string& str){
        int len = str.size();
        this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));
        this->outfile.write(str.c_str(), len);
    }

    // void save_string(int index, const std::string& str) {
    //     // Open the file for reading and writing
    //     std::fstream file(this->filename, std::ios::binary | std::ios::in | std::ios::out);

    //     // Skip metadata
    //     file.seekg(METADATA_SIZE);

    //     // Calculate the position in the index table for the given index
    //     int index_position = METADATA_SIZE + index * sizeof(int);

    //     // Seek to that position in the index table
    //     file.seekg(index_position);

    //     // Read the position of the string
    //     int string_position;
    //     file.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

    //     // Seek to the position of the string in the file
    //     file.seekg(string_position);

    //     // Write the string's length and content to the file
    //     int len = str.size();
    //     file.write(reinterpret_cast<char*>(&len), sizeof(len));
    //     file.write(str.c_str(), len);

    //     file.close();
    // }   

    void write_close() {
        this->outfile.close();
    }

    void read_open() {
        this->infile=std::ifstream(filename, std::ios::binary);
    }

    void read_close() {
        this->infile.close();
    }

    void load_all_strings() {
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(METADATA_SIZE);  // Skip the metadata
        
        while (!infile.eof()) {
            int len = 0;
            infile.read(reinterpret_cast<char*>(&len), sizeof(len));
            if (infile.eof()) break;  // Exit loop if end of file reached

            char* buffer = new char[len + 1];
            infile.read(buffer, len);
            buffer[len] = '\0';
            
            this->strings.push_back(std::string(buffer));
            delete[] buffer;
        }
        this->loaded=true;
    }

    std::string get_string(int index) {
        // Open the file for reading
        std::ifstream infile(this->filename, std::ios::binary);

        // Skip metadata
        infile.seekg(METADATA_SIZE);

        // Calculate the position in the index table for the given index
        int index_position = METADATA_SIZE + index * sizeof(int);

        // Seek to that position in the index table
        infile.seekg(index_position);

        // Read the position of the string
        int string_position;
        infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

        // Seek to the position of the string
        infile.seekg(string_position);

        // Read the length of the string (assuming we stored lengths before the strings, as in the earlier code)
        int len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));

        // Read the string itself
        char* buffer = new char[len + 1];
        infile.read(buffer, len);
        buffer[len] = '\0';

        std::string result(buffer);
        delete[] buffer;

        return result;
    }

    std::vector<std::string> get_strings(int from, int to) {
        std::vector<std::string> results;

        // Open the file for reading
        std::ifstream infile(this->filename, std::ios::binary);
        
        // Skip metadata
        infile.seekg(METADATA_SIZE);

        // Calculate the position in the index table for the 'from' index
        int from_index_position = METADATA_SIZE + from * sizeof(int);

        // Seek to that position in the index table
        infile.seekg(from_index_position);

        // Read the position of the 'from' string
        int from_string_position;
        infile.read(reinterpret_cast<char*>(&from_string_position), sizeof(from_string_position));

        // Calculate the position in the index table for the 'to' index
        int to_index_position = METADATA_SIZE + to * sizeof(int);

        // Seek to that position in the index table
        infile.seekg(to_index_position);

        // Read the position of the 'to' string
        int to_string_position;
        infile.read(reinterpret_cast<char*>(&to_string_position), sizeof(to_string_position));

        // Now, calculate the total length of the chunk of strings we need to read
        int chunk_length = to_string_position - from_string_position;

        // Read the entire chunk of strings
        char* buffer = new char[chunk_length];
        infile.seekg(from_string_position);
        infile.read(buffer, chunk_length);

        // Split the chunk into individual strings based on their lengths and store them in the results vector
        char* ptr = buffer;
        char* end = buffer + chunk_length;
        while (ptr < end) {
            int len;
            memcpy(&len, ptr, sizeof(len));
            ptr += sizeof(len);

            std::string str(ptr, ptr + len);
            results.push_back(str);
            ptr += len;
        }

        delete[] buffer;
        return results;
    }

    std::string get_substring(int index, uint64_t pos, uint64_t size) {
        // Open the file for reading
        std::ifstream infile(this->filename, std::ios::binary);

        // Skip metadata
        infile.seekg(METADATA_SIZE);

        // Calculate the position in the index table for the given index
        int index_position = METADATA_SIZE + index * sizeof(int);

        // Seek to that position in the index table
        infile.seekg(index_position);

        // Read the position of the string
        int string_position;
        infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

        // Seek to the position of the string's length
        infile.seekg(string_position);

        // Read the length of the string
        int len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        
        // Check if the requested substring is valid
        if (pos < 0 || pos + size > len) {
            throw std::out_of_range("Requested substring is out of range");
        }

        // Calculate the position of the substring within the string and seek to it
        int substring_position = string_position + sizeof(int) + pos; // +sizeof(int) to account for the string's length
        infile.seekg(substring_position);

        // Read the substring
        char* buffer = new char[size + 1];
        infile.read(buffer, size);
        buffer[size] = '\0';

        std::string result(buffer);
        delete[] buffer;

        return result;
    }

};

class MemorySequenceAccess: public AbstractSequenceAccess{
    std::vector<std::string>& strings;
public:
    MemorySequenceAccess(std::vector<std::string>& strings): strings(strings){}

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    std::string get_string(int index) {
        return this->strings[index];
    }

    std::string get_substring(int index, uint64_t i, uint64_t size) {
        return this->strings[index].substr(i, size);
    }

};

#endif