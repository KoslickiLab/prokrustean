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
    bool loaded=false;
    vector<string> strings;
    DiskSequenceAccess(std::string filename): filename(filename){}

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    std::string get_metadata() {
        std::ifstream infile(filename, std::ios::binary);
        char metadata[METADATA_SIZE + 1];
        infile.read(metadata, METADATA_SIZE);
        metadata[METADATA_SIZE] = '\0';
        return std::string(metadata);
    }

    bool verify(){
        return EVIDENCE == get_metadata();
    }

    void save_strings(const std::vector<std::string>& strings) {
        this->outfile=std::ofstream(filename, std::ios::binary);

        // Write metadata
        this->outfile.write(EVIDENCE.c_str(), METADATA_SIZE);

        // Write the strings
        for (const std::string& str : strings) {
            int len = str.size();
            this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));
            this->outfile.write(str.c_str(), len);
        }
        this->outfile.close();
    }

    void save_activate() {
        this->outfile=std::ofstream(filename, std::ios::binary);

        // Write metadata
        this->outfile.write(EVIDENCE.c_str(), METADATA_SIZE);
    }

    void save_single(const std::string& str){
        int len = str.size();
        this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));
        this->outfile.write(str.c_str(), len);
    }

    void save_deactivate() {
        this->outfile.close();
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
        if(this->loaded){
            return this->strings[index];
        }
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(METADATA_SIZE); // Skip metadata

        int len;
        for (int i = 0; i <= index; ++i) {
            infile.read(reinterpret_cast<char*>(&len), sizeof(len));
            if (i == index) {
                char* buffer = new char[len + 1];
                infile.read(buffer, len);
                buffer[len] = '\0';
                std::string result(buffer);
                delete[] buffer;
                return result;
            } else {
                infile.seekg(len, std::ios_base::cur);
            }
        }

        return ""; // Invalid index
    }

    std::string get_substring(int index, uint64_t pos, uint64_t size) {
        if(this->loaded){
            return this->strings[index].substr(pos, size);
        }
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(METADATA_SIZE); // Skip metadata

        int len;
        for (int i = 0; i <= index; ++i) {
            infile.read(reinterpret_cast<char*>(&len), sizeof(len));
            if (i == index) {
                if (pos + size > len || pos < 0) return ""; // Invalid range

                infile.seekg(pos, std::ios_base::cur); // Move to the start position of substring
                char* buffer = new char[size + 1];
                infile.read(buffer, size);
                buffer[size] = '\0';
                std::string result(buffer);
                delete[] buffer;
                return result;
            } else {
                infile.seekg(len, std::ios_base::cur);
            }
        }
        return ""; // Invalid index
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