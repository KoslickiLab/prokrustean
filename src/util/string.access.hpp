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
        std::ofstream outfile(filename, std::ios::binary);

        // Write metadata
        outfile.write(EVIDENCE.c_str(), METADATA_SIZE);

        // Write the strings
        for (const std::string& str : strings) {
            int len = str.size();
            outfile.write(reinterpret_cast<char*>(&len), sizeof(len));
            outfile.write(str.c_str(), len);
        }
    }

    std::string get_string(int index) {
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

    // void save_strings(const std::vector<std::string>& strings) {
    //     std::ofstream file(this->filename, std::ios::binary);

    //     // Write metadata.
        
    //     EVIDENCE.resize(METADATA_SIZE, '\0'); // Resize to fixed size, padding with zeros.
    //     file.write(EVIDENCE.c_str(), METADATA_SIZE);

    //     // Write all strings consecutively.
    //     for (const auto& str : strings) {
    //         file.write(str.c_str(), str.size());
    //     }

    //     file << "\n";

    //     // Record the start position of the index.
    //     // Write the index.
    //     uint64_t start = 0;
    //     for (const auto& str : strings) {
    //         uint64_t end = start + str.size() - 1;
    //         file << start << " " << end << "\n";
    //         start = end + 1;
    //     }

    //     file.close();
    // }

    // std::string read_metadata() {
    //     std::ifstream file(this->filename, std::ios::binary);
    //     std::string metadata(METADATA_SIZE, '\0');
    //     file.read(&metadata[0], METADATA_SIZE);
    //     return metadata;
    // }

    // std::string get_string(int index) {
    //     std::ifstream file(filename);

    //     // Skip metadata.
    //     file.ignore(METADATA_SIZE + 1); // +1 for the newline character.

    //     // Move to the indexed position.
    //     for (int k = 0; k < index; k++) {
    //         file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    //     }

    //     uint64_t start, end;
    //     file >> start >> end;

    //     uint64_t string_length = end - start + 1;
    //     char buffer[string_length + 1];

    //     file.seekg(METADATA_SIZE + 1 + start);
    //     file.read(buffer, string_length);
    //     buffer[string_length] = '\0';

    //     return std::string(buffer);
    // }


    // std::string get_substring(int index, uint64_t i, uint64_t j) {
    //     std::ifstream file(filename);

    //     // Skip metadata.
    //     file.ignore(METADATA_SIZE + 1); // +1 for the newline character.

    //     // Jump to the correct index.
    //     for (int k = 0; k < index; k++) {
    //         file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    //     }

    //     uint64_t start, end;
    //     file >> start >> end;

    //     // Adjust the start and end positions based on the requested substring.
    //     if (i > j || (start + j) > end) {
    //         return ""; // Invalid range
    //     }

    //     uint64_t substring_start = METADATA_SIZE + 1 + start + i;
    //     uint64_t substring_length = j - i + 1;
    //     char buffer[substring_length + 1];

    //     file.seekg(substring_start);
    //     file.read(buffer, substring_length);
    //     buffer[substring_length] = '\0';

    //     return std::string(buffer);
    // }


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