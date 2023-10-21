#ifndef UTIL_STRING_ACCESS_HPP_
#define UTIL_STRING_ACCESS_HPP_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../prokrustean.hpp"

using namespace std;

/* for both disk-based, memory-based */
class AbstractSequenceAccess{
public:
	virtual void read_seq(SeqId index, std::string &string)=0;
    virtual void read_seq_substr(SeqId index, Pos i, Pos size, std::string &string)=0;
};

std::string PROKRUSTEAN_EVIDENCE = std::string("I am Prokrustes. Do not confuse me with the greek creature - Procrustes. I hate him.", 256);

struct ProkrusteanSequenceMetadata {
    bool loaded;
    std::string evidence;
    std::string prokrustean_file_name;
    int lmin;
    uint64_t sequence_count;
    uint64_t strata_count;
    // sequence_indics[i] -> starting position of sequence[i]
    std::streampos streampos_sequence_indices; 
    std::streampos streampos_sequence;
    std::streampos streampos_strata;
};

class DiskSequenceAccess: public AbstractSequenceAccess{
public:
    std::string filename;
    std::ofstream outfile;
    std::ifstream infile;
    bool loaded=false;
    vector<string> strings;
    ProkrusteanSequenceMetadata metadata;
    vector<streampos> sequence_start_positions;
    DiskSequenceAccess(std::string filename): filename(filename){
    }

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    // std::string get_metadata() {
    //     std::ifstream tempfile(filename, std::ios::binary);
    //     char metadata[METADATA_SIZE + 1];
    //     tempfile.read(metadata, METADATA_SIZE);
    //     metadata[METADATA_SIZE] = '\0';
    //     return std::string(metadata);
    // }

    void verify(){
        if(std::string(PROKRUSTEAN_EVIDENCE.c_str(), 256)!=this->metadata.evidence){ 
            assert(false);
        }
        // if (PROKRUSTEAN_EVIDENCE != get_metadata()){
        //     throw std::invalid_argument("the file does not exist or seems like not the sequence index generated with the prokrustean. "+filename);
        // }
    }

    void write_metadata(Prokrustean prokrustean) {
        assert(this->outfile.is_open());
        ////////////////////////////////////////////////
        // std::string evidence;
        // std::string prokrustean_file_name;
        // int lmin;
        // uint64_t sequence_count;
        // uint64_t strata_count;
        ////////////////////////////////////////////////
        this->metadata.loaded=true;
        this->metadata.evidence=PROKRUSTEAN_EVIDENCE;
        this->metadata.prokrustean_file_name=prokrustean.file_name;
        this->metadata.sequence_count=prokrustean.sequence_count;
        this->metadata.lmin=prokrustean.lmin;
        this->metadata.strata_count=prokrustean.stratum_count;

        int metadata_fixed_size = 256 + 256 + sizeof(metadata.sequence_count) + sizeof(metadata.lmin) + sizeof(metadata.strata_count) + sizeof(std::streampos) + sizeof(std::streampos) + sizeof(std::streampos);
        std::streampos start_position_of_sequence_positions = metadata_fixed_size; 
        this->metadata.streampos_sequence_indices=start_position_of_sequence_positions;

        // sequence content
        int metadata_and_sequence_meta = metadata_fixed_size + /* seq_pos*/ + prokrustean.sequences__size.size()*sizeof(streampos);
        std::streampos start_position_sequence=metadata_and_sequence_meta;
        this->metadata.streampos_sequence=start_position_sequence;

        // strata content
        uint64_t metadata_and_sequence_meta_and_content = metadata_and_sequence_meta + std::accumulate(prokrustean.sequences__size.begin(), prokrustean.sequences__size.end(), 0);
        std::streampos start_position_strata = metadata_and_sequence_meta_and_content;
        this->metadata.streampos_strata=start_position_strata;

        // sequence positions
        std::streampos seq_position = start_position_sequence;
        for (auto size : prokrustean.sequences__size) {
            this->sequence_start_positions.push_back(seq_position);
            seq_position += size + sizeof(SequenceSize);
        }

        outfile.write(PROKRUSTEAN_EVIDENCE.c_str(), 256);
        prokrustean.file_name.resize(256);
        outfile.write(prokrustean.file_name.c_str(), 256);
        outfile.write(reinterpret_cast<const char*>(&prokrustean.sequence_count), sizeof(prokrustean.sequence_count));
        outfile.write(reinterpret_cast<const char*>(&prokrustean.lmin), sizeof(prokrustean.lmin));
        outfile.write(reinterpret_cast<const char*>(&prokrustean.stratum_count), sizeof(prokrustean.stratum_count));
        outfile.write(reinterpret_cast<const char*>(&start_position_of_sequence_positions), sizeof(start_position_of_sequence_positions));
        outfile.write(reinterpret_cast<const char*>(&start_position_sequence), sizeof(start_position_sequence));
        outfile.write(reinterpret_cast<const char*>(&start_position_strata), sizeof(start_position_strata));
        for(auto pos: this->sequence_start_positions){
            outfile.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
        }
    }

    void load_metadata() {
        auto &metadata=this->metadata;
        std::ifstream _infile(this->filename, std::ios::binary);
        
        char evidence_buffer[256];
        _infile.read(evidence_buffer, 256);
        metadata.evidence = std::string(evidence_buffer, 256);
        
        char prokrustean_file_name_buffer[256];
        _infile.read(prokrustean_file_name_buffer, 256);
        metadata.prokrustean_file_name = std::string(prokrustean_file_name_buffer, 256);

        _infile.read(reinterpret_cast<char*>(&metadata.sequence_count), sizeof(metadata.sequence_count));
        _infile.read(reinterpret_cast<char*>(&metadata.lmin), sizeof(metadata.lmin));
        _infile.read(reinterpret_cast<char*>(&metadata.strata_count), sizeof(metadata.strata_count));
        _infile.read(reinterpret_cast<char*>(&metadata.streampos_sequence_indices), sizeof(metadata.streampos_sequence_indices));
        _infile.read(reinterpret_cast<char*>(&metadata.streampos_sequence), sizeof(metadata.streampos_sequence));
        _infile.read(reinterpret_cast<char*>(&metadata.streampos_strata), sizeof(metadata.streampos_strata));
        // cout << " metadata.sequence_count " << metadata.sequence_count << endl;
        // cout << " metadata.lmin " << metadata.lmin << endl;
        // cout << " metadata.strata_count " << metadata.strata_count << endl;
        // cout << " metadata.streampos_sequence_indices " << metadata.streampos_sequence_indices << endl;
        // cout << " metadata.streampos_sequence " << metadata.streampos_sequence << endl;
        // cout << " metadata.streampos_strata " << metadata.streampos_strata << endl;
        // if later becomes a burden, do not load but fetch ondemand.
        this->sequence_start_positions.resize(metadata.sequence_count);
        for (int i = 0; i < metadata.sequence_count; i++) {
            _infile.read(reinterpret_cast<char*>(&this->sequence_start_positions[i]), sizeof(streampos));
        }
        _infile.close();
        metadata.loaded=true;
    }


    // void write_meta(const std::vector<int>& string_sizes) {
    //     // Ensure the metadata string is the expected length
    //     this->outfile.write(PROKRUSTEAN_EVIDENCE.c_str(), METADATA_SIZE);

    //     // Compute and write the index table
    //     int position = METADATA_SIZE + sizeof(int) * string_sizes.size(); // meta + index table size + all strings size
    //     for (int size : string_sizes) {
    //         this->outfile.write(reinterpret_cast<char*>(&position), sizeof(position));
    //         position += size + sizeof(int); // move by string size + size of the integer that stores the length
    //     }
        
    //     // Pre-allocate space for all the strings
    //     outfile.seekp(position - 1);  // Move to the end of the pre-allocated space
    //     char dummy_byte = 0;
    //     outfile.write(&dummy_byte, 1);  // Write a single byte to make the file system allocate space up to this position
    // }

    void write_strings(const std::vector<std::string>& strings) {
        assert(metadata.streampos_sequence == this->sequence_start_positions[0]);
        // this->outfile.seekp(metadata.streampos_sequence);
        // Write the strings
        for (const std::string& str : strings) {
            SequenceSize len = str.size();
            this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));  // this could be optional now, as we have sizes saved
            this->outfile.write(str.c_str(), len);
        }
    }

    void write_open() {
        this->outfile=std::ofstream(filename, std::ios::binary);
    }

    // void write_single_sequence(SeqId id, const std::string& str){
    //     assert(metadata.loaded);
    //     assert(id<metadata.sequence_count);
    //     assert(metadata.streampos_sequence<=this->sequence_start_positions[id]);
    //     //do I need to recheck with the metadata?
    //     SequenceSize len = str.size();
    //     // this->outfile.seekp(this->sequence_start_positions[id]);
    //     this->outfile.write(reinterpret_cast<char*>(&len), sizeof(len));
    //     this->outfile.write(str.c_str(), len);
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

    void read_seq(SeqId id, std::string &sequence){
        assert(metadata.loaded);
        assert(id<metadata.sequence_count);
        // Seek to the position of the string
        infile.seekg(this->sequence_start_positions[id]);
        // Read the length of the string (assuming we stored lengths before the strings, as in the earlier code)
        SequenceSize len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        // Read the string itself
        char* buffer = new char[len + 1];
        infile.read(buffer, len);
        buffer[len] = '\0';

        sequence=std::string(buffer);
        delete[] buffer;
    }

    void read_seq_substr(SeqId id, Pos pos, Pos size, std::string &string){
        assert(metadata.loaded);
        assert(id<metadata.sequence_count);
        // Seek to that position in the index table
        infile.seekg(this->sequence_start_positions[id]);

        // Read the length of the string
        SequenceSize len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        
        // Check if the requested substring is valid
        if (pos < 0 || pos + size > len) {
            throw std::out_of_range("Requested substring is out of range");
        }

        // Calculate the position of the substring within the string and seek to it
        streampos substring_position =  sizeof(len) + pos; //
        infile.seekg(this->sequence_start_positions[id]+substring_position);

        // Read the substring
        char* buffer = new char[size + 1];
        infile.read(buffer, size);
        buffer[size] = '\0';

        string=std::string(buffer);
        delete[] buffer;
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

    // std::string get_string(int index) {
    //     if(this->loaded){
    //         return this->strings[index];
    //     }
    //     // Open the file for reading
    //     std::ifstream infile(this->filename, std::ios::binary);

    //     // Skip metadata
    //     infile.seekg(METADATA_SIZE);

    //     // Calculate the position in the index table for the given index
    //     int index_position = METADATA_SIZE + index * sizeof(int);

    //     // Seek to that position in the index table
    //     infile.seekg(index_position);

    //     // Read the position of the string
    //     int string_position;
    //     infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

    //     // Seek to the position of the string
    //     infile.seekg(string_position);

    //     // Read the length of the string (assuming we stored lengths before the strings, as in the earlier code)
    //     int len;
    //     infile.read(reinterpret_cast<char*>(&len), sizeof(len));

    //     // Read the string itself
    //     char* buffer = new char[len + 1];
    //     infile.read(buffer, len);
    //     buffer[len] = '\0';

    //     std::string result(buffer);
    //     delete[] buffer;

    //     return result;
    // }

    // std::vector<std::string> get_strings(int from, int to) {

    //     std::vector<std::string> results;

    //     // Open the file for reading
    //     std::ifstream infile(this->filename, std::ios::binary);
        
    //     // Skip metadata
    //     infile.seekg(METADATA_SIZE);

    //     // Calculate the position in the index table for the 'from' index
    //     int from_index_position = METADATA_SIZE + from * sizeof(int);

    //     // Seek to that position in the index table
    //     infile.seekg(from_index_position);

    //     // Read the position of the 'from' string
    //     int from_string_position;
    //     infile.read(reinterpret_cast<char*>(&from_string_position), sizeof(from_string_position));

    //     // Calculate the position in the index table for the 'to' index
    //     int to_index_position = METADATA_SIZE + to * sizeof(int);

    //     // Seek to that position in the index table
    //     infile.seekg(to_index_position);

    //     // Read the position of the 'to' string
    //     int to_string_position;
    //     infile.read(reinterpret_cast<char*>(&to_string_position), sizeof(to_string_position));

    //     // Now, calculate the total length of the chunk of strings we need to read
    //     int chunk_length = to_string_position - from_string_position;

    //     // Read the entire chunk of strings
    //     char* buffer = new char[chunk_length];
    //     infile.seekg(from_string_position);
    //     infile.read(buffer, chunk_length);

    //     // Split the chunk into individual strings based on their lengths and store them in the results vector
    //     char* ptr = buffer;
    //     char* end = buffer + chunk_length;
    //     while (ptr < end) {
    //         int len;
    //         memcpy(&len, ptr, sizeof(len));
    //         ptr += sizeof(len);

    //         std::string str(ptr, ptr + len);
    //         results.push_back(str);
    //         ptr += len;
    //     }

    //     delete[] buffer;
    //     return results;
    // }

    // std::string get_substring(int index, uint64_t pos, uint64_t size) {
    //     if(this->loaded){
    //         return this->strings[index].substr(pos, size);
    //     }
    //     // Open the file for reading
    //     std::ifstream infile(this->filename, std::ios::binary);

    //     // Skip metadata
    //     infile.seekg(METADATA_SIZE);

    //     // Calculate the position in the index table for the given index
    //     int index_position = METADATA_SIZE + index * sizeof(int);

    //     // Seek to that position in the index table
    //     infile.seekg(index_position);

    //     // Read the position of the string
    //     int string_position;
    //     infile.read(reinterpret_cast<char*>(&string_position), sizeof(string_position));

    //     // Seek to the position of the string's length
    //     infile.seekg(string_position);

    //     // Read the length of the string
    //     int len;
    //     infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        
    //     // Check if the requested substring is valid
    //     if (pos < 0 || pos + size > len) {
    //         throw std::out_of_range("Requested substring is out of range");
    //     }

    //     // Calculate the position of the substring within the string and seek to it
    //     int substring_position = string_position + sizeof(int) + pos; // +sizeof(int) to account for the string's length
    //     infile.seekg(substring_position);

    //     // Read the substring
    //     char* buffer = new char[size + 1];
    //     infile.read(buffer, size);
    //     buffer[size] = '\0';

    //     std::string result(buffer);
    //     delete[] buffer;

    //     return result;
    // }

};

class MemorySequenceAccess: public AbstractSequenceAccess{
    std::vector<std::string>& strings;
public:
    MemorySequenceAccess(std::vector<std::string>& strings): strings(strings){}

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    void read_seq(SeqId index, std::string &string){
        string=this->strings[index];
    }
    void read_seq_substr(SeqId index, Pos i, Pos size, std::string &string){
        string=this->strings[index].substr(i, size);
    }
};

#endif