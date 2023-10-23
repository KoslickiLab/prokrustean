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

    void print(){
        cout << " meta: " << endl;
        cout << " evidence " << evidence << endl;
        cout << " prokrustean_file_name " << prokrustean_file_name << endl;
        cout << " lmin " << lmin << endl;
        cout << " sequence_count " << sequence_count << endl;
        cout << " strata_count " << strata_count << endl;
        cout << " streampos_sequence_indices " << streampos_sequence_indices << endl;
    }
};

class DiskSequenceAccess: public AbstractSequenceAccess{
public:
    std::string filename;
    std::ofstream writefile;
    std::ofstream updatefile;
    std::ifstream infile;
    bool loaded=false;
    vector<string> strings;
    ProkrusteanSequenceMetadata metadata;
    vector<streampos> sequence_start_positions;
    DiskSequenceAccess(std::string filename): filename(filename){}

    const uint64_t METADATA_SIZE = 256; // 256 bytes for metadata

    void verify(){
        if(std::string(PROKRUSTEAN_EVIDENCE.c_str(), 256)!=this->metadata.evidence){ 
            assert(false);
        }
        // if (PROKRUSTEAN_EVIDENCE != get_metadata()){
        //     throw std::invalid_argument("the file does not exist or seems like not the sequence index generated with the prokrustean. "+filename);
        // }
    }

    void write_metadata(vector<SequenceSize> &seq_sizes, optional<Prokrustean> prokrustean=nullopt) {
        assert(this->writefile.is_open());
        ////////////////////////////////////////////////
        // std::string evidence;
        // std::string prokrustean_file_name;
        // int lmin;
        // uint64_t sequence_count;
        // uint64_t strata_count;
        ////////////////////////////////////////////////
        this->metadata.loaded=true;
        this->metadata.evidence=PROKRUSTEAN_EVIDENCE;
        this->metadata.prokrustean_file_name=prokrustean.has_value()? prokrustean.value().file_name: "not set";
        this->metadata.sequence_count=seq_sizes.size();
        this->metadata.lmin=prokrustean.has_value()? prokrustean.value().lmin: 0;
        this->metadata.strata_count=prokrustean.has_value()? prokrustean.value().stratum_count: 0;

        // int metadata_fixed_size = 256 + 256 + sizeof(metadata.sequence_count) + sizeof(metadata.lmin) + sizeof(metadata.strata_count) + sizeof(std::streampos) + sizeof(std::streampos) + sizeof(std::streampos);
        // std::streampos start_position_of_sequence_positions = metadata_fixed_size; 
        // this->metadata.streampos_sequence_indices=start_position_of_sequence_positions;

        // // sequence content
        // int metadata_and_sequence_meta = metadata_fixed_size + /* seq_pos*/ + seq_sizes.size()*sizeof(streampos);
        // std::streampos start_position_sequence=metadata_and_sequence_meta;
        // this->metadata.streampos_sequence=start_position_sequence;

        // // strata content
        // uint64_t metadata_and_sequence_meta_and_content = metadata_and_sequence_meta + std::accumulate(seq_sizes.begin(), seq_sizes.end(), 0);
        // std::streampos start_position_strata = metadata_and_sequence_meta_and_content;
        // this->metadata.streampos_strata=start_position_strata;

        // sequence positions
        // std::streampos seq_position = start_position_sequence;
        // for (auto size : seq_sizes) {
        //     this->sequence_start_positions.push_back(seq_position);
        //     seq_position += size + sizeof(SequenceSize);
        // }

        writefile.write(PROKRUSTEAN_EVIDENCE.c_str(), 256);
        this->metadata.prokrustean_file_name.resize(256);
        writefile.write(this->metadata.prokrustean_file_name.c_str(), 256);
        writefile.write(reinterpret_cast<const char*>(&this->metadata.sequence_count), sizeof(this->metadata.sequence_count));
        writefile.write(reinterpret_cast<const char*>(&this->metadata.lmin), sizeof(this->metadata.lmin));
        writefile.write(reinterpret_cast<const char*>(&this->metadata.strata_count), sizeof(this->metadata.strata_count));
        std::streampos startPos = writefile.tellp();
        startPos+=sizeof(std::streampos); // point to right after the variable
        writefile.write(reinterpret_cast<const char*>(&startPos), sizeof(startPos));
        cout << " startPos " << startPos << endl;
        // writefile.write(reinterpret_cast<const char*>(&start_position_of_sequence_positions), sizeof(start_position_of_sequence_positions));
        // writefile.write(reinterpret_cast<const char*>(&start_position_sequence), sizeof(start_position_sequence));
        // writefile.write(reinterpret_cast<const char*>(&start_position_strata), sizeof(start_position_strata));
        // for(auto pos: this->sequence_start_positions){
        //     writefile.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
        // }
    }



    void write_strings(const std::vector<std::string>& strings) {
        // this->outfile.seekp(metadata.streampos_sequence);
        // Write the strings
        // Placeholder for string positions
        std::streampos startPos = writefile.tellp();
        for (uint64_t i = 0; i < strings.size(); ++i) {
            std::streampos placeholder = 0;
            writefile.write(reinterpret_cast<const char*>(&placeholder), sizeof(placeholder));
        }

        std::vector<std::streampos> positions;
        for (const auto& str : strings) {
            positions.push_back(writefile.tellp());

            SequenceSize len = str.size();
            writefile.write(reinterpret_cast<const char*>(&len), sizeof(len));
            writefile.write(str.c_str(), len);
        }

        // Go back and write the actual positions
        writefile.seekp(startPos);
        for (const auto& pos : positions) {
            writefile.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
        }

        // streampos seq_pos=this->writefile.tellp();
        // for (const std::string& str : strings) {
        //     this->writefile.write(reinterpret_cast<const char*>(&seq_pos), sizeof(seq_pos));
        //     streampos added=sizeof(SequenceSize)+str.size()*sizeof();
        //     seq_pos+=added;
        // }

        // for (const std::string& str : strings) {
        //     SequenceSize len = str.size();
        //     this->writefile.write(reinterpret_cast<const char*>(&len), sizeof(len));  // this could be optional now, as we have sizes saved
        //     this->writefile.write(str.c_str(), len);
        //     // cout << " save len " << (int)len << endl;
        // }
        
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
        // _infile.read(reinterpret_cast<char*>(&metadata.streampos_sequence), sizeof(metadata.streampos_sequence));
        // _infile.read(reinterpret_cast<char*>(&metadata.streampos_strata), sizeof(metadata.streampos_strata));
        
        // if later becomes a burden, do not load but fetch ondemand.
        // this->sequence_start_positions.resize(metadata.sequence_count);
        // for (int i = 0; i < metadata.sequence_count; i++) {
        //     _infile.read(reinterpret_cast<char*>(&this->sequence_start_positions[i]), sizeof(streampos));
        // }
        // this->sequence_start_positions.resize(metadata.sequence_count);
        // for (int i = 0; i < metadata.sequence_count; i++) {
        //     this->sequence_start_positions[i]=_infile.tellg();
        //     SequenceSize len;
        //     _infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        //     _infile.seekg(len, std::ios::cur);
        // }
        _infile.close();
        metadata.loaded=true;
    }
    void write_open() {
        this->writefile=std::ofstream(filename, std::ios::binary);
    }

    void write_close() {
        this->writefile.close();
    }

    void update_open() {
        this->updatefile=std::ofstream(filename, std::ios::binary | std::ios::in | std::ios::out);
    }

    void update_single_sequence(SeqId id, const std::string& str){
        assert(metadata.loaded);
        assert(id<metadata.sequence_count);
        assert(metadata.streampos_sequence<=this->sequence_start_positions[id]);
        assert(updatefile.is_open());
        //do I need to recheck with the metadata?
        SequenceSize len = str.size();
        this->updatefile.seekp(this->sequence_start_positions[id]);
        this->updatefile.write(reinterpret_cast<char*>(&len), sizeof(len));
        this->updatefile.write(str.c_str(), len);
    }
    void update_close() {
        this->updatefile.close();
    }

    void read_open() {
        this->infile=std::ifstream(filename, std::ios::binary);
    }

    void read_close() {
        this->infile.close();
    }

    void read_seq(SeqId id, std::string &sequence){
        if(this->loaded){
            sequence=this->strings[id];
            return;
        }
        assert(metadata.loaded);
        assert(id<metadata.sequence_count);
        // Seek to the position of the string
        streampos relative_pos = id * sizeof(std::streampos);
        infile.seekg(this->metadata.streampos_sequence_indices+relative_pos);
        // now at the position of the position pointing to the string
        streampos sequnece_pos;
        infile.read(reinterpret_cast<char*>(&sequnece_pos), sizeof(sequnece_pos));
        infile.seekg(sequnece_pos);
        SequenceSize len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        sequence.resize(len);
        infile.read(&sequence[0], len);
    }

    void read_seq_substr(SeqId id, Pos pos, Pos size, std::string &string){
        if(this->loaded){
            string=this->strings[id].substr(pos,size);
            return;
        }
        assert(metadata.loaded);
        assert(id<metadata.sequence_count);
        // Seek to that position in the index table
        streampos relative_pos = id * sizeof(std::streampos);
        infile.seekg(this->metadata.streampos_sequence_indices+relative_pos);
        // now at the position of the position pointing to the string
        streampos sequnece_pos;
        infile.read(reinterpret_cast<char*>(&sequnece_pos), sizeof(sequnece_pos));
        infile.seekg(sequnece_pos);
        // infile.seekg(this->sequence_start_positions[id]);

        // Read the length of the string
        SequenceSize len;
        infile.read(reinterpret_cast<char*>(&len), sizeof(len));
        
        // Check if the requested substring is valid
        if (pos < 0 || pos + size > len) {
            cout << "pos: " << pos << " size " << size << " but len " << len << " id " << id << endl;
            throw std::out_of_range("Requested substring is out of range");
        }

        // Calculate the position of the substring within the string and seek to it
        // streampos substring_position =  sizeof(len) + pos; //
        infile.seekg(pos, std::ios::cur);
        // infile.seekg(this->sequence_start_positions[id]+substring_position, std::ios_base::cu);

        string.resize(size);
        infile.read(&string[0], size);
        // // Read the substring
        // char* buffer = new char[size + 1];
        // infile.read(buffer, size);
        // buffer[size] = '\0';

        // string=std::string(buffer);
        // delete[] buffer;
    }

    void load_all_strings() {
        std::ifstream _infile(filename, std::ios::binary);
        _infile.seekg(metadata.streampos_sequence_indices);  // Skip the metadata
        streampos first_sequnece_pos;
        _infile.read(reinterpret_cast<char*>(&first_sequnece_pos), sizeof(first_sequnece_pos));
        _infile.seekg(first_sequnece_pos);
        
        this->strings.resize(this->metadata.sequence_count);
        uint64_t i=0;
        while (!_infile.eof()&& i<this->metadata.sequence_count) {
            SequenceSize len;
            _infile.read(reinterpret_cast<char*>(&len), sizeof(len));
            if (_infile.eof()) break;

            this->strings[i].resize(len);
            _infile.read(&this->strings[i][0], len);
            i++;
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

    void read_seq(SeqId index, std::string &string){
        string=this->strings[index];
    }
    void read_seq_substr(SeqId index, Pos i, Pos size, std::string &string){
        string=this->strings[index].substr(i, size);
    }
};

#endif