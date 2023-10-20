#ifndef UTIL_DATA_STORE_HPP_
#define UTIL_DATA_STORE_HPP_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class AbstractStringDataStore{
public:
    virtual void store(std::string string)=0;
};

class DiskStringDataStore: public  AbstractStringDataStore{
    std::ofstream outfile;
    std::vector<std::string> buffer;
    const uint64_t batchSize=10000;
    
    void writeBuffer() {
        for (const auto& str : buffer) {
            outfile << str << '\n';
        }
        buffer.clear();
    }
public:
    DiskStringDataStore(string filename){
        outfile.open(filename);
        if (!outfile.is_open()) {
            throw std::runtime_error("Failed to open file for writing.");
        }
        buffer.reserve(batchSize);
    }

     ~DiskStringDataStore() {
        this->close();
    }

    void close(){
        if (!buffer.empty()) {
            writeBuffer();
        }
        outfile.close();
    }
    void store(std::string string){
        buffer.push_back(string);
        if (buffer.size() >= batchSize) {
            writeBuffer();
        }
    }
};


class MemoryStringDataStore: public AbstractStringDataStore{
public:
    vector<string> strings;

    void store(std::string string){
        this->strings.push_back(string);
        if(this->strings.size()%1000==0){
            cout << "string " << this->strings.size() <<endl;
        }
    }
    void reset(){
        this->strings.clear();
    }
};

#endif