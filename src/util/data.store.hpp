#ifndef UTIL_DATA_STORE_HPP_
#define UTIL_DATA_STORE_HPP_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../data_types.hpp"

using namespace std;

class AbstractStringDataStore{
public:
    virtual void store(std::string string)=0;
    virtual void chop_and_store(std::string string, int k)=0;
};

class DiskStringDataStore: public  AbstractStringDataStore{
    std::ofstream outfile;
    std::vector<std::string> buffer;
    const uint64_t batchSize=10000;
    SpinLock lock;
    
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
        this->lock.lock();
            outfile << string << '\n';
            // buffer.push_back(string);
            // if (buffer.size() >= batchSize) {
            //     writeBuffer();
            // }
        this->lock.unlock();
    }
    
    void chop_and_store(std::string string, int k){
        if(string.size()<k){
            return;
        }
        this->lock.lock();
            for(int p=0; p<string.size()-(k-1); p++){
                outfile << string.substr(p, k) << '\n';
            }
        this->lock.unlock();
    }
};


class MemoryStringDataStore: public AbstractStringDataStore{
public:
    vector<string> strings;

    void store(std::string string){
        this->strings.push_back(string);
    }
    void chop_and_store(std::string string, int k){
        if(string.size()<k){
            return;
        } else {
            for(int p=0; p<string.size()-(k-1); p++){
                this->strings.push_back(string.substr(p, k));
            }
        }
    }
    void reset(){
        this->strings.clear();
    }
};

#endif