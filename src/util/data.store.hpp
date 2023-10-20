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
    string filename;
    std::ofstream* file;
public:
    DiskStringDataStore(string filename):filename(filename){}
    void activate(){
        this->file=new std::ofstream(this->filename);
    }
    void deactivate(){
        this->file->close();
        delete this->file;
    }
    void store(std::string string){
        assert(this->file!=nullptr);
        (*this->file) << string << endl;
    }
};

class MemoryStringDataStore: public AbstractStringDataStore{
public:
    vector<string> strings;

    void store(std::string string){
        this->strings.push_back(string);
    }
    void reset(){
        this->strings.clear();
    }
};

#endif