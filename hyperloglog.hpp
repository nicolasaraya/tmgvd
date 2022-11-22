#ifndef HyperLogLog_h
#define HyperLogLog_h
#include <iostream> 
#include <math.h>
#include <fstream>
#include <omp.h>
#include <mutex>
#include <sdsl/vectors.hpp>
#include "Utils.hpp"

using namespace std; 


class HLL{
    public:
        HLL(const string, unsigned char);
        ~HLL();
        double compute();
        int* getSketch();
        double estimate();
        bool unionHLL(int*);
    private:
        
        double alpha;
        mutex* mtx;

        void* M_ptr;

        int* M;
        hash<string> h;
        string pathFile; 
        unsigned char option; 
        /*********/
        void update_alpha();
        int p(uint64_t);
        void add(string); 
};



#endif



