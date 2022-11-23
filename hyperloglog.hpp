#ifndef HyperLogLog_h
#define HyperLogLog_h
#include <iostream> 
#include <math.h>
#include <fstream>
#include <sdsl/vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <vector>
#include <omp.h>
#include <mutex>
#include "Utils.hpp"

using namespace std; 


class HLL{
    public:
        HLL(const string);
        ~HLL();
        double compute();
        double estimate();
        double entropy();
    private:
        
        double alpha;
        mutex mtx[m];

        int* M = NULL; 
        sdsl::int_vector<> M_int_vector; 
        vector<int> freq; 
        uint64_t N;
        //sdsl::enc_vector M_enc_vector;


        hash<string> h;
        string pathFile; 
        /*********/
        void update_alpha();
        int p(uint64_t);
        void add(string); 
};



#endif



