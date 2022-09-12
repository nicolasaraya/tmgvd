#ifndef PCSA_H
#define PCSA_H
#include <iostream> 
#include <fstream>
#include <functional>
#include <math.h>
#include <vector> 
#include <bitset>
#include <unistd.h>
#include <omp.h>
#include <mutex> 
using namespace std; 


class PCSA{
    public:
        PCSA(const string, const unsigned char);
        ~PCSA();
        uint64_t compute();
    private:
        uint64_t R(uint64_t);
        void update(string);
        uint64_t estimation();
        /*****************/
        const double phi = 0.77351; 
        const double error = 0.05;
        const double m =(0.78)/error ; 
        const unsigned int M = ceil(m); 
        mutex* mtx;
        const unsigned char desplazamiento = 64 - log2(M);
        uint64_t* sketch;
        hash<string> h1;
        string pathFile; 
        unsigned char k;
};

#endif