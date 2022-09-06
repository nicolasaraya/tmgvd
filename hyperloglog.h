#ifndef HyperLogLog_h
#define HyperLogLog_h
#include <iostream> 
#include <math.h>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std; 


class HLL{
    public:
        HLL(const string, const unsigned char);
        ~HLL();
        double compute();
        bool unionHLL(int*); 
        int* getSketch();
    private:
        const int b = 6;
        const int m = pow(2,b);
        double alpha;
        int* M;
        hash<string> h;
        string pathFile; 
        unsigned char k;
        /*********/
        void update_alpha();
        int p(uint64_t);
        void add(string); 
        double estimate();
};
#endif