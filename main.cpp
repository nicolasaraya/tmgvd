//g++ main.cpp -g -O3 -lm -std=c++11
#include <iostream>
#include <functional>
#include <fstream>
#include <omp.h>
#include "hyperloglog.h"
#include "PCSA.h"
#include "hyperloglog.h"
#include "metrictime2.hpp"

using namespace std; 

const unsigned char k = 31; 
//const string pathFile = "/home/fly/bigdata/files/GCF_000001405.39_GRCh38.p13_genomic.fna";
const string pathFile = "./files/GCF_000308155.1_EptFus1.0_genomic.fna";
//const string pathFile = "./files/test.fna";


int main(int argc, char const *argv[]){
    uint64_t resPCSA; 
    double resHLL; 
    TIMERSTART(_PCSA);
    PCSA* p = new PCSA(pathFile, k);
    resPCSA = p->compute();
    TIMERSTOP(_PCSA);
    delete(p);

    TIMERSTART(_HLL);
    HLL* h = new HLL(pathFile, k);
    resHLL = h->compute();
    TIMERSTOP(_HLL);
    delete(h);

    cout << "PCSA estimation: " << resPCSA <<endl;
    cout << "HyperLogLog estimation: " << resHLL << endl;

    return 0;
}

