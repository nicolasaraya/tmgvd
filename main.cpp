//g++ main.cpp -g -O3 -lm -std=c++11
#include <iostream>
#include <functional>
#include <fstream>
#include <omp.h>
#include "PCSA.h"
#include "hyperloglog.h"
#include "metrictime2.hpp"

using namespace std; 

const unsigned char k = 31; 
//const string pathFile = "/home/fly/bigdata/files/GCF_000001405.39_GRCh38.p13_genomic.fna";
const string pathFile = "./files/GCF_000308155.1_EptFus1.0_genomic.fna";
//const string pathFile = "./files/test2.fna";


double Jaccard(HLL* hll1, HLL* hll2){
    double A = hll1->estimate(); 
    double B = hll2->estimate();
    hll1->unionHLL(hll2->getSketch()); 
    double AuB = hll1->estimate();
    double AnB = A + B - AuB; 
    double jac = AnB/AuB;  
    return jac; 
}


int main(int argc, char const *argv[]){
    uint64_t resPCSA; 
    double resHLL; 
    
    TIMERSTART(_PCSA);
    PCSA* p = new PCSA(pathFile, k);
    resPCSA = p->compute();
    TIMERSTOP(_PCSA);
    if(p !=NULL) delete(p);

    TIMERSTART(_HLL);
    HLL* h = new HLL(pathFile, k);
    resHLL = h->compute();
    TIMERSTOP(_HLL);
    if(h !=NULL) delete(h);

    cout << "PCSA estimation: " << resPCSA <<endl;
    cout << "HyperLogLog estimation: " << resHLL << endl;

    return 0;
}

