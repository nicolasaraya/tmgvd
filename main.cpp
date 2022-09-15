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
const string pathFile2 = "./files/GCF_000001405.39_GRCh38.p13_genomic.fna";

//const string pathFile = "./files/test2.fna";


double JaccardHLL(HLL* hll1, HLL* hll2){
    double A = hll1->estimate(); 
    double B = hll2->estimate();
    hll1->unionHLL(hll2->getSketch()); 
    double AuB = hll1->estimate();
    double AnB = A + B - AuB; 
    double jac = AnB/AuB;  
    return jac; 
}

double JaccardPCSA(PCSA* p, PCSA* p2){
    double A = p->estimation(); 
    double B = p2->estimation();
    p->unionPCSA(p2->getSketch()); 
    double AuB = p->estimation();
    double AnB = A + B - AuB; 
    double jac = AnB/AuB;  
    return jac; 
}


int main(int argc, char const *argv[]){
    uint64_t resPCSA; 
    double resHLL, resHLL2; 
    
    TIMERSTART(_PCSA);
    PCSA* p = new PCSA(pathFile, k);
    resPCSA = p->compute();
    TIMERSTOP(_PCSA);
    TIMERSTART(_PCSA2);
    PCSA* p2 = new PCSA(pathFile, k);
    resPCSA = p2->compute();
    TIMERSTOP(_PCSA2);

    TIMERSTART(_HLL);
    HLL* h = new HLL(pathFile, k);
    resHLL = h->compute();
    TIMERSTOP(_HLL);
    TIMERSTART(_HLL2);
    HLL* h2 = new HLL(pathFile2, k);
    resHLL2 = h2->compute();
    TIMERSTOP(_HLL2);

    cout << "PCSA estimation: " << resPCSA <<endl;
    cout << "HyperLogLog estimation: " << resHLL << ", " <<  resHLL2 << endl;
    cout << "JaccardPCSA: " << JaccardPCSA(p,p2) << endl; 
    cout << "JaccardHLL: " << JaccardHLL(h,h2) << endl; 


    
    if(h != NULL) delete(h);
    if(h2 != NULL) delete(h2);
    if(p != NULL) delete(p);
    if(p2 != NULL) delete(p2);
    return 0;
}

