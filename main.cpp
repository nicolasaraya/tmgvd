//g++ main.cpp -g -O3 -lm -std=c++11
#include <iostream>
#include <functional>
#include <fstream>
#include <omp.h>
#include "hyperloglog.h"
#include "PCSA.h"
#include "metrictime2.hpp"

using namespace std; 

const unsigned char k = 31; 
const string pathFile = "./files/GCF_000001405.39_GRCh38.p13_genomic.fna";
//const string pathFile = "./files/GCF_000308155.1_EptFus1.0_genomic.fna";


int main(int argc, char const *argv[]){

    TIMERSTART(PCSA);
    pcsa(pathFile, k);
    TIMERSTOP(PCSA);

    //cout << "si: " <<  __builtin_clz(128) << endl; 
    return 0;
}

