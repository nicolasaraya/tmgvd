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
const string pathFile = "/home/fly/bigdata/files/GCF_000001405.39_GRCh38.p13_genomic.fna";
//const string pathFile = "/home/fly/bigdata/files/GCF_000308155.1_EptFus1.0_genomic.fna";


int main(int argc, char const *argv[]){

    TIMERSTART(PCSA);
    //pcsa(pathFile, k);
    hyper_loglog(pathFile,k);
    //double asd = double(5)/double(2);
    //cout << asd << endl;
    TIMERSTOP(PCSA);

    //cout << "si: " <<  __builtin_clz(128) << endl; 
    return 0;
}

