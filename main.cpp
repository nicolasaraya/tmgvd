//g++ main.cpp -g -O3 -lm -std=c++11
#include <iostream>
#include <functional>
#include <fstream>
#include <omp.h>
#include "hyperloglog.h"
#include "PCSA.h"

using namespace std; 

const unsigned char k = 31; 
const string pathFile = "./files/GCF_000001405.39_GRCh38.p13_genomic.fna";
//const string pathFile = "./files/GCF_000308155.1_EptFus1.0_genomic.fna";



int main(){
    ifstream file;
    file.open(pathFile);

    pcsa(&file, k);


        /*
    cout << "si: " <<  __builtin_clz(128) << endl; 
  


*/
}