#include <iostream>
#include <functional>
#include <fstream>
#include <omp.h>
#include "hyperloglog.hpp"
#include "metrictime2.hpp"
#include <sdsl/suffix_arrays.hpp>

using namespace std; 

const unsigned char k = 31; 
const string pathFile = "./data/file.fna";


double JaccardHLL(HLL* hll1, HLL* hll2){
    double A = hll1->estimate(); 
    double B = hll2->estimate();
    hll1->unionHLL(hll2->getSketch()); 
    double AuB = hll1->estimate();
    double AnB = A + B - AuB; 
    double jac = AnB/AuB;  
    return jac; 
}


int main(int argc, char const *argv[]){
    /*
    HLL* h;


    TIMERSTART(_HLL);
    h = new HLL(pathFile, k);
    TIMERSTOP(_HLL);

    if(h != NULL) delete(h);
    return 0;

    */

    sdsl::csa_wt<> fm_index;
    sdsl::construct_im(fm_index, "mississippi!", 1);
    cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
    sdsl::store_to_file(fm_index,"fm_index-file.sdsl");
    ofstream out("fm_index-file.sdsl.html");
    sdsl::write_structure<sdsl::HTML_FORMAT>(fm_index,out);
    
}

