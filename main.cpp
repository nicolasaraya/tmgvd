#include <iostream>
#include "hyperloglog.hpp"
#include "metrictime2.hpp"

using namespace std; 


int main(int argc, char const *argv[]){
    
    const string pathFile = argv[1];
    HLL* h;
    TIMERSTART(_HLL);
    h = new HLL(pathFile);
    double a =  h->compute();
    double b = h->entropy();
    cout << "Estimation: " << a << endl;
    cout << "Entropia: " << b << endl;
    TIMERSTOP(_HLL);

    //if(h != NULL) delete(h);
    return 0;

    
    /*
    sdsl::csa_wt<> fm_index;
    sdsl::construct_im(fm_index, "mississippi!", 1);
    cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
    sdsl::store_to_file(fm_index,"fm_index-file.sdsl");
    ofstream out("fm_index-file.sdsl.html");
    sdsl::write_structure<sdsl::HTML_FORMAT>(fm_index,out);
    */
    
}

