#ifndef UTILS_HPP
#define UTILS_HPP

const unsigned char k = 31; 
const double error = 0.05;
const int temp_m = int(  double( pow(1.04, 2) ) / double( pow(error, 2) ) ) ; 
//const int m = 512;
const int m = pow(2,ceil(log2(temp_m)));

constexpr int b = floor(log2(m));
const int NUM_THREADS = 1;


#endif 