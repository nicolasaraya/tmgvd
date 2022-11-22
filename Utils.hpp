#ifndef UTILS_HPP
#define UTILS_HPP

const unsigned char k = 31; 
const double error = 0.10;
const int m = int(  double( pow(1.04, 2) ) / double( pow(error, 2) ) ) ; 
const int b = ceil(log2(m));


#endif 