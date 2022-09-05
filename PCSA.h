#ifndef PCSA_H
#define PCSA_H
#include <iostream> 
#include <fstream>
#include <functional>
#include <math.h>
#include <vector> 
#include <bitset>
#include <unistd.h>
#include <omp.h>

using namespace std; 

const double phi = 0.77351; 
const double error = 0.05;
const double m = (0.78)/error; 
///
const int M = ceil(m); 
const int desplazamiento = 64 - log2(M);
long long int* sketch = new long long int[M];
hash<string> h1;

long R(long x){
    return (~x & (x+1)); 
}

void update(string kmer){
    unsigned long long int x = h1(kmer); 
    bitset<64> x_bit(x);

    int y = x >> desplazamiento;

    /*
    bitset<64> y_bit(y);
    cout << "x: " << x << " / kmer: " << kmer << endl;
    cout << x_bit << endl;
    cout << y_bit << endl;
    cout << y << endl;
    cout << "*******" << endl;
    */
    sketch[y] = sketch[y] | R(x);
}

long estimation(long long int* sketch){
    int sum = 0;
    for(int i = 0; i < M; i++) {
        if(R(sketch[i]) - 1 >= 0){
            sum += log2(R(sketch[i]));
        }
    }
    double media = 1.0*sum/M;
    return M * pow(2,media) / phi; 
}

void pcsa(const string pathFile, const unsigned char k){
    ifstream file;
    file.open(pathFile);
    string line; 

    for(int i = 0; i < M; i++) sketch[i] = 0;

    omp_set_num_threads(1);

    while(file >> line){
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
             
            if(line.size() <= k) {
                update(line);
                continue;
            }
            //cout << "line: " << line << endl;
            #pragma omp parallel for
            for(int i = 0; i <= line.size() - k; i++){
                string kmer; 
                int l = 0; 
                for(int j = 0; j <= k; j++){

                   if(i+j+l < line.size()){
                            char aux = line[i+j+l];
                            if(aux == 'C' || aux == 'A' || aux == 'G' || aux == 'T'){
                                kmer += aux ;
                            }
                            else if(aux - 32 == 'C' || aux - 32== 'A' || aux -32 == 'G' || aux -32 == 'T'){
                                kmer += aux - 32 ;
                            }
                            else{
                                j--;
                                l++;
                            }
                    }
                    else break;
                }
                update(kmer);
                //cout << "k-mer " << i << ": " << kmer << endl;
                kmer.clear();
            }
       }
    }
    cout << "Cardinalidad estimada: " << estimation(sketch) << endl; 
}

#endif