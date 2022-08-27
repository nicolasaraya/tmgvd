#ifndef PCSA_H
#define PCSA_H
#include <iostream> 
#include <fstream>
#include <functional>
#include <math.h>
#include <vector> 
#include <bitset>
#include <unistd.h>



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
    long long int x = h1(kmer); 
    int y = x >> desplazamiento;
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



void pcsa(ifstream* file, const unsigned char k){
    string line; 
    for(int i = 0; i < M; i++) sketch[i] = 0;

    int countLines = 0; 
    while(std::getline(*file, line)){
        if(countLines >= 2883199){
            cout << "line: " << line << endl;
            cout << line.size() << endl;
            cout << countLines << endl;
        }
        if(countLines%100 == 0) cout << countLines << endl;
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        
        else{ //linea valida
            
            //cout << line << endl;
            if(line.size() - k <= 0 ){
                continue;
            }    
            //grande pacifico conchetumare txiuuu
            for(int i = 0; i < line.size() - k; i++){
                string kmer; 

                
                
                for(int j = 1; j < k; j++){
                    if(i+j < line.size()){
                        kmer += line[i+j]; 
                    }
                    else{
                        break;
                    }
                
                }
                
                //update(kmer);
                
            }
        }
        countLines++;
        line.clear();
    }
    cout << "ok" << endl;
    cout << "Cardinalidad estimada: " << estimation(sketch) << endl; 



}

#endif