#ifndef HyperLogLog_h
#define HyperLogLog_h
#include <iostream> 
#include <math.h>
using namespace std; 
#include <bits/stdc++.h>
const int b = 6;
const int m_ = pow(2,b);
double alpha;

int* M_ = new int[m_];
hash<string> h;

void update_alpha(){
    switch(m_){
        case 16: alpha = 0.673;
            break;
        case 32: alpha = 0.697;
            break;
        case 64: alpha = 0.709;
            break;
        default: alpha = 0.7213 / (1.0 + 1.079 / double(m_));
    }   
}

int p(unsigned long long int x){
    int pos = 0;
    while (x > 0) {
        x = x >> 1;
        pos++;
    }
    return 64-pos+1-b;
    
    //bitset<64> x_bit(x);
    //cout << x_bit << endl;
    //cout << __builtin_ctzll(x) << endl;
    //cout << "*  " << endl;
    return __builtin_ctzll(x);
}

double compute(int* M_){
    cout << "computando " << endl;
    double sum = 0;
    int V = 0;
    for(int j=0 ; j < m_ ; j++){
        sum+= pow(2,-M_[j]);
        cout << "sum: " << sum << endl;
        if(M_[j] == 0) V++;
    }
    cout << "alpha: " << alpha << " ; m_^2: " << pow(m_,2) << endl;
    double E = alpha*pow(m_,2)*pow(sum,-1);
    //double E = 2;
    cout << "resultado: " << E << endl;
    if( (E <= (5/2)*m_) && V!= 0){
        cout << "primera condicion" << endl;
        E = m_ * log2(m_/V);
    }
    if(E > ((double)1/(double)30)*(pow(2,64))){
        cout << "segunda condicion" << endl;
        E = -1 * pow(2,64)*log2( (1-E) / pow(2,64) );
    }
    return E;
}
void add(string kmer){
    unsigned long long int x = h(kmer);
    //int  j = (~((1UL<< (64-b))-1) & x) >> (64-b); 
    int j = x >> (64-b);
    //unsigned long long int w = ((1UL<< (64-b))-1) & x;
    unsigned long long int w = (x << b) >> b ;

    M_[j] = max(M_[j],p(w));
}

void hyper_loglog(const string pathFile, const unsigned char k){

    update_alpha();
    ifstream file;
    file.open(pathFile);
    string line; 

    for(int i = 0; i < m_; i++) M_[i] = 0;
    int i = 0;
    omp_set_num_threads(4);

    while(file >> line){
        i++;
        if(i % 100000==0){
            cout << i << endl;
        }
        //cout << line << endl;
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
            if(line.size() <= k) {
                add(line);
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
                add(kmer);
                //cout << "k-mer " << i << ": " << kmer << endl;
                kmer.clear();
            }
        }
    }
    cout << "Cardinalidad estimada: " << compute(M_) << endl; 

}
#endif