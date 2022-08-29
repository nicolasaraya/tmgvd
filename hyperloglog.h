#ifndef HyperLogLog_h
#define HyperLogLog_h
#include <iostream> 
#include <math.h>
using namespace std; 

const int b = 4;
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

int p(long long int x){
    int pos = __builtin_clz(x);
    return pos + 1;
}

double compute(int* M_){
    double sum = 0;
    for(int j=0 ; j < m_ ; j++){
        sum+= pow(2,-M_[j]);
    }
    //double E = alpha*pow(m_,2)*pow(sum,-1);
    double E = 2;
    return E;
}

void hyper_loglog(const string pathFile, const unsigned char k){

    update_alpha();
    ifstream file;
    file.open(pathFile);
    string line; 

    int i = 0;
    while(file >> line){
        if(i == 300){
            break;
        }
        cout << line << endl;
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
            
            //cout << "line: " << line << endl;
            for(int i = 0; i <= line.size() - k; i++){
                string kmer; 
                int l = 0; 
                for(int j = 0; j < k; j++){
                    kmer += line[i+j];
                }
                cout << "kmer: " << kmer << endl;

                long long int x = h(kmer);

                //j = b
                //w =  64-b
                //M_[j] = max(M_[j],p(w))

            }
        }
    }
    compute(M_);

}
#endif