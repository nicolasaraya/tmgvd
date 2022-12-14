#include "hyperloglog.h"

using namespace std; 

HLL::HLL(const string pathFile, const unsigned char k){
    this->pathFile = pathFile;
    this->k = k; 
    M = new int[m];
    mtx = new mutex[m];
    for(int i = 0; i < m; i++) M[i] = 0;

    update_alpha();

}

HLL::~HLL(){
    delete(M);
    delete(mtx);
}

void HLL::update_alpha(){
    switch(m){
        case 16: alpha = 0.673;
            break;
        case 32: alpha = 0.697;
            break;
        case 64: alpha = 0.709;
            break;
        default: alpha = 0.7213 / (1.0 + 1.079 / double(m));
    }   
}

int HLL::p(uint64_t x){
    int pos = 0;
    while (x > 0) {
        x = x >> 1;
        pos++;
    }
    return 64-pos+1-b;
    return __builtin_ctzll(x);
}

void HLL::add(string kmer){
    uint64_t x = h(kmer);
    //int  j = (~((1UL<< (64-b))-1) & x) >> (64-b); 
    uint32_t j = x >> (64-b);
    //unsigned long long int w = ((1UL<< (64-b))-1) & x;
    uint64_t w = (x << b) >> b ;
    mtx[j].lock();
    M[j] = max(M[j],p(w));
    mtx[j].unlock();
}


double HLL::estimate(){
    //cout << "computando " << endl;
    double sum = 0;
    int V = 0;
    for(int j=0 ; j < m ; j++){
        sum+= pow(2,-M[j]);
        //cout << "sum: " << sum << endl;
        if(M[j] == 0) V++;
    }
    //cout << "alpha: " << alpha << " ; m_^2: " << pow(m,2) << endl;
    double E = alpha*pow(m,2)*pow(sum,-1);
 
    //cout << "resultado: " << E << endl;
    if( (E <= 2.5*m) && V!= 0){
        //cout << "primera condicion" << endl;
        E = m * log2(m/V);
    }
    if(E > (1.0/30.0)*(pow(2,64))){
        //cout << "segunda condicion" << endl;
        E = -1 * pow(2,64)*log2( 1- ( E / pow(2,64)) );
    }
    return E;
}

double HLL::compute(){
    ifstream file;
    file.open(pathFile);
    string line; 
    uint64_t count = 0;
    uint64_t countKmers = 0;
    omp_set_num_threads(7);

    while(file >> line){
        //if(count++ % 100000 == 0) cout << count << endl;
        
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
            if(line.size() <= k) {
                add(line);
                continue;
            }
            //cout << "line: " << line << endl;            countKmers += line.size() - k;
            countKmers += line.size() - k;
            #pragma omp parallel for
            for(size_t i = 0; i <= line.size() - k; i++){
                string kmer; 
                int l = 0; 
                for(size_t j = 0; j <= k; j++){

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
    cout << "Complejidad HLL" << m * log2(log2(countKmers)) << endl;


    return estimate();

}


bool HLL::unionHLL(int* M2){
    for(int i = 0; i < m; i++){
        M[i] = max(M[i], M2[i]);
    }
    return true;
}

int* HLL::getSketch(){
    return M;
}