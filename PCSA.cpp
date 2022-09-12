#include "PCSA.h"

PCSA::PCSA(const string pathFile, const unsigned char k){
    this->pathFile = pathFile; 
    this->k = k; 
    sketch = new uint64_t[M];
    for(size_t i = 0; i < M; i++) sketch[i] = 0;

}

PCSA::~PCSA(){
    delete(sketch);
}

uint64_t PCSA::R(uint64_t x){
     return (~x & (x+1)); 
}

void PCSA::update(string kmer){
    uint64_t x = h1(kmer); 
    unsigned char y = x >> desplazamiento;
    sketch[y] = sketch[y] | R(x);
}

uint64_t PCSA::estimation(){
    uint32_t sum = 0;
    for(size_t i = 0; i < M; i++) {
        if(R(sketch[i]) - 1 >= 0){
            sum += log2(R(sketch[i]));
        }
    }
    double media = 1.0*sum/M;
    return M * pow(2,media) / phi; 
}

uint64_t PCSA::compute(){
    cout << m << ","  <<  M << endl;
    ifstream file;
    file.open(pathFile);
    string line; 

    omp_set_num_threads(7);
    while(file >> line){
        //if(count++ % 100000==0) cout << count << endl;
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
             
            if(line.size() <= k) {
                update(line);
                continue;
            }
            //cout << "line: " << line << endl;
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
                update(kmer);
                //cout << "k-mer " << i << ": " << kmer << endl;
                kmer.clear();
            }
       }
    }
    return estimation();
}