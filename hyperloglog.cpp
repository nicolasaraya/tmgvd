#include "hyperloglog.hpp"

using namespace std; 

HLL::HLL(const string pathFile){
    cout << "pathFile = " << pathFile << endl;
    this->pathFile = pathFile;
    update_alpha();
    cout << "HLL params [ m: " << m << ", b: " << b << ", error: " << error << " ] " << endl; 
    
    freq.assign(m,0);
    N=0;
    for(int i = 0; i < m; i++) M_int_vector[i] = 0;
    
            
           
}

HLL::~HLL(){
    if(M != NULL) delete(M);
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
    return 64-pos-b;
    //return __builtin_ctzll(x);
}

void HLL::add(string kmer){
    N++;
    uint64_t x = h(kmer);
    //int  j = (~((1UL<< (64-b))-1) & x) >> (64-b); 
    uint32_t j = x >> (64-b);
    //unsigned long long int w = ((1UL<< (64-b))-1) & x;
    uint64_t w = (x << b) >> b ;
    mtx[j].lock();
    M_int_vector[j] = max(int(M_int_vector[j]),p(w));
    freq[j]++;
    mtx[j].unlock();
}

double HLL::estimate(){
    //cout << "computando " << endl;
    sdsl::util::bit_compress(M_int_vector);

    cout << "int_vector: " << size_in_mega_bytes(M_int_vector) << endl;
    
    sdsl::enc_vector<> M_enc(M_int_vector);
    cout << "enc_vector: " << size_in_mega_bytes(M_enc) << endl;
    
    sdsl::wt_int<sdsl::rrr_vector<b>> M_wt_int;
    sdsl::construct_im(M_wt_int, M_int_vector);
    cout << "wt_int: " << size_in_mega_bytes(M_wt_int) << endl;
    
    sdsl::wt_huff<sdsl::sd_vector<>> M_wt_huff;
    sdsl::construct_im(M_wt_huff, M_int_vector);
    cout << "wt_huff: " << size_in_mega_bytes(M_wt_huff) << endl;
    
    sdsl::wm_int<> M_wm_int;
    sdsl::construct_im(M_wt_int, M_int_vector);    
    cout << "wm_int: " << size_in_mega_bytes(M_wm_int) << endl;


    double sum = 0;
    int V = 0;
    for(int j=0 ; j < m ; j++){
        //cout << "sumando " <<M_int_vector[j] << endl;
        sum += pow(2,- int(M_int_vector[j]));
        if(M_int_vector[j] == 0) V++;
    }
    cout << "alpha: " << alpha << " ; m_^2: " << pow(m,2) << " ; suma: " << pow(sum,-1) << endl;
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

double HLL::entropy(){
    double entr = 0;
    //N = estimate();
    cout << N << endl;
    for(int i = 0; i < m; i++){
        double h = double(freq[i]) / double(N);
        entr -=  (h)*log2(h);
    }
    return entr;
}

double HLL::compute(){
    ifstream file;
    file.open(pathFile);
    if(!file.is_open()){
        cout << "Archivo no existe, " << pathFile << endl;
        exit(0); 
    }
    string line; 
    uint64_t countKmers = 0;
    //omp_set_num_threads(NUM_THREADS);
    //int count = 0;
    while(file >> line){
        //cout << "i: " << count++ <<", line " << line << endl; 
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
            if(line.size() <= k) {
                add(line);
                continue;
            }
            //cout << "line: " << line << endl;            countKmers += line.size() - k;
            countKmers += line.size() - k;
            //#pragma omp parallel for
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
                    else{
                        break;
                    }
                }
                add(kmer);
                //cout << "k-mer " << i << ": " << kmer << endl;
                kmer.clear();
            }
        }
    }
    file.close();
    //return 0;
    return estimate();
    
}
