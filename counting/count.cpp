#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const string pathFile = "../files/preProcessed.fna";
//const string pathFile = "../files/file1.fna";
const unsigned char k = 31;


int main(){
	ifstream file;
    file.open(pathFile);
    string line; 
    map<string, int> mp; 

    uint64_t count = 0;
    uint64_t countKmer = 0;
    while(file >> line){
        if(count++ % 100000==0) cout << count << endl;
        if(line[0] != 'A' && line[0] != 'a' && line[0] != 'T' && line[0] != 't' && line[0] != 'C' && line[0] != 'c' && line[0] != 'G' && line[0] != 'g') continue; //linea no valida
        else if(line[1] != 'A' && line[1] != 'a' && line[1] != 'T' && line[1] != 't' && line[1] != 'C' && line[1] != 'c' && line[1] != 'G' && line[1] != 'g') continue; //linea no valida

        else{ //linea valida
            if(line.size() <= k) {
                continue;
            }
            //cout << "line: " << line << endl;
            #pragma omp parallel for
            for(size_t i = 0; i < line.size() - k; i++){
                string kmer; 
                for(size_t j = 0; j  < k; j++){
                   kmer += line[i+j];
                }   

                if(mp.count(kmer) > 0){ //si existe
                    //cout << kmer << " repetido " << countKmer << endl;
                }
                else{
                    mp.insert(make_pair(kmer, 1));
                    countKmer++;
                    //cout << kmer << "nuevo" << countKmer << endl;
                    
                }
                //cout << "k-mer " << i << ": " << kmer << endl;
                kmer.clear();
            }
            //cout << "nueva linea" << endl;
       }
    }

    cout << "Cardinalidad: " << countKmer << endl;

	
	return 0;
}