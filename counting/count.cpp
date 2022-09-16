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
                
            }
            kmer.clear();
        }
    }
    cout << "Cardinalidad: " << countKmer << endl;
	return 0;
}