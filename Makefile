CC=g++
CFLAGS=-g -O0 -fopenmp -lm -Wall -std=c++11 -DNDEBUG -I sdsl_lib/include/  -L sdsl_lib/lib/ 
LSDSLFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 
EXECUTABLE=main.o
CONFIGURE=clear mk
OBJECTS = main.cpp hyperloglog.cpp hyperloglog.hpp metrictime2.hpp

all: $(EXECUTABLE)

configure: $(CONFIGURE)

mk:
	rm -rf data
	mkdir -p ./data
	wget  -O ./data/file.fna.gz http://www.inf.udec.cl/~chernand/datalabs/genomas/bacterias/GCF_001522075.1_ASM152207v1_genomic.fna.gz
	cd data && gzip -d file.fna.gz && cd ..

main.o: 
	$(CC) $(CFLAGS) $(OBJECTS) $(LSDSLFLAGS) -o hll.o

clear:
	rm -fr hll
	rm -rf data
	rm -rf sdsl_lib
	rm -rf sdsl-lite
	rm -r hll.o

install-sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	mkdir sdsl_lib
	./sdsl-lite/install.sh ./sdsl_lib/
	rm -rf sdsl-lite
