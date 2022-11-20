CC=g++
CFLAGS=-g -O0 -fopenmp -lm -Wall -std=c++11 -DNDEBUG -I ./sdsl_lib/include -L ./sdsl_lib/lib -lsdsl -ldivsufsort -ldivsufsort64 
EXECUTABLE=main.o
CONFIGURE=clear mk
OBJECTS = main.cpp hyperloglog.cpp hyperloglog.hpp metrictime2.hpp

all: $(EXECUTABLE)

configure: $(CONFIGURE)

mk:
	mkdir -p ./data
	wget  -O ./data/file.fna.gz http://www.inf.udec.cl/~chernand/datalabs/genomas/bacterias/GCF_000717965.1_ASM71796v1_genomic.fna.gz
	cd data && gzip -d file.fna.gz && cd ..

main.o: 
	$(CC) $(CFLAGS) -o hll $(OBJECTS)

clear:
	rm -fr hll
	rm -rf data

install-sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	mkdir sdsl_lib
	./sdsl-lite/install.sh ./sdsl_lib/
	rm -rf sdsl-lite
