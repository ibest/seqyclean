CXX = g++
CFLAGS = -Wall -g -L.
AR       = ar cr
SRC = src/
BIN = bin/
OBJ = obj/
LIBRARY := ${OBJ}lgzstream.a

all:  mkobj mkbin gzstream.o libgzstream.a abi.o poly.o sff.o sffreader.o ascii.o util.o Read.o QualTrim.o Report.o iz_SSAHA.o pairwise.o Dictionary.o KMerRoutine.o MainPipeLine.o Illumina.o Roche.o dup.o flash.o main.o seqyclean 
		
#					
seqyclean :   $(OBJ)main.o $(OBJ)flash.o $(OBJ)dup.o $(OBJ)Roche.o $(OBJ)Illumina.o $(OBJ)MainPipeLine.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)pairwise.o $(OBJ)iz_SSAHA.o $(OBJ)Report.o $(OBJ)QualTrim.o $(OBJ)sffreader.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sffreader.o $(OBJ)sff.o $(OBJ)abi.o $(OBJ)gzstream.o
	$(CXX) $(CFLAGS) -O3 -o  $(BIN)seqyclean $(OBJ)main.o $(OBJ)flash.o $(OBJ)dup.o $(OBJ)Roche.o $(OBJ)Illumina.o $(OBJ)MainPipeLine.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)pairwise.o $(OBJ)iz_SSAHA.o $(OBJ)Report.o $(OBJ)QualTrim.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sffreader.o $(OBJ)sff.o $(OBJ)poly.o $(OBJ)abi.o $(OBJ)gzstream.o -I$(LIBRARY) -lpthread -lz
	
main.o :  
	$(CXX) -Wall -g -O3  -c -o $(OBJ)main.o $(SRC)main.cpp 
		
flash.o : 
	$(CXX) -Wall -g -O3 -c -o $(OBJ)flash.o $(SRC)flash.cpp 

dup.o : 
	$(CXX) -Wall -g -O3  -c -o $(OBJ)dup.o $(SRC)dup.cpp 
	
Roche.o : 
	$(CXX) -Wall -g -O3  -c -o $(OBJ)Roche.o $(SRC)Roche.cpp 
	
Illumina.o : 
	$(CXX) -Wall -g -O3 -c -o $(OBJ)Illumina.o $(SRC)Illumina.cpp 

MainPipeLine.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)MainPipeLine.o $(SRC)MainPipeLine.cpp 
	
KMerRoutine.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)KMerRoutine.o $(SRC)KMerRoutine.cpp
	
Dictionary.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)Dictionary.o $(SRC)Dictionary.cpp

pairwise.o :  
	$(CXX) -Wall -g -O3  -c -o $(OBJ)pairwise.o $(SRC)pairwise.cpp
	
iz_SSAHA.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)iz_SSAHA.o $(SRC)iz_SSAHA.cpp
	
Report.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)Report.o $(SRC)Report.cpp
	
	
QualTrim.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)QualTrim.o $(SRC)QualTrim.cpp
	
Read.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)Read.o $(SRC)Read.cpp
	
util.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)util.o $(SRC)util.cpp
	
ascii.o :
	$(CXX) -Wall -g -O3 -c -o $(OBJ)ascii.o $(SRC)ascii.cpp
	
sffreader.o: $(SRC)sffreader.cpp $(SRC)sff.h
	g++ -g -I $(SRC) -O3 -c -o $(OBJ)sffreader.o $(SRC)sffreader.cpp
	
sff.o: 
	g++ -g -O3 -c -o $(OBJ)sff.o $(SRC)sff.c
	
poly.o: 
	g++ -g -O3 -c -o $(OBJ)poly.o $(SRC)poly.c
	
abi.o: 
	g++ -g -O3 -c -o $(OBJ)abi.o $(SRC)abi.c
	
libgzstream.a : $(OBJ)gzstream.o $(SRC)gzstream.h
	${AR} $(OBJ)libgzstream.a $(OBJ)gzstream.o
	
gzstream.o : $(SRC)gzstream.C $(SRC)gzstream.h
	#$(CXX) -I $(SRC) -O -Wall -c -o $(OBJ)gzstream.o $(SRC)gzstream.C 
	gcc -I $(SRC) -O -Wall -c -o $(OBJ)gzstream.o $(SRC)gzstream.C 

mkobj :
	rm -rf ${OBJ}
	mkdir ${OBJ}

mkbin :
	rm -rf ${BIN}
	mkdir ${BIN}
	

clean:
	rm -rf ${BIN}
	rm -rf ${OBJ}
	
	
test :  
	@echo "Test 454..."
	@bin/./seqyclean -454 test_data/in.sff -qual -o unit_test/test454 > /dev/null
	@diff test_data/test454.sff unit_test/test454.sff > /dev/null || echo "Test 454 failed"
	@echo "Test PE Illumina..."
	@bin/./seqyclean -1 test_data/artif_pe1.fastq -2 test_data/artif_pe2.fastq -qual -o unit_test/testIlluminaPE > /dev/null
	@diff test_data/testIlluminaPE_PE1.fastq unit_test/testIlluminaPE_PE1.fastq > /dev/null || echo "Test PE Illumina failed"
	@diff test_data/testIlluminaPE_PE2.fastq unit_test/testIlluminaPE_PE2.fastq > /dev/null || echo "Test PE Illumina failed"
	@diff test_data/testIlluminaPE_SE.fastq unit_test/testIlluminaPE_SE.fastq > /dev/null || echo "Test PE Illumina failed"
	@echo "Test SE Illumina..."
	@bin/./seqyclean -U test_data/artif_pe1.fastq -qual -o unit_test/testIlluminaSE > /dev/null
	@diff test_data/testIlluminaSE_SE.fastq unit_test/testIlluminaSE_SE.fastq > /dev/null || echo "Test SE Illumina failed" 
	@echo "Test removing of duplicates..."
	@bin/./seqyclean --dup -1 test_data/R1.fastq -2 test_data/R2.fastq -o unit_test/testIlluminaDup > /dev/null
	@diff test_data/testIlluminaDup_PE1.fastq unit_test/testIlluminaDup_PE1.fastq > /dev/null || echo "Test removing of duplicates failed" 
	@diff test_data/testIlluminaDup_PE2.fastq unit_test/testIlluminaDup_PE2.fastq > /dev/null || echo "Test removing of duplicates failed" 
	@diff test_data/testIlluminaDup_SE.fastq unit_test/testIlluminaDup_SE.fastq > /dev/null || echo "Test removing of duplicates failed"
	@echo "Test Illumina poly A/T PE..."
	@bin/./seqyclean -1 test_data/artif_pat_pe1.fastq -2 test_data/artif_pat_pe2.fastq -qual -polyat -o unit_test/test_polAT_IlluminaPE > /dev/null
	@diff test_data/test_polAT_IlluminaPE_PE1.fastq unit_test/test_polAT_IlluminaPE_PE1.fastq > /dev/null || echo "Test poly A/T Illumina PE failed"
	@diff test_data/test_polAT_IlluminaPE_PE2.fastq unit_test/test_polAT_IlluminaPE_PE2.fastq > /dev/null || echo "Test poly A/T Illumina PE failed"
	@diff test_data/test_polAT_IlluminaPE_SE.fastq unit_test/test_polAT_IlluminaPE_SE.fastq > /dev/null || echo "Test poly A/T Illumina PE failed"
	@echo "Test Illumina poly A/T SE..."
	@bin/./seqyclean -U test_data/artif_pat_pe1.fastq -qual -polyat -o unit_test/test_polAT_IlluminaSE > /dev/null
	@diff test_data/test_polAT_IlluminaSE_SE.fastq unit_test/test_polAT_IlluminaSE_SE.fastq > /dev/null || echo "Test poly A/T Illumina SE failed"
	@echo "Done."
	@rm -r unit_test