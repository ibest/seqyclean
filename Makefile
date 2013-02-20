CXX = g++
CFLAGS = -Wall -g -L. -${OBJ}lgzstream  
AR       = ar cr
SRC = src/
BIN = bin/
OBJ = obj/

all:  mkobj mkbin gzstream.o libgzstream.a abi.o poly.o sff.o sffreader.o ascii.o util.o Read.o QualTrim.o Report.o iz_SSAHA.o pairwise.o Dictionary.o KMerRoutine.o MainPipeLine.o Illumina.o main.o seqyclean 
		
#					
seqyclean :   $(OBJ)main.o $(OBJ)Illumina.o $(OBJ)MainPipeLine.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)pairwise.o $(OBJ)iz_SSAHA.o $(OBJ)Report.o $(OBJ)QualTrim.o $(OBJ)sffreader.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sffreader.o $(OBJ)sff.o $(OBJ)abi.o $(OBJ)gzstream.o
	$(CXX) $(CFLAGS)  -o  $(BIN)seqyclean $(OBJ)main.o $(OBJ)Illumina.o $(OBJ)MainPipeLine.o $(OBJ)KMerRoutine.o $(OBJ)Dictionary.o $(OBJ)pairwise.o $(OBJ)iz_SSAHA.o $(OBJ)Report.o $(OBJ)QualTrim.o $(OBJ)Read.o $(OBJ)util.o $(OBJ)ascii.o $(OBJ)sffreader.o $(OBJ)sff.o $(OBJ)poly.o $(OBJ)abi.o $(OBJ)gzstream.o -lpthread -Xlinker  -zmuldefs -lz
	
main.o :  
	$(CXX) -Wall -g  -c -o $(OBJ)main.o $(SRC)main.cpp 
	
Illumina.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Illumina.o $(SRC)Illumina.cpp 
	
MainPipeLine.o :
	$(CXX) -Wall -g  -c -o $(OBJ)MainPipeLine.o $(SRC)MainPipeLine.cpp 
	
KMerRoutine.o :
	$(CXX) -Wall -g  -c -o $(OBJ)KMerRoutine.o $(SRC)KMerRoutine.cpp
	
Dictionary.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Dictionary.o $(SRC)Dictionary.cpp

pairwise.o :  
	$(CXX) -Wall -g  -c -o $(OBJ)pairwise.o $(SRC)pairwise.cpp
	
iz_SSAHA.o :
	$(CXX) -Wall -g  -c -o $(OBJ)iz_SSAHA.o $(SRC)iz_SSAHA.cpp
	
Report.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Report.o $(SRC)Report.cpp
	
	
QualTrim.o :
	$(CXX) -Wall -g  -c -o $(OBJ)QualTrim.o $(SRC)QualTrim.cpp
	
Read.o :
	$(CXX) -Wall -g  -c -o $(OBJ)Read.o $(SRC)Read.cpp
	
util.o :
	$(CXX) -Wall -g  -c -o $(OBJ)util.o $(SRC)util.cpp
	
ascii.o :
	$(CXX) -Wall -g  -c -o $(OBJ)ascii.o $(SRC)ascii.cpp
	
sffreader.o: $(SRC)sffreader.cpp $(SRC)sff.h
	g++ -g -I $(SRC) -c -o $(OBJ)sffreader.o $(SRC)sffreader.cpp
	
sff.o: $(SRC)sff.h $(SRC)sff.c
	g++ -g -I $(SRC) -c -o $(OBJ)sff.o $(SRC)sff.c
	
poly.o: 
	g++ -g -I $(SRC) -c -o $(OBJ)poly.o $(SRC)poly.c
	
abi.o: $(SRC)abi.c $(SRC)abi.h
	g++ -g -I $(SRC) -c -o $(OBJ)abi.o $(SRC)abi.c
	
libgzstream.a : $(OBJ)gzstream.o $(SRC)gzstream.h
	${AR} $(OBJ)libgzstream.a $(OBJ)gzstream.o
	
gzstream.o : $(SRC)gzstream.C $(SRC)gzstream.h
	$(CXX) -I $(SRC) -O -Wall -c -o $(OBJ)gzstream.o $(SRC)gzstream.C 

mkobj :
	rm -rf ${OBJ}
	mkdir ${OBJ}

mkbin :
	rm -rf ${BIN}
	mkdir ${BIN}
	

clean:
	rm -rf ${BIN}
	rm -rf ${OBJ}
