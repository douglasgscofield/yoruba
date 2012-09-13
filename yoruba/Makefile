BAMTOOLS_ROOT = /Users/douglasgscofield/tools/bamtools
BAMTOOLS_INCLUDE_DIR = $(BAMTOOLS_ROOT)/include
BAMTOOLS_LIB_DIR = $(BAMTOOLS_ROOT)/lib

CXX=		g++
# CXXFLAGS=	-Wall -O3
CXXFLAGS=	-Wall -g -I$(BAMTOOLS_INCLUDE_DIR)
PROG=		ibeji sefibo
LIBS=		-lz -lbamtools

all: $(PROG)

ibejiAlignment.o: ibejiAlignment.h

utils.o: utils.h

processReadPair.o: processReadPair.h utils.h

# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

ibeji.o: ibejiAlignment.h processReadPair.h utils.h SimpleOpt.h

ibeji: ibeji.o processReadPair.o utils.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ ibeji.o processReadPair.o utils.o ibejiAlignment.o -L$(BAMTOOLS_LIB_DIR) $(LIBS)

sefibo.o: ibejiAlignment.h utils.h SimpleOpt.h

sefibo: sefibo.o processReadPair.o utils.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ sefibo.o processReadPair.o utils.o ibejiAlignment.o -L$(BAMTOOLS_LIB_DIR) $(LIBS)

clean:
	rm -fr gmon.out *.o *.a $(PROG) a.out *~
