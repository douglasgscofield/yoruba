BAMTOOLS_ROOT = ../bamtools
BAMTOOLS_INCLUDE_DIR = $(BAMTOOLS_ROOT)/include
BAMTOOLS_LIB_DIR = $(BAMTOOLS_ROOT)/lib

CXX=		g++
# CXXFLAGS=	-Wall -O3
# CXXFLAGS=	-Wall -g -D_WITH_DEBUG -D_STANDALONE -I$(BAMTOOLS_INCLUDE_DIR)
CXXFLAGS=	-Wall -g -D_WITH_DEBUG -I$(BAMTOOLS_INCLUDE_DIR)
# PROG=		ibeji sefibo kojopodipo
# PROG=		yoruba_kojopodipo
PROG=		yoruba
LIBS=		-lbamtools -lz

all: $(PROG)

yoruba: yoruba.o yoruba_kojopodipo.o yoruba_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

yoruba.o: yoruba.h yoruba_kojopodipo.h

ibejiAlignment.o: ibejiAlignment.h

yoruba_utils.o: yoruba_utils.h

processReadPair.o: processReadPair.h utils.h

# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

ibeji.o: ibejiAlignment.h processReadPair.h utils.h SimpleOpt.h

yoruba_kojopodipo.o: yoruba_kojopodipo.h yoruba_utils.h SimpleOpt.h

yoruba_kojopodipo: yoruba_kojopodipo.o 
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

ibeji: ibeji.o processReadPair.o utils.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

sefibo.o: ibejiAlignment.h utils.h SimpleOpt.h

sefibo: sefibo.o processReadPair.o utils.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

clean:
	rm -fr gmon.out *.o *.a $(PROG) a.out *~
