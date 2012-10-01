BAMTOOLS_ROOT = ../bamtools
BAMTOOLS_INCLUDE_DIR = $(BAMTOOLS_ROOT)/include
BAMTOOLS_LIB_DIR = $(BAMTOOLS_ROOT)/lib

CXX=		g++
# CXXFLAGS=	-Wall -O3
# CXXFLAGS=	-Wall -g -D_WITH_DEBUG -D_STANDALONE -I$(BAMTOOLS_INCLUDE_DIR)
# CXXFLAGS=	-Wall --pedantic-errors -g -D_WITH_DEBUG -D_BAMTOOLS_EXTENSION -I$(BAMTOOLS_INCLUDE_DIR)
CXXFLAGS=	-Wall -g -D_WITH_DEBUG -D_BAMTOOLS_EXTENSION -I$(BAMTOOLS_INCLUDE_DIR)
PROG=		yoruba
LIBS=		-lbamtools -lz

all: $(PROG)

yoruba: yoruba.o yoruba_gbagbe.o yoruba_inu.o yoruba_kojopodipo.o yoruba_util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

yoruba.o: yoruba.h yoruba_kojopodipo.h

ibejiAlignment.o: ibejiAlignment.h

processReadPair.o: processReadPair.h util.h

# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

yoruba_ibeji.o: ibejiAlignment.h processReadPair.h util.h SimpleOpt.h

yoruba_gbagbe.o: yoruba_gbagbe.h yoruba_util.h SimpleOpt.h

yoruba_inu.o: yoruba_inu.h yoruba_util.h SimpleOpt.h

yoruba_kojopodipo.o: yoruba_kojopodipo.h yoruba_util.h SimpleOpt.h

yoruba_util.o: yoruba_util.h

#yoruba_gbagbe: yoruba_gbagbe.c yoruba_util.o yoruba_inu.h SimpleOpt.h
#	$(CXX) $(CXXFLAGS) -c -D_STANDALONE yoruba_gbagbe.c
#	$(CXX) $(CXXFLAGS) -o $@ yoruba_gbagbe.o -L$(BAMTOOLS_LIB_DIR) $(LIBS)

#yoruba_inu: yoruba_inu.c yoruba_util.o yoruba_inu.h SimpleOpt.h
#	$(CXX) $(CXXFLAGS) -c -D_STANDALONE yoruba_inu.c
#	$(CXX) $(CXXFLAGS) -o $@ yoruba_inu.o -L$(BAMTOOLS_LIB_DIR) $(LIBS)

#yoruba_kojopodipo: yoruba_kojopodipo.c yoruba_util.o yoruba_kojopodipo.h SimpleOpt.h
#	$(CXX) $(CXXFLAGS) -c -D_STANDALONE yoruba_kojopodopo.c
#	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

sefibo.o: ibejiAlignment.h yoruba_util.h SimpleOpt.h

sefibo: sefibo.o processReadPair.o yoruba_util.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

clean:
	rm -fr gmon.out *.o *.a $(PROG) yoruba_kojopodipo a.out *~

