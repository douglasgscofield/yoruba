BAMTOOLS_ROOT = ../bamtools
BAMTOOLS_BUILD_DIR = $(BAMTOOLS_ROOT)/build
BAMTOOLS_INCLUDE_DIR = $(BAMTOOLS_ROOT)/include
BAMTOOLS_LIB_DIR = $(BAMTOOLS_ROOT)/lib

CXX=		g++
CXXFLAGS=	-Wall -pg -g -D_WITH_DEBUG -D_BAMTOOLS_EXTENSION -I$(BAMTOOLS_INCLUDE_DIR)
#CXXFLAGS=	-Wall -g -D_WITH_DEBUG -D_BAMTOOLS_EXTENSION -I$(BAMTOOLS_INCLUDE_DIR)

PROG=		yoruba

LIBS=		-lbamtools -lz

OBJS=		yoruba.o \
			yoruba_gbagbe.o \
			yoruba_inu.o \
			yoruba_kojopodipo.o \
			yoruba_seda.o \
			yoruba_util.o

HEAD_COMM=  yoruba_util.h SimpleOpt.h

HEAD=		$(HEAD_COMM) \
			yoruba.h \
			yoruba_gbagbe.h \
			yoruba_inu.h \
			yoruba_kojopodipo.h \
			yoruba_seda.h


#---------------------------  Main program


all: $(PROG)

yoruba: bamtools-headers $(OBJS) bamtools-static-library
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -L$(BAMTOOLS_LIB_DIR) $(LIBS)


#---------------------------  Individual object files


# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

# rebuild everything if the common headers change
$(OBJS): $(HEAD_COMM)

# rebuild the main file if any header changes
yoruba.o: $(HEAD)

yoruba_gbagbe.o: yoruba_gbagbe.h 

yoruba_inu.o: yoruba_inu.h 

yoruba_kojopodipo.o: yoruba_kojopodipo.h 

# seda (mark/remove duplicates) is not yet read for alpha
yoruba_seda.o: yoruba_seda.h 

yoruba_util.o: yoruba_util.h

yoruba_ibeji.o: ibejiAlignment.h processReadPair.h 


#---------------------------  Other targets


bamtools-headers:
	( cd $(BAMTOOLS_BUILD_DIR) ; make SharedHeaders )
	( cd $(BAMTOOLS_BUILD_DIR) ; make APIHeaders )

bamtools-static-library:
	( cd $(BAMTOOLS_BUILD_DIR) ; make BamTools-static/fast )

bamtools-clean:
	( cd $(BAMTOOLS_BUILD_DIR) ; make clean )

clean:
	rm -f gmon.out *.o $(PROG)

clean-all: clean bamtools-clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


sefibo.o: ibejiAlignment.h yoruba_util.h SimpleOpt.h

sefibo: sefibo.o processReadPair.o yoruba_util.o ibejiAlignment.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

ibejiAlignment.o: ibejiAlignment.h

processReadPair.o: processReadPair.h util.h

