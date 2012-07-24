ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib

CXX           = g++ -m64
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++ -m64
LDFLAGS       = -g
SOFLAGS       = -shared

#PG da qui per macosx
#PG -----------------

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

NGLIBS         = $(ROOTGLIBS) 
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
CXXFLAGS      += $(FASTJETFLAGS)
LIBS           = $(ROOTLIBS)

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        +=  -lTreePlayer
NGLIBS        += $(FASTJETLIBS)
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
INCLUDEDIRCOMMON = ../
SRCDIR           = ./
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = ./lib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/
OUTLIBEGAMMA	 = $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

all:  lumiproc

lumiproc:  $(SRCDIR)lumiproc.C 
	$(CXX) $(CXXFLAGS) -o lumiproc  $(GLIBS) $ $<


VecbosApp.clean:
	rm -f VecbosApp

clean:
	rm -f $(OUTLIB)*.o 

