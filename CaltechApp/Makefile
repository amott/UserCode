ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lMLP -lXMLIO -lTMVA
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lTMVA -lMLP -lXMLIO

FASTJETFLAGS = $(shell FASTJET/bin/fastjet-config --cxxflags)
FASTJETLIBS  = $(shell FASTJET/bin/fastjet-config --libs --plugins)

CXX           = g++ -m64
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs 
LD            = g++ -m64
LDFLAGS       = -g
SOFLAGS       = -shared

#PG da qui per macosx
#PG -----------------

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
CXXFLAGS      += $(FASTJETFLAGS)
LIBS           = $(ROOTLIBS)
NGLIBS         = $(ROOTGLIBS) 
#NGLIBS        += -lMinuit -lTMVA.1 -lMLP -lTreePlayer
NGLIBS        += -lMinuit -lMLP -lTreePlayer
NGLIBS        += $(FASTJETLIBS)
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./include
INCLUDEDIRCOMMON = ./
#INCLUDEDIRTMVA   = ./TMVA/include
SRCDIR           = ./src/
#CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I$(INCLUDEDIRTMVA) -I.
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = ./lib/
H2GLIB           = ./h2glib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/
OUTLIBEGAMMA	 = $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

all:  lib CaltechReducerApp

lib: 	$(OUTLIB)CaltechReducer.o \
	$(OUTLIB)HggVertexing.o \
	$(OUTLIB)HggEGEnergyCorrector.o \
	VecbosObjects \
	commonlib \
	egammalib

.PHONE: lib all 
# analysis functions
CaltechReducerApp: $(SRCDIR)CaltechReducerApp.C \
		$(OUTLIB)CaltechReducer.o \
		$(OUTLIB)ArgParser.o
		$(CXX) $(CXXFLAGS) -o CaltechReducerApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*.o $(OUTLIBEGAMMA)/*.o  $(GLIBS) $ $<	


$(OUTLIB)HggEGEnergyCorrector.o: $(SRCDIR)HggEGEnergyCorrector.cc \
				$(OUTLIB)Vecbos.o \
				$(OUTLIB)GBRTree.o \
				$(OUTLIB)GBRForest.o\
				$(OUTLIB)MitDict.o \
				$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggEGEnergyCorrector.o $<

$(OUTLIB)HggVertexing.o: $(SRCDIR)HggVertexingNew.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggVertexing.o $<

$(OUTLIB)HggEnergyScale.o: $(SRCDIR)HggEnergyScale.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)MitDict.o \
			$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggEnergyScale.o $<

$(OUTLIB)CaltechReducer.o: $(SRCDIR)CaltechReducer.cc \
		$(OUTLIB)HggEGEnergyCorrector.o \
		$(OUTLIB)HggVertexing.o \
		VecbosObjects \
		$(OUTLIB)VecbosDict.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CaltechReducer.o $<


## -------------------------
## VECBOS OBJECTS
## -------------------------
$(OUTLIB)VecbosBaseObject.o: $(SRCDIR)VecbosBaseObject.cc
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosBaseObject.o $<

$(OUTLIB)VecbosGen.o: $(SRCDIR)VecbosGen.cc \
		$(OUTLIB)VecbosBaseObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosGen.o $<

$(OUTLIB)VecbosPhysicsObject.o: $(SRCDIR)VecbosPhysicsObject.cc \
		$(OUTLIB)VecbosGen.o \
		$(OUTLIB)VecbosBaseObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosPhysicsObject.o $<

$(OUTLIB)VecbosVtx.o: $(SRCDIR)VecbosVtx.cc \
		$(OUTLIB)VecbosBaseObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosVtx.o $<

$(OUTLIB)VecbosTrack.o: $(SRCDIR)VecbosTrack.cc \
		$(OUTLIB)VecbosBaseObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosTrack.o $<

$(OUTLIB)VecbosJet.o: $(SRCDIR)VecbosJet.cc \
		$(OUTLIB)VecbosPhysicsObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosJet.o $<

$(OUTLIB)VecbosMET.o: $(SRCDIR)VecbosMET.cc \
		$(OUTLIB)VecbosPhysicsObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosMET.o $<

$(OUTLIB)VecbosMu.o: $(SRCDIR)VecbosMu.cc \
		$(OUTLIB)VecbosPhysicsObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosMu.o $<


$(OUTLIB)VecbosEGObject.o: $(SRCDIR)VecbosEGObject.cc \
		$(OUTLIB)VecbosPhysicsObject.o
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosEGObject.o $<

VecbosObjects:  $(OUTLIB)VecbosEGObject.o \
		$(OUTLIB)VecbosMu.o \
		$(OUTLIB)VecbosMET.o \
		$(OUTLIB)VecbosJet.o \
		$(OUTLIB)VecbosVtx.o \
		$(OUTLIB)VecbosTrack.o \
		$(OUTLIB)VecbosPhysicsObject.o \
		$(OUTLIB)VecbosGen.o \
		$(OUTLIB)VecbosBaseObject.o

.PHONY: VecbosObjects

VecbosDict.cc: 	$(INCLUDEDIR)/VecbosObjectLinkDef.h \
		VecbosObjects
	rootcint -l -f VecbosDict.cc -c -I$(INCLUDEDIR) -p $(INCLUDEDIR)/VecbosEGObject.hh $(INCLUDEDIR)/VecbosGen.hh $(INCLUDEDIR)/VecbosVtx.hh $(INCLUDEDIR)/VecbosJet.hh $(INCLUDEDIR)/VecbosMET.hh $(INCLUDEDIR)/VecbosMu.hh  \
	$(INCLUDEDIR)/VecbosVtx.hh $(INCLUDEDIR)/VecbosTrack.hh $<

$(OUTLIB)VecbosDict.o: VecbosDict.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosDict.o $<

$(OUTLIB)GBRForest.o: $(SRCDIR)GBRForest.cxx $(OUTLIB)GBRTree.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GBRForest.o $<

$(OUTLIB)GBRTree.o: $(SRCDIR)GBRTree.cxx
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GBRTree.o $<

MitPhysicsDict.cc: $(INCLUDEDIR)/MitPhysicsUtilsLinkDef.h $(OUTLIB)GBRForest.o $(OUTLIB)GBRTree.o
	rootcint -l -f MitPhysicsDict.cc -c -I$(INCLUDEDIR) -p $(INCLUDEDIR)/GBRForest.h $(INCLUDEDIR)/GBRTree.h MitPhysicsUtilsLinkDef.h

$(OUTLIB)MitDict.o: MitPhysicsDict.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MitDict.o $<

$(OUTLIB)ArgParser.o: $(SRCDIR)ArgParser.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ArgParser.o $<


$(OUTLIB)JetCorrectorParameters.o: $(SRCDIR)JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectorParameters.o $<

$(OUTLIB)SimpleJetCorrectionUncertainty.o: $(SRCDIR)SimpleJetCorrectionUncertainty.cc \
	$(OUTLIB)JetCorrectorParameters.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrectionUncertainty.o $<

$(OUTLIB)JetCorrectionUncertainty.o: $(SRCDIR)JetCorrectionUncertainty.cc \
	$(OUTLIB)SimpleJetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectionUncertainty.o $<

$(OUTLIB)Jet.o: $(SRCDIR)Jet.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Jet.o $<

$(OUTLIB)MET.o: $(SRCDIR)MET.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MET.o $<

$(OUTLIB)CaloTower.o: $(SRCDIR)CaloTower.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CaloTower.o $<

# auxiliary functions to compute selections/efficiencies
$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Conditions.o $<
$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Utils.o $<
$(OUTLIBCOMMON)Skimmer.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Skimmer.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Skimmer.o $<
$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Counters.o $<
$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Selection.o $<
$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)TriggerMask.o $<
$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<
$(OUTLIBCOMMON)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CutBasedEleIDSelector.o $<
$(OUTLIBCOMMON)EcalCleaner.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/EcalCleaner.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EcalCleaner.o $<
$(OUTLIBEGAMMA)ElectronTrackerIsolation.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronTrackerIsolation.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronTrackerIsolation.o $<
$(OUTLIBEGAMMA)ElectronCaloIsolation.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronCaloIsolation.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronCaloIsolation.o $<
$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronBestCandidateSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronBestCandidateSelector.o $<
$(OUTLIBEGAMMA)LikelihoodPdf.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdf.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodPdf.o $<
$(OUTLIBEGAMMA)LikelihoodSpecies.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodSpecies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodSpecies.o $<
$(OUTLIBEGAMMA)LikelihoodPdfProduct.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdfProduct.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodPdfProduct.o $<
$(OUTLIBEGAMMA)ElectronLikelihood.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronLikelihood.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronLikelihood.o $<

commonlib: 	$(OUTLIBCOMMON)Conditions.o \
		$(OUTLIBCOMMON)Utils.o \
		$(OUTLIBCOMMON)Skimmer.o \
		$(OUTLIBCOMMON)Counters.o \
		$(OUTLIBCOMMON)Selection.o \
		$(OUTLIBCOMMON)TriggerMask.o \
		$(OUTLIBCOMMON)EfficiencyEvaluator.o \
		$(OUTLIBCOMMON)CutBasedEleIDSelector.o \
		$(OUTLIBCOMMON)EcalCleaner.o

egammalib:	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
		$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
		$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
		$(OUTLIBEGAMMA)LikelihoodPdf.o \
		$(OUTLIBEGAMMA)LikelihoodSpecies.o \
		$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
		$(OUTLIBEGAMMA)ElectronLikelihood.o

.PHONY: egammalib commonlib
$(OUTLIB)VecbosBase.o: $(SRCDIR)VecbosBase.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosBase.o $<


$(OUTLIB)Vecbos.o: $(SRCDIR)Vecbos.cc \
		$(OUTLIB)VecbosBase.o \
		$(OUTLIB)CaloTower.o \
		$(OUTLIB)Jet.o \
		$(OUTLIB)JetCorrectionUncertainty.o \
		$(OUTLIB)MET.o \
		$(OUTLIBCOMMON)Conditions.o \
		$(OUTLIBCOMMON)Utils.o \
		$(OUTLIBCOMMON)Skimmer.o \
		$(OUTLIBCOMMON)Counters.o \
		$(OUTLIBCOMMON)Selection.o \
		$(OUTLIBCOMMON)TriggerMask.o \
		$(OUTLIBCOMMON)EfficiencyEvaluator.o \
		$(OUTLIBCOMMON)CutBasedEleIDSelector.o \
		$(OUTLIBCOMMON)EcalCleaner.o \
		$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
		$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
		$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
		$(OUTLIBEGAMMA)LikelihoodPdf.o \
		$(OUTLIBEGAMMA)LikelihoodSpecies.o \
		$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
		$(OUTLIBEGAMMA)ElectronLikelihood.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Vecbos.o $<
VecbosApp.clean:
	rm -f VecbosApp

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o $(OUTLIBEGAMMA)*.o
	rm -f VecbosDict.cc
	rm -f VecbosDict.h
	rm -f MitPhysicsDict.cc
	rm -f MitPhysicsDict.h
	rm -f VecbosApp
	rm -f CaltechReducerApp
