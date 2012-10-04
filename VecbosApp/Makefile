ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) 

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
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/
OUTLIBEGAMMA	 = $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

all:  lib VecbosApp CountEvents GenVecbosApp

AlpgenValidation:  $(SRCDIR)AlpgenValidation.C \
	$(OUTLIB)Vecbos.o \
	$(OUTLIB)ThiagoAnalysis.o
	$(CXX) $(CXXFLAGS) -o AlpgenValidation $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<

CountEvents:  $(SRCDIR)CountEvents.C	
	$(CXX) $(CXXFLAGS) -o CountEvents $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<

OptimizeElectronID: $(SRCDIR)OptimizeElectronID.C
	$(CXX) $(CXXFLAGS) -o OptimizeElectronID $(GLIBS) $ $<

OptimizeElectronIsolation: $(SRCDIR)OptimizeElectronIsolation.C
	$(CXX) $(CXXFLAGS) -o OptimizeElectronIsolation $(GLIBS) $ $<

AnalyzeTTbarMC: $(SRCDIR)AnalyzeTTbarMC.C \
	$(OUTLIB)VecbosTTCounter.o
	$(CXX) $(CXXFLAGS) -o AnalyzeTTbarMC $(OUTLIB)/VecbosBase.o $(OUTLIB)/VecbosTTCounter.o $(GLIBS) $ $<

GenVecbosApp: $(SRCDIR)GenVecbosApp.C \
	$(OUTLIB)VecbosBase.o \
	$(OUTLIB)Vecbos.o \
	$(OUTLIB)GenVecbos.o \
	$(OUTLIB)GenWjets.o \
	$(OUTLIB)GenZjets.o 
	$(CXX) $(CXXFLAGS) -o GenVecbosApp $(OUTLIB)/*.o $(GLIBS) $(OUTLIBEGAMMA)/*o $(OUTLIBCOMMON)/*o $ $<

VecbosApp:  $(SRCDIR)VecbosApp.C \
	$(OUTLIB)Vecbos.o \
	$(OUTLIB)DiJet.o \
	$(OUTLIB)Razor.o \
	$(OUTLIB)RazorMultiB.o \
	$(OUTLIB)RazorDiMuB.o \
	$(OUTLIB)MonoJet.o \
	$(OUTLIB)CandleCalib.o \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)RedVecbosTree.o \
	$(OUTLIB)RedVecbosPFTree.o \
	$(OUTLIB)RedTopTree.o \
	$(OUTLIB)RedEleIdIsolPFTree.o \
	$(OUTLIB)RedVecbosMcTree.o \
	$(OUTLIB)RedVecbosVertexTree.o \
	$(OUTLIB)RedIsolationVtxOptimTree.o \
	$(OUTLIB)RedEleIDOptimTree.o \
	$(OUTLIB)LQ3Helper.o \
	$(OUTLIB)VBTFLeptEff.o \
	$(OUTLIB)LQ3Analysis.o \
	$(OUTLIB)SUSYNLO.o 	
	$(CXX) $(CXXFLAGS) -o VecbosApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<
# fix for 2012 data format
#	$(OUTLIB)SUSYMultiTop.o \
#	$(OUTLIB)RazorDiPhoton.o \
#	$(OUTLIB)RazorBoostedTop.o \
#	$(OUTLIB)SUSYMultiB.o \
#	$(OUTLIB)CreateWJetDataset.o \

FindDuplicatedApp:  $(SRCDIR)FindDuplicatedApp.C \
	$(OUTLIB)FindDuplicatedEvents.o 
	$(CXX) $(CXXFLAGS) -o FindDuplicatedApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<

lib: 	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIBCOMMON)CutBasedEleIDSelector.o \
	$(OUTLIBCOMMON)EcalCleaner.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIB)RedVecbosTree.o \
	$(OUTLIB)RedVecbosPFTree.o \
	$(OUTLIB)RedTopTree.o \
	$(OUTLIB)RedEleIdIsolPFTree.o \
	$(OUTLIB)RedVecbosMcTree.o \
	$(OUTLIB)RedVecbosVertexTree.o \
	$(OUTLIB)RedIsolationVtxOptimTree.o \
	$(OUTLIB)RedEleIDOptimTree.o \
	$(OUTLIB)VecbosTTCounter.o \
	$(OUTLIB)TestAnalysis.o \
	$(OUTLIB)DiJet.o \
	$(OUTLIB)Razor.o \
	$(OUTLIB)RazorLeptons.o \
	$(OUTLIB)RazorMultiB.o \
	$(OUTLIB)RazorDiMuB.o \
	$(OUTLIB)MonoJet.o \
	$(OUTLIB)CandleCalib.o \
	$(OUTLIB)LQ3Helper.o \
	$(OUTLIB)VBTFLeptEff.o \
	$(OUTLIB)SUSYNLO.o 	
#	$(OUTLIB)VecbosMuMuSelection.o \
#	$(OUTLIB)LFJetControlSample.o \
#	$(OUTLIB)CandleCalib_ee.o 
#	$(OUTLIB)EventShapeAnalysis.o \
#	$(OUTLIB)VtxStudy.o \
#	$(OUTLIB)GenJetAnalysis.o \
#	$(OUTLIB)RSAnalysis.o \
#	$(OUTLIB)RSZZAnalysis.o \
#	$(OUTLIB)SUSYAnalysis.o \
#	$(OUTLIB)ThiagoAnalysis.o \
#	$(OUTLIB)VecbosAnalysis.o   \
# fix for 2012 data format
#	$(OUTLIB)SUSYMultiTop.o \
#	$(OUTLIB)RazorDiPhoton.o \
#	$(OUTLIB)RazorBoostedTop.o \
#	$(OUTLIB)SUSYMultiB.o \
#	$(OUTLIB)VecbosEESelection.o \
#	$(OUTLIB)TopControlSample.o \
#	$(OUTLIB)VecbosPFEESelection.o \
#	$(OUTLIB)CreateWJetDataset.o \
#	$(OUTLIB)LQ3Analysis.o \

# analysis functions
$(OUTLIB)EventShapeAnalysis.o: $(SRCDIR)EventShapeAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)EventShapeAnalysis.o $<

$(OUTLIB)GenJetAnalysis.o: $(SRCDIR)GenJetAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenJetAnalysis.o $<

$(OUTLIB)VtxStudy.o: $(SRCDIR)VtxStudy.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VtxStudy.o $<

$(OUTLIB)FindDuplicatedEvents.o: $(SRCDIR)FindDuplicatedEvents.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FindDuplicatedEvents.o $<

$(OUTLIB)SUSYAnalysis.o: $(SRCDIR)SUSYAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SUSYAnalysis.o $<

$(OUTLIB)CandleCalib.o: $(SRCDIR)CandleCalib.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CandleCalib.o $<

$(OUTLIB)GammaPlusJet.o: $(SRCDIR)GammaPlusJet.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GammaPlusJet.o $<

$(OUTLIB)CandleCalib_ee.o: $(SRCDIR)CandleCalib_ee.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CandleCalib_ee.o $<

$(OUTLIB)DiJet.o: $(SRCDIR)DiJet.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)DiJet.o $<

$(OUTLIB)RazorHiggsBB.o: $(SRCDIR)RazorHiggsBB.cc $(OUTLIB)Vecbos.o $(OUTLIB)SUSYNLO.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorHiggsBB.o $<

$(OUTLIB)RazorBoostedTop.o: $(SRCDIR)RazorBoostedTop.cc $(OUTLIB)Vecbos.o $(OUTLIB)SUSYNLO.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorBoostedTop.o $<

$(OUTLIB)RazorDiPhoton.o: $(SRCDIR)RazorDiPhoton.cc $(OUTLIB)Vecbos.o $(OUTLIB)SUSYNLO.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorDiPhoton.o $<

$(OUTLIB)Razor.o: $(SRCDIR)Razor.cc $(OUTLIB)Vecbos.o $(OUTLIB)SUSYNLO.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Razor.o $<

$(OUTLIB)RazorLeptons.o: $(SRCDIR)RazorLeptons.cc $(OUTLIB)Vecbos.o $(OUTLIB)SUSYNLO.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorLeptons.o $<

$(OUTLIB)SUSYTau.o: $(SRCDIR)SUSYTau.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SUSYTau.o $<

$(OUTLIB)SUSYMultiTop.o: $(SRCDIR)SUSYMultiTop.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SUSYMultiTop.o $<

$(OUTLIB)SUSYMultiB.o: $(SRCDIR)SUSYMultiB.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SUSYMultiB.o $<

$(OUTLIB)RazorMultiB.o: $(SRCDIR)RazorMultiB.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorMultiB.o $<

$(OUTLIB)RazorDiMuB.o: $(SRCDIR)RazorDiMuB.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RazorDiMuB.o $<

$(OUTLIB)MonoJet.o: $(SRCDIR)MonoJet.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MonoJet.o $<

$(OUTLIB)SF_Filler.o: $(SRCDIR)SF_Filler.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SF_Filler.o $<

$(OUTLIB)H4b.o: $(SRCDIR)H4b.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)H4b.o $<

$(OUTLIB)RSAnalysis.o: $(SRCDIR)RSAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RSAnalysis.o $<

$(OUTLIB)RSZZAnalysis.o: $(SRCDIR)RSZZAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RSZZAnalysis.o $<

$(OUTLIB)TestAnalysis.o: $(SRCDIR)TestAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TestAnalysis.o $<

$(OUTLIB)ThiagoAnalysis.o: $(SRCDIR)ThiagoAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ThiagoAnalysis.o $<

$(OUTLIB)VecbosAnalysis.o: $(SRCDIR)VecbosAnalysis.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosAnalysis.o $<

$(OUTLIB)RedVecbosTree.o: $(SRCDIR)RedVecbosTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedVecbosTree.o $<

$(OUTLIB)RedTopTree.o: $(SRCDIR)RedTopTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedTopTree.o $<

$(OUTLIB)RedVecbosPFTree.o: $(SRCDIR)RedVecbosPFTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedVecbosPFTree.o $<

$(OUTLIB)RedEleIdIsolPFTree.o: $(SRCDIR)RedEleIdIsolPFTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIdIsolPFTree.o $<

$(OUTLIB)RedVecbosMcTree.o: $(SRCDIR)RedVecbosMcTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedVecbosMcTree.o $<

$(OUTLIB)RedVecbosVertexTree.o: $(SRCDIR)RedVecbosVertexTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedVecbosVertexTree.o $<

$(OUTLIB)RedIsolationVtxOptimTree.o: $(SRCDIR)RedIsolationVtxOptimTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedIsolationVtxOptimTree.o $<

$(OUTLIB)RedEleIDOptimTree.o: $(SRCDIR)RedEleIDOptimTree.cc $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDOptimTree.o $<

$(OUTLIB)VecbosEESelection.o: $(SRCDIR)VecbosEESelection.cc $(OUTLIB)Vecbos.o \
	$(OUTLIB)JetCounter.o $(OUTLIB)McTruthEvent.o $(OUTLIB)CutBasedSelectorEE.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosEESelection.o $<

$(OUTLIB)TopControlSample.o: $(SRCDIR)TopControlSample.cc $(OUTLIB)Vecbos.o \
	$(OUTLIB)JetCounter.o $(OUTLIB)McTruthEvent.o $(OUTLIB)CutBasedSelectorEE.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TopControlSample.o $<

$(OUTLIB)LFJetControlSample.o: $(SRCDIR)LFJetControlSample.cc $(OUTLIB)Vecbos.o \
	$(OUTLIB)JetCounter.o $(OUTLIB)McTruthEvent.o $(OUTLIB)CutBasedSelectorEE.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LFJetControlSample.o $<

$(OUTLIB)VecbosPFEESelection.o: $(SRCDIR)VecbosPFEESelection.cc $(OUTLIB)Vecbos.o \
	$(OUTLIB)JetCounter.o $(OUTLIB)McTruthEvent.o $(OUTLIB)CutBasedSelectorEE.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosPFEESelection.o $<

$(OUTLIB)VecbosTTCounter.o: $(SRCDIR)VecbosTTCounter.cc $(OUTLIB)VecbosBase.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosTTCounter.o $<

$(OUTLIB)CreateWJetDataset.o: $(SRCDIR)CreateWJetDataset.cc $(OUTLIB)VecbosBase.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CreateWJetDataset.o $<

$(OUTLIB)LQ3Helper.o: $(SRCDIR)LQ3Helper.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LQ3Helper.o $<

$(OUTLIB)LQ3Analysis.o: $(SRCDIR)LQ3Analysis.cc $(OUTLIB)Vecbos.o $(OUTLIB)LQ3Helper.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LQ3Analysis.o $<

$(OUTLIB)SUSYNLO.o: $(SRCDIR)SUSYNLO.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SUSYNLO.o $<

$(OUTLIB)VBTFLeptEff.o: $(SRCDIR)VBTFLeptEff.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VBTFLeptEff.o $<

#$(OUTLIB)VecbosMuMuSelection.o: $(SRCDIR)VecbosMuMuSelection.cc $(OUTLIB)Vecbos.o \
#	$(OUTLIB)JetCounter.o $(OUTLIB)McTruthEvent.o $(OUTLIB)CutBasedSelectorEE.o
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosMuMuSelection.o $<

# main VecBos functions
$(OUTLIB)Vecbos.o: $(SRCDIR)Vecbos.cc $(OUTLIB)VecbosBase.o $(OUTLIB)Jet.o  \
	$(OUTLIB)MET.o $(OUTLIB)CaloTower.o $(OUTLIB)CoolTools.o $(OUTLIB)JetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Vecbos.o $<

# base class with ntuple structure
$(OUTLIB)VecbosBase.o: $(SRCDIR)VecbosBase.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosBase.o $<

# auxiliary functions (MET, JEt, CaloTower) for offline reclustering
$(OUTLIB)CoolTools.o: $(SRCDIR)CoolTools.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CoolTools.o $<

$(OUTLIB)Jet.o: $(SRCDIR)Jet.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Jet.o $<

$(OUTLIB)JetCounter.o: $(SRCDIR)JetCounter.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCounter.o $<

$(OUTLIB)JetCorrectorParameters.o: $(SRCDIR)JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectorParameters.o $<

$(OUTLIB)SimpleJetCorrectionUncertainty.o: $(SRCDIR)SimpleJetCorrectionUncertainty.cc \
	$(OUTLIB)JetCorrectorParameters.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrectionUncertainty.o $<

$(OUTLIB)JetCorrectionUncertainty.o: $(SRCDIR)JetCorrectionUncertainty.cc \
	$(OUTLIB)SimpleJetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectionUncertainty.o $<

$(OUTLIB)McTruthEvent.o: $(SRCDIR)McTruthEvent.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)McTruthEvent.o $<

$(OUTLIB)CutBasedSelectorEE.o: $(SRCDIR)CutBasedSelectorEE.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedSelectorEE.o $<

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

# GEN VECBOS STUSY
$(OUTLIB)GenVecbos.o: $(SRCDIR)GenVecbos.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenVecbos.o $<

$(OUTLIB)GenWjets.o: $(SRCDIR)GenWjets.C $(OUTLIB)GenVecbos.o $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenWjets.o $<

$(OUTLIB)GenZjets.o: $(SRCDIR)GenZjets.C $(OUTLIB)GenVecbos.o $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenZjets.o $<

VecbosApp.clean:
	rm -f VecbosApp

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o $(OUTLIBEGAMMA)*.o
	rm -f VecbosApp
