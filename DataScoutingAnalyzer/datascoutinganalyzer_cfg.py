import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_39.root',
    )
)

process.demo = cms.EDAnalyzer('DataScoutingAnalyzer',
                              jets = cms.InputTag("ak5CaloJets"),
                              rho  = cms.InputTag("kt6CaloJets","rhos"),
                              jetMatchingThreshold = cms.double(20),
                              met = cms.InputTag("met"),
                              electrons = cms.InputTag(""),
                              muons = cms.InputTag(""),
                              noise = cms.InputTag(""),
                              outputFile = cms.string("test.root")

)


process.p = cms.Path(process.demo)
