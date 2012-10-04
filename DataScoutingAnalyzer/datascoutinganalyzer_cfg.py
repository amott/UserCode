import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'

process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFNoPUJetsProducerSequence_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")
process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "offlinePrimaryVertices"

# --- noise filters ---
process.load("HiggsAnalysis.HiggsToWW2e.METOptionalFilterFlags_cff")
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.electronAndPhotonSuperClusters_cff")
process.load("HiggsAnalysis.HiggsToWW2e.seedBasicClusters_cff")

# --- good vertex filter ---
process.load("HiggsAnalysis.HiggsToWW2e.vertexFiltering_cff")

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_0.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_1.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_10.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_11.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_12.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_13.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_14.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_15.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_16.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_17.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_18.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_19.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_2.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_20.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_21.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_22.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_23.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_24.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_25.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_26.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_27.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_28.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_29.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_3.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_30.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_31.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_32.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_33.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_34.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_35.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_36.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_37.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_38.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_39.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_4.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_40.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_41.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_42.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_43.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_44.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_45.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_46.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_47.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_5.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_6.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_7.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_8.root',
'/store/cmst3/user/wreece/CMG/JetHT/Run2012B-v1/RAW/V00-01-03/AOD/outputPhysicsDST_HLTplusAOD_9.root',
)
)

ak5CaloL2Relative = cms.ESSource('LXXXCorrectionService',
                                 era = cms.string('Jec11V12'),
                                 section   = cms.string(''),
                                 level     = cms.string('L2Relative'),
                                 # the above 3 elements are needed only when the service is initialized from local txt files
                                 algorithm = cms.string('AK5'),
                                 # the 'algorithm' tag is also the name of the DB payload
                                 useCondDB = cms.untracked.bool(True)
                                 )


process.correctedJets = cms.EDProducer('CaloJetCorrectionProducer',
                                       src = cms.InputTag('ak5CaloJets'),
                                       correctors = cms.vstring( 'ak5CaloJetsL2L3')
                                       )

process.demo = cms.EDAnalyzer('PFJetScoutingAnalyzer',
                              #jets = cms.InputTag("ak5PFL1FastL2L3"),
                              jets = cms.InputTag("ak5PFJets"),
                              jetCorrections = cms.string("ak5PFL1FastL2L3Residual"),
                              rho  = cms.InputTag("kt6PFJets","rho"),
                              jetThreshold = cms.double(20),
                              met = cms.InputTag("met"),
                              electrons = cms.InputTag(""),
                              muons = cms.InputTag(""),
                              noise = cms.InputTag(""),
                              outputFile = cms.string("test.root")

)


process.p = cms.Path(process.goodPrimaryVertices * process.metOptionalFilterSequence * process.pfIsolationAllSequence * process.ourJetSequenceDataReduced * process.demo)
