import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

process.demo = cms.EDAnalyzer('DataScoutingAnalyzer',
                              jets = cms.InputTag("ak5CaloJets"),
                              rho  = cms.InputTag("kt6CaloJets","rho"),
                              jetThreshold = cms.double(30),
                              met = cms.InputTag("met"),
                              electrons = cms.InputTag(""),
                              muons = cms.InputTag(""),
                              noise = cms.InputTag(""),
                              outputFile = cms.string("test.root")

)


process.p = cms.Path(process.demo)
