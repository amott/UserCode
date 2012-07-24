#!/bin/sh

prefix=$1

echo "submitting signals..." 
python cmst3_submit_manyfilesperjob_eeAnalysis.py WPYTHIA WToENu_TuneZ2 10 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py WPYTHIA WToMuNu_TuneZ2 10 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py WPYTHIA WToTauNu_TuneZ2_7TeV-pythia6-tauola 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py ZPYTHIA DYToEE_M-20_CT10_TuneZ2_PU 10 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py ZPYTHIA DYToMuMu_M-20_CT10_TuneZ2_PU 10 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py ZPYTHIA DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola 10 VecbosApp 4 8nh $prefix
echo "done signals" 

echo "submitting QCD di-jets..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-20to30_EMEnriched_TuneZ2 30 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-30to80_EMEnriched_TuneZ2 64 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-80to170_EMEnriched_TuneZ2 17 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-20to30_BCtoE_TuneZ2 8 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-30to80_BCtoE_TuneZ2 8 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py QCD QCD_Pt-80to170_BCtoE_TuneZ2 8 VecbosApp 4 8nh $prefix
echo "done QCD di-jets."

echo "SUBMITTING PHOTON + JETS..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G1Jet_Pt-120to180_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G1Jet_Pt-180to240_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G1Jet_Pt-240to300_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G1Jet_Pt-300to5000_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G1Jet_Pt-60to120_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-120to180_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-180to240_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-20to60_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-240to300_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-300to5000_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G2Jets_Pt-60to120_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-120to180_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-180to240_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-20to60_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-240to300_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-300to5000_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G3Jets_Pt-60to120_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-120to180_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-180to240_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-20to60_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-240to300_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-300to5000_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py PhotonJet G4Jets_Pt-60to120_TuneZ2_7TeV 5 VecbosApp 4 8nh $prefix
echo "DONE PHOTON + JETS."

echo "SUBMITTING TTBAR..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py TTbar TTJets_TuneD6T 10 VecbosApp 4 8nh $prefix
echo "DONE TTBAR."

echo "SUBMITTING SINGLE TOP..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py SingleTop TToBLNu_TuneZ2_s-channel_7TeV-madgraph 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py SingleTop TToBLNu_TuneZ2_t-channel_7TeV-madgraph 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py SingleTop TToBLNu_TuneZ2_tW-channel_7TeV-madgraph 5 VecbosApp 4 8nh $prefix
echo "DONE SINGLE TOP."

echo "submitting dibosons..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py DiBosons WWTo2L2Nu_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py DiBosons WZTo3LNu_TuneZ2 5 VecbosApp 4 8nh $prefix
python cmst3_submit_manyfilesperjob_eeAnalysis.py DiBosons ZZtoAnything_TuneZ2 5 VecbosApp 4 8nh $prefix
echo "done with dibosons."

