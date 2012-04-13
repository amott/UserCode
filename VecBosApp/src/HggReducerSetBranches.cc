void HggReducer::setOutputBranches(){
  /////Photon
  outTree->Branch("nPhoton",&nPhoton, "nPhoton/I");

  //for MC
  outTree->Branch("photonenergySmeared",photonenergySmeared,"photonenergySmeared[nPhoton]/I");
  //for Data
  outTree->Branch("photonenergyScaled",photonenergyScaled,"photonenergyScaled[nPhoton]/I");
  

  //the default one ( from cmssw energy correction)
  outTree->Branch("photonenergydefault",photonenergydefault,"photonenergydefault[nPhoton]/F");
  

  outTree->Branch("photonpt",photonpt,"photonpt[nPhoton]/F");
  outTree->Branch("photonenergy",photonenergy,"photonenergy[nPhoton]/F");
  outTree->Branch("photonenergy0",photonenergy0,"photonenergy0[nPhoton]/F");
  outTree->Branch("photoneta",photoneta,"photoneta[nPhoton]/F");
  outTree->Branch("photonphi",photonphi,"photonphi[nPhoton]/F");
  outTree->Branch("photonvertexx",photonvertexx,"photonvertexx[nPhoton]/F");
  outTree->Branch("photonvertexy",photonvertexy,"photonvertexy[nPhoton]/F");
  outTree->Branch("photonvertexz",photonvertexz,"photonvertexz[nPhoton]/F");
  outTree->Branch("photonhasPixelSeed",photonhasPixelSeed,"photonhasPixelSeed[nPhoton]/I");
  outTree->Branch("photonhasConversionTracks",photonhasConversionTracks,"photonhasConversionTracks[nPhoton]/I");
  outTree->Branch("photonscrawEnergy",photonscrawEnergy,"photonscrawEnergy[nPhoton]/F");
  outTree->Branch("photonscpreshowerEnergy",photonscpreshowerEnergy,"photonscpreshowerEnergy[nPhoton]/F");
  outTree->Branch("photonscphiWidth",photonscphiWidth,"photonscphiWidth[nPhoton]/F");
  outTree->Branch("photonscetaWidth",photonscetaWidth,"photonscetaWidth[nPhoton]/F");
  outTree->Branch("photonscenergy",photonscenergy,"photonscenergy[nPhoton]/F");
  
  
  outTree->Branch("photonsceta",photonsceta,"photonsceta[nPhoton]/F");
  outTree->Branch("photoncaloeta",photoncaloeta,"photoncaloeta[nPhoton]/F");
  outTree->Branch("photonscphi",photonscphi,"photonscphi[nPhoton]/F");
  outTree->Branch("photone3x3",photone3x3,"photone3x3[nPhoton]/F");
  outTree->Branch("photone1x5",photone1x5,"photone1x5[nPhoton]/F");
  outTree->Branch("photone2x5",photone2x5,"photone2x5[nPhoton]/F");
  outTree->Branch("photone5x5",photone5x5,"photone5x5[nPhoton]/F");
  outTree->Branch("photonmaxEnergyXtal",photonmaxEnergyXtal,"photonmaxEnergyXtal[nPhoton]/F");
  outTree->Branch("photonr9",photonr9,"photonr9[nPhoton]/F");
  outTree->Branch("photonr90",photonr90,"photonr90[nPhoton]/F"); ///the un-scaled one
  outTree->Branch("photonr9Scaled",photonr9Scaled,"photonr9Scaled[nPhoton]/I");
  

  outTree->Branch("photonsigmaIetaIeta",photonsigmaIetaIeta,"photonsigmaIetaIeta[nPhoton]/F");
  outTree->Branch("photonhadronicOverEm",photonhadronicOverEm,"photonhadronicOverEm[nPhoton]/F");
  outTree->Branch("photonecalRecHitSumEtConeDR03",photonecalRecHitSumEtConeDR03,"photonecalRecHitSumEtConeDR03[nPhoton]/F");
  outTree->Branch("photonhcalDepth1TowerSumEtConeDR03",photonhcalDepth1TowerSumEtConeDR03,"photonhcalDepth1TowerSumEtConeDR03[nPhoton]/F");
  outTree->Branch("photonhcalDepth2TowerSumEtConeDR03",photonhcalDepth2TowerSumEtConeDR03,"photonhcalDepth2TowerSumEtConeDR03[nPhoton]/F");
  outTree->Branch("photonhcalTowerSumEtConeDR03",photonhcalTowerSumEtConeDR03,"photonhcalTowerSumEtConeDR03[nPhoton]/F");
  outTree->Branch("photontrkSumPtHollowConeDR03",photontrkSumPtHollowConeDR03,"photontrkSumPtHollowConeDR03[nPhoton]/F");
  outTree->Branch("photontrkSumPtSolidConeDR03",photontrkSumPtSolidConeDR03,"photontrkSumPtSolidConeDR03[nPhoton]/F");
  outTree->Branch("photonnTrkHollowConeDR03",photonnTrkHollowConeDR03,"photonnTrkHollowConeDR03[nPhoton]/I");
  outTree->Branch("photonnTrkSolidConeDR03",photonnTrkSolidConeDR03,"photonnTrkSolidConeDR03[nPhoton]/I");
  outTree->Branch("photonecalRecHitSumEtConeDR04",photonecalRecHitSumEtConeDR04,"photonecalRecHitSumEtConeDR04[nPhoton]/F");
  outTree->Branch("photonhcalDepth1TowerSumEtConeDR04",photonhcalDepth1TowerSumEtConeDR04,"photonhcalDepth1TowerSumEtConeDR04[nPhoton]/F");
  outTree->Branch("photonhcalDepth2TowerSumEtConeDR04",photonhcalDepth2TowerSumEtConeDR04,"photonhcalDepth2TowerSumEtConeDR04[nPhoton]/F");
  outTree->Branch("photonhcalTowerSumEtConeDR04",photonhcalTowerSumEtConeDR04,"photonhcalTowerSumEtConeDR04[nPhoton]/F");
  outTree->Branch("photontrkSumPtHollowConeDR04",photontrkSumPtHollowConeDR04,"photontrkSumPtHollowConeDR04[nPhoton]/F");
 outTree->Branch("photontrkSumPtSolidConeDR04",photontrkSumPtSolidConeDR04,"photontrkSumPtSolidConeDR04[nPhoton]/F");
 outTree->Branch("photonnTrkHollowConeDR04",photonnTrkHollowConeDR04,"photonnTrkHollowConeDR04[nPhoton]/I");
 outTree->Branch("photonnTrkSolidConeDR04",photonnTrkSolidConeDR04,"photonnTrkSolidConeDR04[nPhoton]/I");
 outTree->Branch("photonseedtime",photonseedtime,"photonseedtime[nPhoton]/F");
 outTree->Branch("photonseedoutOfTimeChi2",photonseedoutOfTimeChi2,"photonseedoutOfTimeChi2[nPhoton]/F");
 outTree->Branch("photonseedchi2",photonseedchi2,"photonseedchi2[nPhoton]/F");
 outTree->Branch("photonseedrecoFlag",photonseedrecoFlag,"photonseedrecoFlag[nPhoton]/I");
 outTree->Branch("photonseedseverityLevel",photonseedseverityLevel,"photonseedseverityLevel[nPhoton]/I");
 outTree->Branch("photonfiducialFlag",photonfiducialFlag,"photonfiducialFlag[nPhoton]/I");
 outTree->Branch("photonconversionsize",photonconversionsize,"photonconversionsize[nPhoton]/I");
 outTree->Branch("photonscindex",photonscindex,"photonscindex[nPhoton]/I");
 
 outTree->Branch("photonhasMatchedPromptElectron",photonhasMatchedPromptElectron,"photonhasMatchedPromptElectron[nPhoton]/I");
 

 if(!_isData){ //only add these on MC
   outTree->Branch("photongenphtmatch",photongenphtmatch,"photongenphtmatch[nPhoton][3]/F");
   outTree->Branch("photongenelematch",photongenelematch,"photongenelematch[nPhoton][3]/F");
   outTree->Branch("photongenphtconv",photongenphtconv,"photongenphtconv[nPhoton][4]/F");
   outTree->Branch("photonPartonmatch",photonPartonmatch,"photonPartonmatch[nPhoton][6]/F");
   outTree->Branch("photonPartonmatchp",photonPartonmatchp,"photonPartonmatchp[nPhoton][3]/F");
 }

 
 outTree->Branch("photonInEB",photonInEB,"photonInEB[nPhoton]/I");
 outTree->Branch("photonieta",photonieta,"photonieta[nPhoton]/I");
 outTree->Branch("photoniphi",photoniphi,"photoniphi[nPhoton]/I");
 outTree->Branch("photonswissCross",photonswissCross,"photonswissCross[nPhoton]/F");
 outTree->Branch("photonE2overE9",photonE2overE9,"photonE2overE9[nPhoton]/F");
 
 outTree->Branch("photoncaloPositionx",photoncaloPositionx,"photoncaloPositionx[nPhoton]/F");
 outTree->Branch("photoncaloPositiony",photoncaloPositiony,"photoncaloPositiony[nPhoton]/F");
 outTree->Branch("photoncaloPositionz",photoncaloPositionz,"photoncaloPositionz[nPhoton]/F");
 
 outTree->Branch("photonsccaloPositionx",photonsccaloPositionx,"photonsccaloPositionx[nPhoton]/F");
 outTree->Branch("photonsccaloPositiony",photonsccaloPositiony,"photonsccaloPositiony[nPhoton]/F");
 outTree->Branch("photonsccaloPositionz",photonsccaloPositionz,"photonsccaloPositionz[nPhoton]/F");
 

 ///Photon's conversion
 outTree->Branch("photonconversionsize",photonconversionsize,"photonconversionsize[nPhoton]/I");
 outTree->Branch("photonconversionVertexx",photonconversionVertexx,"photonconversionVertexx[nPhoton]/F");
 outTree->Branch("photonconversionVertexy",photonconversionVertexy,"photonconversionVertexy[nPhoton]/F");
 outTree->Branch("photonconversionVertexz",photonconversionVertexz,"photonconversionVertexz[nPhoton]/F");
 
 outTree->Branch("photonconversionrefittedPairMomentumx",photonconversionrefittedPairMomentumx,"photonconversionrefittedPairMomentumx[nPhoton]/F");
 outTree->Branch("photonconversionrefittedPairMomentumy",photonconversionrefittedPairMomentumy,"photonconversionrefittedPairMomentumy[nPhoton]/F");
 outTree->Branch("photonconversionrefittedPairMomentumz",photonconversionrefittedPairMomentumz,"photonconversionrefittedPairMomentumz[nPhoton]/F");
 
 outTree->Branch("photonconversionpairInvariantMass",photonconversionpairInvariantMass,"photonconversionpairInvariantMass[nPhoton]/F");
 outTree->Branch("photonconversionpairCotThetaSeparation",photonconversionpairCotThetaSeparation,"photonconversionpairCotThetaSeparation[nPhoton]/F");
 outTree->Branch("photonconversionEoverPrefittedTracks",photonconversionEoverPrefittedTracks,"photonconversionEoverPrefittedTracks[nPhoton]/F");
 outTree->Branch("photonconversionzOfPrimaryVertexFromTracks",photonconversionzOfPrimaryVertexFromTracks,"photonconversionzOfPrimaryVertexFromTracks[nPhoton]/F");
 outTree->Branch("photonconversiondistOfMinimumApproach",photonconversiondistOfMinimumApproach,"photonconversiondistOfMinimumApproach[nPhoton]/F");
 outTree->Branch("photonconversiondPhiTracksAtVtx",photonconversiondPhiTracksAtVtx,"photonconversiondPhiTracksAtVtx[nPhoton]/F");
 outTree->Branch("photonconversiondPhiTracksAtEcal",photonconversiondPhiTracksAtEcal,"photonconversiondPhiTracksAtEcal[nPhoton]/F");
 outTree->Branch("photonconversiondEtaTracksAtEcal",photonconversiondEtaTracksAtEcal,"photonconversiondEtaTracksAtEcal[nPhoton]/F");
 outTree->Branch("photonconversionnTracks",photonconversionnTracks,"photonconversionnTracks[nPhoton]/I");
 outTree->Branch("photonconversionMVAout",photonconversionMVAout,"photonconversionMVAout[nPhoton]/F");
 outTree->Branch("photonconversionVertexisValid",photonconversionVertexisValid,"photonconversionVertexisValid[nPhoton]/I");
 outTree->Branch("photonconversionVertexchi2",photonconversionVertexchi2,"photonconversionVertexchi2[nPhoton]/F");
 outTree->Branch("photonconversionChiSquaredProbability",photonconversionChiSquaredProbability,"photonconversionChiSquaredProbability[nPhoton]/F");
 outTree->Branch("photonconversion_track1_dz",photonconversion_track1_dz,"photonconversion_track1_dz[nPhoton]/F");
 outTree->Branch("photonconversion_track1_dzError",photonconversion_track1_dzError,"photonconversion_track1_dzError[nPhoton]/F");
 outTree->Branch("photonconversion_track1_charge",photonconversion_track1_charge,"photonconversion_track1_charge[nPhoton]/I");
 outTree->Branch("photonconversion_track1_d0",photonconversion_track1_d0,"photonconversion_track1_d0[nPhoton]/F");
 outTree->Branch("photonconversion_track1_tracksPout",photonconversion_track1_tracksPout,"photonconversion_track1_tracksPout[nPhoton]/F");
 outTree->Branch("photonconversion_track1_tracksPin",photonconversion_track1_tracksPin,"photonconversion_track1_tracksPin[nPhoton]/F");
 outTree->Branch("photonconversion_track1_algo",photonconversion_track1_algo,"photonconversion_track1_algo[nPhoton]/I");


 outTree->Branch("photonconversion_track2_dz",photonconversion_track2_dz,"photonconversion_track2_dz[nPhoton]/F");
 outTree->Branch("photonconversion_track2_dzError",photonconversion_track2_dzError,"photonconversion_track2_dzError[nPhoton]/F");
 outTree->Branch("photonconversion_track2_charge",photonconversion_track2_charge,"photonconversion_track2_charge[nPhoton]/I");
 outTree->Branch("photonconversion_track2_d0",photonconversion_track2_d0,"photonconversion_track2_d0[nPhoton]/F");
 outTree->Branch("photonconversion_track2_tracksPout",photonconversion_track2_tracksPout,"photonconversion_track2_tracksPout[nPhoton]/F");
 outTree->Branch("photonconversion_track2_tracksPin",photonconversion_track2_tracksPin,"photonconversion_track2_tracksPin[nPhoton]/F");
 outTree->Branch("photonconversion_track2_algo",photonconversion_track2_algo,"photonconversion_track2_algo[nPhoton]/I");
 

 //more added
 outTree->Branch("photonscnhits",photonscnhits,"photonscnhits[nPhoton]/I");
 outTree->Branch("photonscclusterSize",photonscclusterSize,"photonscclusterSize[nPhoton]/I");
 

 //explicitly give the class name for these
 outTree->Branch("photonscbclusterenergy", "std::vector<std::vector<float> >", &photonscbclusterenergy);
 outTree->Branch("photonscbclusterposx", "std::vector<std::vector<float> >", &photonscbclusterposx);
 outTree->Branch("photonscbclusterposy", "std::vector<std::vector<float> >", &photonscbclusterposy);
 outTree->Branch("photonscbclusterposz", "std::vector<std::vector<float> >", &photonscbclusterposz);
 
 
 outTree->Branch("photonscbcseedind",photonscbcseedind,"photonscbcseedind[nPhoton]/I");
 outTree->Branch("photonscbc2ind",photonscbc2ind,"photonscbc2ind[nPhoton]/I");
 outTree->Branch("photonscbclastind",photonscbclastind,"photonscbclastind[nPhoton]/I");
 outTree->Branch("photonscbclast2ind",photonscbclast2ind,"photonscbclast2ind[nPhoton]/I");
 outTree->Branch("photonbcseedenergy",photonbcseedenergy,"photonbcseedenergy[nPhoton]/F");
 outTree->Branch("photonbcseedeta",photonbcseedeta,"photonbcseedeta[nPhoton]/F");
 outTree->Branch("photonbcseedphi",photonbcseedphi,"photonbcseedphi[nPhoton]/F");
 outTree->Branch("photonbcseedeMax",photonbcseedeMax,"photonbcseedeMax[nPhoton]/F");
 outTree->Branch("photonbcseede2nd",photonbcseede2nd,"photonbcseede2nd[nPhoton]/F");
 outTree->Branch("photonbcseedeLeft",photonbcseedeLeft,"photonbcseedeLeft[nPhoton]/F");
 outTree->Branch("photonbcseedeRight",photonbcseedeRight,"photonbcseedeRight[nPhoton]/F");
 outTree->Branch("photonbcseedeTop",photonbcseedeTop,"photonbcseedeTop[nPhoton]/F");
 outTree->Branch("photonbcseedeBottom",photonbcseedeBottom,"photonbcseedeBottom[nPhoton]/F");
 outTree->Branch("photonbcseede3x3",photonbcseede3x3,"photonbcseede3x3[nPhoton]/F");
 outTree->Branch("photonbcseede5x5",photonbcseede5x5,"photonbcseede5x5[nPhoton]/F");
 outTree->Branch("photonbcseedCovIEtaIEta",photonbcseedCovIEtaIEta,"photonbcseedCovIEtaIEta[nPhoton]/F");
 outTree->Branch("photonbcseedCovIEtaIPhi",photonbcseedCovIEtaIPhi,"photonbcseedCovIEtaIPhi[nPhoton]/F");
 outTree->Branch("photonbcseedCovIPhiIPhi",photonbcseedCovIPhiIPhi,"photonbcseedCovIPhiIPhi[nPhoton]/F");
 outTree->Branch("photonbcseedetacry",photonbcseedetacry,"photonbcseedetacry[nPhoton]/F");
 outTree->Branch("photonbcseedphicry",photonbcseedphicry,"photonbcseedphicry[nPhoton]/F");
 outTree->Branch("photonbcseedieta",photonbcseedieta,"photonbcseedieta[nPhoton]/I");
 outTree->Branch("photonbcseediphi",photonbcseediphi,"photonbcseediphi[nPhoton]/I");
 

 outTree->Branch("photonbc2energy",photonbc2energy,"photonbc2energy[nPhoton]/F");
 outTree->Branch("photonbc2eta",photonbc2eta,"photonbc2eta[nPhoton]/F");
 outTree->Branch("photonbc2phi",photonbc2phi,"photonbc2phi[nPhoton]/F");
 outTree->Branch("photonbc2eMax",photonbc2eMax,"photonbc2eMax[nPhoton]/F");
 outTree->Branch("photonbc2e2nd",photonbc2e2nd,"photonbc2e2nd[nPhoton]/F");
 outTree->Branch("photonbc2eLeft",photonbc2eLeft,"photonbc2eLeft[nPhoton]/F");
 outTree->Branch("photonbc2eRight",photonbc2eRight,"photonbc2eRight[nPhoton]/F");
 outTree->Branch("photonbc2eTop",photonbc2eTop,"photonbc2eTop[nPhoton]/F");
 outTree->Branch("photonbc2eBottom",photonbc2eBottom,"photonbc2eBottom[nPhoton]/F");
 outTree->Branch("photonbc2e3x3",photonbc2e3x3,"photonbc2e3x3[nPhoton]/F");
 outTree->Branch("photonbc2e5x5",photonbc2e5x5,"photonbc2e5x5[nPhoton]/F");
 outTree->Branch("photonbc2CovIEtaIEta",photonbc2CovIEtaIEta,"photonbc2CovIEtaIEta[nPhoton]/F");
 outTree->Branch("photonbc2CovIEtaIPhi",photonbc2CovIEtaIPhi,"photonbc2CovIEtaIPhi[nPhoton]/F");
 outTree->Branch("photonbc2CovIPhiIPhi",photonbc2CovIPhiIPhi,"photonbc2CovIPhiIPhi[nPhoton]/F");
 outTree->Branch("photonbc2ieta",photonbc2ieta,"photonbc2ieta[nPhoton]/I");
 outTree->Branch("photonbc2iphi",photonbc2iphi,"photonbc2iphi[nPhoton]/I");
 outTree->Branch("photonbc2etacry",photonbc2etacry,"photonbc2etacry[nPhoton]/F");
 outTree->Branch("photonbc2phicry",photonbc2phicry,"photonbc2phicry[nPhoton]/F");
 
 
 outTree->Branch("photonbclastenergy",photonbclastenergy,"photonbclastenergy[nPhoton]/F");
 outTree->Branch("photonbclasteta",photonbclasteta,"photonbclasteta[nPhoton]/F");
 outTree->Branch("photonbclastphi",photonbclastphi,"photonbclastphi[nPhoton]/F");
 outTree->Branch("photonbclasteMax",photonbclasteMax,"photonbclasteMax[nPhoton]/F");
 outTree->Branch("photonbclaste2nd",photonbclaste2nd,"photonbclaste2nd[nPhoton]/F");
 outTree->Branch("photonbclasteLeft",photonbclasteLeft,"photonbclasteLeft[nPhoton]/F");
 outTree->Branch("photonbclasteRight",photonbclasteRight,"photonbclasteRight[nPhoton]/F");
 outTree->Branch("photonbclasteTop",photonbclasteTop,"photonbclasteTop[nPhoton]/F");
 outTree->Branch("photonbclasteBottom",photonbclasteBottom,"photonbclasteBottom[nPhoton]/F");
 outTree->Branch("photonbclaste3x3",photonbclaste3x3,"photonbclaste3x3[nPhoton]/F");
 outTree->Branch("photonbclaste5x5",photonbclaste5x5,"photonbclaste5x5[nPhoton]/F");
 outTree->Branch("photonbclastCovIEtaIEta",photonbclastCovIEtaIEta,"photonbclastCovIEtaIEta[nPhoton]/F");
 outTree->Branch("photonbclastCovIEtaIPhi",photonbclastCovIEtaIPhi,"photonbclastCovIEtaIPhi[nPhoton]/F");
 outTree->Branch("photonbclastCovIPhiIPhi",photonbclastCovIPhiIPhi,"photonbclastCovIPhiIPhi[nPhoton]/F");
 
 outTree->Branch("photonbclast2energy",photonbclast2energy,"photonbclast2energy[nPhoton]/F");
 outTree->Branch("photonbclast2eta",photonbclast2eta,"photonbclast2eta[nPhoton]/F");
 outTree->Branch("photonbclast2phi",photonbclast2phi,"photonbclast2phi[nPhoton]/F");
 outTree->Branch("photonbclast2eMax",photonbclast2eMax,"photonbclast2eMax[nPhoton]/F");
 outTree->Branch("photonbclast2e2nd",photonbclast2e2nd,"photonbclast2e2nd[nPhoton]/F");
 outTree->Branch("photonbclast2eLeft",photonbclast2eLeft,"photonbclast2eLeft[nPhoton]/F");
 outTree->Branch("photonbclast2eRight",photonbclast2eRight,"photonbclast2eRight[nPhoton]/F");
 outTree->Branch("photonbclast2eTop",photonbclast2eTop,"photonbclast2eTop[nPhoton]/F");
 outTree->Branch("photonbclast2eBottom",photonbclast2eBottom,"photonbclast2eBottom[nPhoton]/F");
 outTree->Branch("photonbclast2e3x3",photonbclast2e3x3,"photonbclast2e3x3[nPhoton]/F");
 outTree->Branch("photonbclast2e5x5",photonbclast2e5x5,"photonbclast2e5x5[nPhoton]/F");
 outTree->Branch("photonbclast2CovIEtaIEta",photonbclast2CovIEtaIEta,"photonbclast2CovIEtaIEta[nPhoton]/F");
 outTree->Branch("photonbclast2CovIEtaIPhi",photonbclast2CovIEtaIPhi,"photonbclast2CovIEtaIPhi[nPhoton]/F");
 outTree->Branch("photonbclast2CovIPhiIPhi",photonbclast2CovIPhiIPhi,"photonbclast2CovIPhiIPhi[nPhoton]/F");
 
 outTree->Branch("photone2nd",photone2nd,"photone2nd[nPhoton]/F");
 outTree->Branch("photoneLeft",photoneLeft,"photoneLeft[nPhoton]/F");
 outTree->Branch("photoneRight",photoneRight,"photoneRight[nPhoton]/F");
 outTree->Branch("photoneBottom",photoneBottom,"photoneBottom[nPhoton]/F");
 outTree->Branch("photoneTop",photoneTop,"photoneTop[nPhoton]/F");

outTree->Branch("photone1x3",photone1x3,"photone1x3[nPhoton]/F");
outTree->Branch("photone3x1",photone3x1,"photone3x1[nPhoton]/F");
outTree->Branch("photone2x2",photone2x2,"photone2x2[nPhoton]/F");
outTree->Branch("photone3x2",photone3x2,"photone3x2[nPhoton]/F");
outTree->Branch("photone4x4",photone4x4,"photone4x4[nPhoton]/F");
outTree->Branch("photone2x5Right",photone2x5Right,"photone2x5Right[nPhoton]/F");
outTree->Branch("photone2x5Left",photone2x5Left,"photone2x5Left[nPhoton]/F");
outTree->Branch("photone2x5Top",photone2x5Top,"photone2x5Top[nPhoton]/F");
outTree->Branch("photone2x5Bottom",photone2x5Bottom,"photone2x5Bottom[nPhoton]/F");
outTree->Branch("photone2x5Max",photone2x5Max,"photone2x5Max[nPhoton]/F");
outTree->Branch("photonlat",photonlat,"photonlat[nPhoton][3]/F");
outTree->Branch("photonCovEtaEta",photonCovEtaEta,"photonCovEtaEta[nPhoton]/F");
outTree->Branch("photonCovEtaPhi",photonCovEtaPhi,"photonCovEtaPhi[nPhoton]/F");
outTree->Branch("photonCovPhiPhi",photonCovPhiPhi,"photonCovPhiPhi[nPhoton]/F");
outTree->Branch("photonCovIEtaIPhi",photonCovIEtaIPhi,"photonCovIEtaIPhi[nPhoton]/F");
outTree->Branch("photonCovIPhiIPhi",photonCovIPhiIPhi,"photonCovIPhiIPhi[nPhoton]/F");
outTree->Branch("photonzernike20",photonzernike20,"photonzernike20[nPhoton]/F");
outTree->Branch("photonzernike42",photonzernike42,"photonzernike42[nPhoton]/F");

outTree->Branch("photon_convp",photon_convp,"photon_convp[nPhoton]/F");
outTree->Branch("photon_convpt",photon_convpt,"photon_convpt[nPhoton]/F");
outTree->Branch("photon_convpttrk1",photon_convpttrk1,"photon_convpttrk1[nPhoton]/F");
outTree->Branch("photon_convpttrk2",photon_convpttrk2,"photon_convpttrk2[nPhoton]/F");
outTree->Branch("photon_deltaphiconvsc",photon_deltaphiconvsc,"photon_deltaphiconvsc[nPhoton]/F");
outTree->Branch("photon_deltaetaconvsc",photon_deltaetaconvsc,"photon_deltaetaconvsc[nPhoton]/F");
outTree->Branch("photon_convrho",photon_convrho,"photon_convrho[nPhoton]/F");
outTree->Branch("photon_convz",photon_convz,"photon_convz[nPhoton]/F");


 outTree->Branch("photondrtotrk",photondrtotrk,"photondrtotrk[nPhoton]/F"); 
 outTree->Branch("photontrkisoselvtxdr03","std::vector<std::vector<float> >", &photontrkisoselvtxdr03);

 outTree->Branch("photontrkisoworstdr04",photontrkisoworstdr04,"photontrkisoworstdr04[nPhoton]/F");
 outTree->Branch("photonindvtxtrkisoworstdr04",photonindvtxtrkisoworstdr04,"photonindvtxtrkisoworstdr04[nPhoton]/I"); 
 outTree->Branch("indvertexSelected_allpairpresel","std::vector<short>",&indvertexSelected_allpairpresel);


 ///information for the vertex
 outTree->Branch("nVertex",&nVertex,"nVertex/I");
 outTree->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
 outTree->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
 outTree->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
 outTree->Branch("vertexchi2",vertexchi2,"vertexchi2[nVertex]/F");
 outTree->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
 outTree->Branch("vertexnormalizedChi2",vertexnormalizedChi2,"vertexnormalizedChi2[nVertex]/F");
 outTree->Branch("vertextrackSize",vertextrackSize,"vertextrackSize[nVertex]/I");
 outTree->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
 outTree->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");

 
 
 //physics declared -- should be set by the JSON, but can't hurt
 outTree->Branch("phyDeclared",&phyDeclared,"phyDeclared/I");
 
  
 outTree->Branch("rho", &rho,"rho/F");
 outTree->Branch("rhoEtaMax44", &rhoEtaMax44,"rhoEtaMax44/F");
 
 
 //Event info
 outTree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
 outTree->Branch("runNumber",&runNumber,"runNumber/I");
 outTree->Branch("evtNumber",&evtNumber,"evtNumber/I");
 outTree->Branch("bunchX",&bunchX,"bunchX/I");
 outTree->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
 outTree->Branch("evtTime",&evtTime,"evtTime/I");
 outTree->Branch("isRealData",&_isData,"isRealData/I");
  
   
 outTree->Branch("pileupBunchX","std::vector<short>", &pileupBunchX);
 outTree->Branch("pileupNInteraction","std::vector<short>", &pileupNInteraction);
 outTree->Branch("pileupTrueNumInterations",&pileupTrueNumInterations);
 

  //trigger -- here we depart a bit from Yong's original code
 for(int i=0;i<triggerNames->size();i++){
   outTree->Branch(triggerNames->at(i), &(triggerBits[i]), Form("%s/I",triggerNames->at(i)) );  // this will produce 1 int per trigger in the output tree
 }
 


  if(!_isData){
    //generator level information

    outTree->Branch("signalProcessID",&signalProcessID,"signalProcessID/I");
    outTree->Branch("qScale",&qScale,"qScale/F");

    ///gen electron, muon,photon
    outTree->Branch("nGenPht",&nGenPht,"nGenPht/I");
    outTree->Branch("etaGenPht",etaGenPht,"etaGenPht[nGenPht]/F");
    outTree->Branch("phiGenPht",phiGenPht,"phiGenPht[nGenPht]/F");
    outTree->Branch("ptGenPht",ptGenPht,"ptGenPht[nGenPht]/F");
    outTree->Branch("vxGenPht",vxGenPht,"vxGenPht[nGenPht]/F");
    outTree->Branch("vyGenPht",vyGenPht,"vyGenPht[nGenPht]/F");
    outTree->Branch("vzGenPht",vzGenPht,"vzGenPht[nGenPht]/F");
    outTree->Branch("pidmomGenPht",pidmomGenPht,"pidmomGenPht[nGenPht]/I");
    outTree->Branch("pidmom2GenPht",pidmom2GenPht,"pidmom2GenPht[nGenPht]/I");
    outTree->Branch("pidmom3GenPht",pidmom3GenPht,"pidmom3GenPht[nGenPht]/I");
    outTree->Branch("statusGenPht",statusGenPht,"statusGenPht[nGenPht]/I");

    outTree->Branch("nGenEle",&nGenEle,"nGenEle/I");
    outTree->Branch("chaGenEle",chaGenEle,"chaGenEle[nGenEle]/I");
    outTree->Branch("etaGenEle",etaGenEle,"etaGenEle[nGenEle]/F");
    outTree->Branch("phiGenEle",phiGenEle,"phiGenEle[nGenEle]/F");
    outTree->Branch("ptGenEle",ptGenEle,"ptGenEle[nGenEle]/F");
    outTree->Branch("pidmomGenEle",pidmomGenEle,"pidmomGenEle[nGenEle]/I");
    outTree->Branch("pidmom2GenEle",pidmom2GenEle,"pidmom2GenEle[nGenEle]/I");
    outTree->Branch("pidmom3GenEle",pidmom3GenEle,"pidmom3GenEle[nGenEle]/I");
    outTree->Branch("statusGenEle",statusGenEle,"statusGenEle[nGenEle]/I");
    outTree->Branch("vxGenEle",vxGenEle,"vxGenEle[nGenEle]/F");
    outTree->Branch("vyGenEle",vyGenEle,"vyGenEle[nGenEle]/F");
    outTree->Branch("vzGenEle",vzGenEle,"vzGenEle[nGenEle]/F");
    
    outTree->Branch("nGenMu",&nGenMu,"nGenMu/I");
    outTree->Branch("etaGenMu",etaGenMu,"etaGenMu[nGenMu]/F");
    outTree->Branch("phiGenMu",phiGenMu,"phiGenMu[nGenMu]/F");
    outTree->Branch("ptGenMu",ptGenMu,"ptGenMu[nGenMu]/F");
    outTree->Branch("pidmomGenMu",pidmomGenMu,"pidmomGenMu[nGenMu]/I");
    outTree->Branch("pidmom2GenMu",pidmom2GenMu,"pidmom2GenMu[nGenMu]/I");
    outTree->Branch("pidmom3GenMu",pidmom3GenMu,"pidmom3GenMu[nGenMu]/I");
    outTree->Branch("statusGenMu",statusGenMu,"statusGenMu[nGenMu]/I");
    outTree->Branch("vxGenMu",vxGenMu,"vxGenMu[nGenMu]/F");
    outTree->Branch("vyGenMu",vyGenMu,"vyGenMu[nGenMu]/F");
    outTree->Branch("vzGenMu",vzGenMu,"vzGenMu[nGenMu]/F");
    outTree->Branch("chaGenMu",chaGenMu,"chaGenMu[nGenMu]/I");

    //generator level higgs
    outTree->Branch("genhiggsm",&genhiggsm,"genhiggsm/F");
    outTree->Branch("genhiggspt",&genhiggspt,"genhiggspt/F");
    outTree->Branch("genhiggseta",&genhiggseta,"genhiggseta/F");
    outTree->Branch("genhiggsphi",&genhiggsphi,"genhiggsphi/F");
    outTree->Branch("genhiggsvx",&genhiggsvx,"genhiggsvx/F");
    outTree->Branch("genhiggsvy",&genhiggsvy,"genhiggsvy/F");
    outTree->Branch("genhiggsvz",&genhiggsvz,"genhiggsvz/F");
  }

}
