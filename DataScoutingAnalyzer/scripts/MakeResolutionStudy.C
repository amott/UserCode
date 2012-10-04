{

  const int nXbins = 6; //ETA
  float xBins[nXbins+1];
  const int yStartBin=4;
  const int yStopBin=8;
  int nYbins = yStopBin-yStartBin; //PT
  float yBins[nYbins+1];
  const int nZbins = 80; // RESOLUTION
  float zBins[nZbins+1];

  for(int i=0;i<nXbins+1;i++) xBins[i] = 0+i/2.;
  for(int i=0;i<nZbins+1;i++) zBins[i] = -0.5+i/80.;
  for(int i=yStartBin;i<yStopBin+1;i++) yBins[i-yStartBin] = 20.+TMath::Power(2,i)*10.;



  TH3F *JetRes = new TH3F("JetRes","",nXbins,xBins,nYbins,yBins,nZbins,zBins);
  DSComp->Project("JetRes","(recoJetPt[dsJetMatchIndex]-dsJetPt)/recoJetPt[dsJetMatchIndex]:dsJetPt:abs(dsJetEta)","dsJetMatchIndex>-1");


  setstyle();
  gStyle->SetOptTitle(1);
  TCanvas cv;
  cv.Divide(nYbins,nXbins);
  for(int iX=0;iX<nXbins;iX++){
    std::cout << xBins[iX] << " < |eta| < " << xBins[iX+1] <<std::endl;
    for(int iY=0; iY<nYbins; iY++){
      std::cout << "\t" << yBins[iY] << " < pT < " << yBins[iY+1] << std::endl;
      cv.cd(iY+iX*nYbins+1);
      TH1D* tmp = JetRes->ProjectionZ(Form("pz_%d_%d",iX,iY),iX,iX+1,iY,iY+1);
      tmp.SetTitle(Form("%0.2f < |#eta| < %0.2f  %0.1f < p_{T} < %0.1f",xBins[iX],xBins[iX+1],yBins[iY],yBins[iY+1]));
      tmp.SetXTitle("(RECO-HLT)/HLT");
      tmp.Fit("gaus");
      tmp.Draw();
		   
    }
  }

}
