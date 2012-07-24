//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef CandleCalib_h
#define CandleCalib_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class CandleCalib : public Vecbos{
public:

  CandleCalib(TTree *tree=0); /// Class Constructor
  virtual ~CandleCalib();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }

private:
  vector<TLorentzVector> Zcand;
  vector<int> iZDaugh1;
  vector<int> iZDaugh2;
  
  void   EraseZ(int iZ);

  double SumPt(int iMu, int iZ);
  int  BestPV(int bestzIdx);
  double pTMuon(int imu);
  vector<TH1D*> CreateHistos1D(string dirname);
  void FillHistos(vector<TH1D*> histos, int iZ);
  double DeltaPhi_PiHalf(double phi1, double phi2);

  double weight;
  vector<int> m_requiredTriggers;

};
#endif
