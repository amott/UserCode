//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef CandleCalibee_h
#define CandleCalibee_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"

using namespace std;

class CandleCalibee : public Vecbos{
public:

  CandleCalibee(TTree *tree=0); /// Class Constructor
  virtual ~CandleCalibee();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void UsePF(bool useit) {usePF = useit;}

private:
  vector<TLorentzVector> Zcand;
  vector<int> iZDaugh1;
  vector<int> iZDaugh2;

  double CalcDxyPV(int iMu, int iPV);
  double CalcErrDxyPV(int iMu, int iPV);
  double CalcDzPV(int iMu, int iPV);
  double CalcErrDzPV(int iMu, int iPV);
  void   EraseZ(int iZ);

  int  BestPV(int bestzIdx);
  double SumPt(int iMu, int iZ);
  double pTEle(int imu);
  vector<TH1D*> CreateHistos1D(string dirname);
  void FillHistos(vector<TH1D*> histos, int iZ, int iPV);
  double DeltaPhi_PiHalf(double phi1, double phi2);
  bool IsItTheElectron(Jet jet, int iEle);

  double weight;
  bool usePF;

  bool b_zz;
  bool b_ww;
  bool b_wz;
  int pTbin;
  int nParton;
  double w_weight[6]; 
  double z_weight[6]; 
  double ww_weight[6]; 
  double wz_weight[6]; 
  double zz_weight[6]; 

};
#endif
