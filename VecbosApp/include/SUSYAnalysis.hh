//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

/// The RoganAnalysis class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef SUSYAnalysis_h
#define SUSYAnalysis_h

class SUSYAnalysis : public Vecbos{
public:

  SUSYAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~SUSYAnalysis();     /// Class Destructor
  /// The function to run on each events
  void SUSYAnalysis::Loop();
    
private:
  
  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TProfile*> CreateHistosT(string dirname); 
  vector<TH1D*>     CreateHistosD(string dirname); 
  vector<TH2D*>     CreateHistos2D(string dirname);

  /// Fill the histograms passes as input
  void FillHistos(vector<TH1D*> h_d, vector<TProfile*> h_p, vector<TH2D*> h_2D);

  double x0; ///< the X coordinate of the new vertex
  double y0; ///< the Y coordinate of the new vertex
  double z0; ///< the Z coordinate of the new vertex

  double ECHF;
  double EEMF;

  bool ISO;

  vector<CaloTower> c_uncorr;
  vector<CaloTower> c_uncorr_v;
  vector<CaloTower> c_fixed;
  vector<CaloTower> c_var;
  
  
  vector<Jet> j_uncorr;
  vector<Jet> j_corr_0;
  vector<Jet> j_corr_1;
  vector<Jet> j_corr_2;
  vector<Jet> j_corr_3;
  
  vector<Jet> j_uncorr_all;
  vector<Jet> j_corr_0_all;
  vector<Jet> j_corr_1_all;
  vector<Jet> j_corr_2_all;
  vector<Jet> j_corr_3_all;

  vector<MET> m_uncorr;
  vector<MET> m_uncorr_v;
  vector<MET> m_fixed;
  vector<MET> m_var;
  vector<MET> m_gen;

  vector<Jet> j_dum;
  vector<Jet> j_dum_all;
  vector<MET> m_dum;
  
};
#endif
