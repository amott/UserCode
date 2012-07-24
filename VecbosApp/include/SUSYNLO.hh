#ifndef SUSYNLO_h
#define SUSYNLO_h

#include <TH2D.h>

class SUSYNLO {
public:
  /// constructor
  SUSYNLO(TString filename, TString label, int iM0, double minM0, double maxM0, int iM12, double minM12, double maxM12);
  /// destructor
  ~SUSYNLO();
  /// Initialize the xsec values from the file
  void SetXsec();
  /// Get Xsec k factors from file
  void SetKfactor();
  /// Get Xsec/kFactor values for a given process and m0-m12 point
  double GetVal(double m0, double m12, int iP1, int iP2);

private:

  bool NeuCha(int iP);
  bool Squark(int iP);
  bool Slepton(int iP);
  bool Gluino(int iP);
  bool Stop(int iP);
  bool Sbottom(int iP);

  TString _filename;

  TH2D* ng;
  TH2D* ns;
  TH2D* nn;
  TH2D* ll;
  TH2D* sb;
  TH2D* ss;
  TH2D* tb;
  TH2D* bb;
  TH2D* gg;
  TH2D* sg;

};

#endif
