#ifndef VecbosTTCounter_h
#define VecbosTTCounter_h

#include "VecbosBase.hh"

class VecbosTTCounter : public VecbosBase{

public:

  VecbosTTCounter(TTree *tree=0);
  virtual ~VecbosTTCounter();

  //! loop over events
  void Loop();

private:

  float n_WbWb;
  float n_ee, n_emu, n_mumu;
  float n_taue, n_taumu, n_tautau;
  float n_tot;

};

#endif

