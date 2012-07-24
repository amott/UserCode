#include <iostream>
#include "include/VecbosTTCounter.hh"

using namespace std;

VecbosTTCounter::VecbosTTCounter(TTree *tree)
  : VecbosBase(tree) { 

  n_WbWb = 0;
  n_ee = 0;
  n_emu = 0;
  n_mumu = 0;
  n_taue = 0;
  n_taumu = 0;
  n_tautau = 0;
  n_tot = 0;

}

VecbosTTCounter::~VecbosTTCounter() { }

void VecbosTTCounter::Loop() {

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Total number of entries in the chain = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // normalization: tt -> Wb Wb
    int Wplus=-1, Wminus=-1;
    int eplus = -1, eminus = -1;
    int muplus = -1, muminus = -1;
    int tauplus = -1, tauminus = -1;
    for(int imc=0; imc<25 && imc<nMc; imc++) {
      // look for the 2 W's
      if(idMc[imc]==-24 && idMc[mothMc[imc]]==-6) Wminus = imc;
      if(idMc[imc]==24 && idMc[mothMc[imc]]==6) Wplus = imc;
      
      // look for leptons
      if(idMc[imc]==-11 && idMc[mothMc[imc]]==24 && idMc[mothMc[mothMc[imc]]]==6) eplus = imc;
      if(idMc[imc]==11 && idMc[mothMc[imc]]==-24 && idMc[mothMc[mothMc[imc]]]==-6) eminus = imc;
      
      if(idMc[imc]==-13 && idMc[mothMc[imc]]==24 && idMc[mothMc[mothMc[imc]]]==6) muplus = imc;
      if(idMc[imc]==13 && idMc[mothMc[imc]]==-24 && idMc[mothMc[mothMc[imc]]]==-6) muminus = imc;
      
      if(idMc[imc]==-15 && idMc[mothMc[imc]]==24 && idMc[mothMc[mothMc[imc]]]==6) tauplus = imc;
      if(idMc[imc]==15 && idMc[mothMc[imc]]==-24 && idMc[mothMc[mothMc[imc]]]==-6) tauminus = imc;
    }
    
    n_tot++;

    if(Wplus>-1 && Wminus>-1) {

      n_WbWb++;

      if(eplus>-1 && eminus>-1) n_ee++;
      if( (eplus>-1 && muminus>-1) || (eminus>-1 && muplus>-1) ) n_emu++;
      if(muplus>-1 && muminus>-1) n_mumu++;
      if(tauplus>-1 && tauminus>-1) n_tautau++;
      if( (tauplus>-1 && eminus>-1) || (tauminus>-1 && eplus>-1) ) n_taue++;
      if( (tauplus>-1 && muminus>-1) || (tauminus>-1 && muplus>-1 ) ) n_taumu++;

    }

  } // end loop over entries

  cout << "--- run over " << n_tot << " events" << endl;
  cout << "-----------------------------------" << endl;
  cout << "--- n_WbWb = " << n_WbWb << " (" << 100 * n_WbWb / n_tot << " %)" << endl;
  cout << "--- n_ee  = " << n_ee << " (" << 100 * n_ee / n_WbWb << " %)" << endl;
  cout << "--- n_emu = " << n_emu << " (" << 100 * n_emu / n_WbWb << " %)" << endl;
  cout << "--- n_mumu = " << n_mumu << " (" << 100 * n_mumu / n_WbWb << "%)" << endl;
  cout << "--- n_taue = " << n_taue << " (" << 100 * n_taue / n_WbWb << "%)" << endl;
  cout << "--- n_taumu = " << n_taumu << " (" << 100 * n_taumu / n_WbWb << "%)" << endl;
  cout << "--- n_tautau = " << n_tautau << " (" << 100 * n_tautau / n_WbWb << "%)" << endl;
  cout << "------------------------------------" << endl;

}
