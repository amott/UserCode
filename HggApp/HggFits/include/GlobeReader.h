//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 22 11:28:21 2013 by ROOT version 5.32/00
// from TTree event/Reduced tree
// found on file: root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/moriond2013_reduction_v1/data/DoublePhoton_Run2012B-13Jul2012-v1/DoublePhoton_Run2012B-13Jul2012-v1_0.root
//////////////////////////////////////////////////////////

#ifndef GlobeReader_h
#define GlobeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class GlobeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           lumis;
   Int_t           run;
   Int_t           bx;
   TClonesArray    *sc_p4;
   TClonesArray    *sc_xyz;
   Int_t           sc_n;
   Float_t         sc_sphi[93];   //[sc_n]
   Float_t         sc_seta[93];   //[sc_n]
   Float_t         sc_raw[93];   //[sc_n]
   Float_t         sc_pre[93];   //[sc_n]
   Int_t           sc_bcseedind[93];   //[sc_n]
   Int_t           sc_nbc[93];   //[sc_n]
   Int_t           ecalhit_n;
   Int_t           ecalhit_detid[2644];   //[ecalhit_n]
   TClonesArray    *ecalhit_p4;
   Int_t           bc_n;
   TClonesArray    *bc_p4;
   Float_t         bc_s25[352];   //[bc_n]
   Int_t           tk_n;
   Float_t         met_pfmet;
   Float_t         met_phi_pfmet;
   Float_t         met_sumet_pfmet;
   Float_t         met_mEtSig_pfmet;
   Float_t         met_significance_pfmet;
   Float_t         met_pfmetType1;
   Float_t         met_phi_pfmetType1;
   Float_t         met_sumet_pfmetType1;
   Float_t         met_mEtSig_pfmetType1;
   Float_t         met_significance_pfmetType1;
   Int_t           el_std_n;
   TClonesArray    *el_std_p4;
   TClonesArray    *el_std_sc;
   Int_t           el_std_scind[7];   //[el_std_n]
   Float_t         el_std_sieie[7];   //[el_std_n]
   Float_t         el_std_sieiesc[7];   //[el_std_n]
   Int_t           el_std_hp_expin[7];   //[el_std_n]
   Float_t         el_std_detain[7];   //[el_std_n]
   Float_t         el_std_dphiin[7];   //[el_std_n]
   Float_t         el_std_eopin[7];   //[el_std_n]
   Float_t         el_std_fbrem[7];   //[el_std_n]
   Float_t         el_std_ip_gsf[7];   //[el_std_n]
   Float_t         el_std_hoe[7];   //[el_std_n]
   Float_t         el_std_eseedopin[7];   //[el_std_n]
   Float_t         el_std_tkiso03[7];   //[el_std_n]
   Float_t         el_std_ecaliso03[7];   //[el_std_n]
   Float_t         el_std_ecaliso04[7];   //[el_std_n]
   Float_t         el_std_hcaliso03[7];   //[el_std_n]
   Float_t         el_std_hcaliso04[7];   //[el_std_n]
   Float_t         el_std_dcot[7];   //[el_std_n]
   Float_t         el_std_dist[7];   //[el_std_n]
   Int_t           el_std_tkind[7];   //[el_std_n]
   Float_t         el_std_d0[7];   //[el_std_n]
   Float_t         el_std_z0[7];   //[el_std_n]
   TClonesArray    *el_std_posvtx;
   Float_t         el_std_D0Vtx[7][100];   //[el_std_n]
   Float_t         el_std_DZVtx[7][100];   //[el_std_n]
   Float_t         el_std_ip_ctf[7];   //[el_std_n]
   Bool_t          el_std_tkdrv[7];   //[el_std_n]
   Bool_t          el_std_ecaldrv[7];   //[el_std_n]
   Float_t         el_std_pfiso_charged[7];   //[el_std_n]
   Float_t         el_std_pfiso_neutral[7];   //[el_std_n]
   Float_t         el_std_pfiso_photon[7];   //[el_std_n]
   Float_t         el_std_conv_vtxProb[7];   //[el_std_n]
   Int_t           el_std_conv[7];   //[el_std_n]
   Float_t         el_std_pin[7];   //[el_std_n]
   Float_t         el_std_mva_nontrig[7];   //[el_std_n]
   Float_t         el_std_mva_trig[7];   //[el_std_n]
   Float_t         el_std_regr_energy[7];   //[el_std_n]
   Float_t         el_std_regr_energyerr[7];   //[el_std_n]
   Int_t           mu_glo_n;
   TClonesArray    *mu_glo_p4;
   Float_t         mu_glo_dof[39];   //[mu_glo_n]
   Float_t         mu_glo_chi2[39];   //[mu_glo_n]
   Int_t           mu_glo_type[39];   //[mu_glo_n]
   Int_t           mu_glo_charge[39];   //[mu_glo_n]
   Int_t           mu_glo_nmatches[39];   //[mu_glo_n]
   TClonesArray    *mu_glo_posvtx;
   Int_t           mu_glo_losthits[39];   //[mu_glo_n]
   Int_t           mu_glo_validhits[39];   //[mu_glo_n]
   Int_t           mu_glo_innerhits[39];   //[mu_glo_n]
   Int_t           mu_glo_pixelhits[39];   //[mu_glo_n]
   Int_t           mu_glo_validChmbhits[39];   //[mu_glo_n]
   Float_t         mu_glo_tkpterr[39];   //[mu_glo_n]
   Float_t         mu_glo_ecaliso03[39];   //[mu_glo_n]
   Float_t         mu_glo_hcaliso03[39];   //[mu_glo_n]
   Float_t         mu_glo_tkiso03[39];   //[mu_glo_n]
   Float_t         mu_glo_dz[39];   //[mu_glo_n]
   Float_t         mu_glo_D0Vtx[39][100];   //[mu_glo_n]
   Float_t         mu_glo_DZVtx[39][100];   //[mu_glo_n]
   Int_t           mu_tkLayers[39];   //[mu_glo_n]
   Float_t         mu_glo_chhadiso04[39];   //[mu_glo_n]
   Float_t         mu_glo_nehadiso04[39];   //[mu_glo_n]
   Float_t         mu_glo_photiso04[39];   //[mu_glo_n]
   Float_t         mu_dbCorr[39];   //[mu_glo_n]
   Int_t           pho_n;
   Float_t         pho_feta[10][5];   //[pho_n]
   Float_t         pho_crackcorr[10];   //[pho_n]
   Float_t         pho_localcorr[10];   //[pho_n]
   Int_t           pho_isEB[10];   //[pho_n]
   Int_t           pho_isEE[10];   //[pho_n]
   Float_t         pho_see[10];   //[pho_n]
   Float_t         pho_sieie[10];   //[pho_n]
   Float_t         pho_sipip[10];   //[pho_n]
   Float_t         pho_sieip[10];   //[pho_n]
   Float_t         pho_e1x5[10];   //[pho_n]
   Float_t         pho_e3x3[10];   //[pho_n]
   Float_t         pho_e5x5[10];   //[pho_n]
   Float_t         pho_emaxxtal[10];   //[pho_n]
   Float_t         pho_hoe[10];   //[pho_n]
   Float_t         pho_r1x5[10];   //[pho_n]
   Float_t         pho_r2x5[10];   //[pho_n]
   Float_t         pho_r9[10];   //[pho_n]
   Int_t           pho_isEBGap[10];   //[pho_n]
   Int_t           pho_isEEGap[10];   //[pho_n]
   Int_t           pho_isEBEEGap[10];   //[pho_n]
   Float_t         pho_zernike20[10];   //[pho_n]
   Float_t         pho_zernike42[10];   //[pho_n]
   Float_t         pho_e2nd[10];   //[pho_n]
   Float_t         pho_e2x5max[10];   //[pho_n]
   Float_t         pho_e2x5right[10];   //[pho_n]
   Float_t         pho_e2x5left[10];   //[pho_n]
   Float_t         pho_e2x5top[10];   //[pho_n]
   Float_t         pho_e2x5bottom[10];   //[pho_n]
   Float_t         pho_eright[10];   //[pho_n]
   Float_t         pho_eleft[10];   //[pho_n]
   Float_t         pho_etop[10];   //[pho_n]
   Float_t         pho_ebottom[10];   //[pho_n]
   Float_t         pho_r19[10];   //[pho_n]
   Float_t         pho_maxoraw[10];   //[pho_n]
   Float_t         pho_cep[10];   //[pho_n]
   Float_t         pho_lambdaratio[10];   //[pho_n]
   Float_t         pho_lambdadivcov[10];   //[pho_n]
   Float_t         pho_etawidth[10];   //[pho_n]
   Float_t         pho_brem[10];   //[pho_n]
   Float_t         pho_smaj[10];   //[pho_n]
   Float_t         pho_e2x2[10];   //[pho_n]
   Float_t         pho_seed_time[10];   //[pho_n]
   Float_t         pho_seed_outoftimechi2[10];   //[pho_n]
   Float_t         pho_seed_chi2[10];   //[pho_n]
   Float_t         pho_seed_recoflag[10];   //[pho_n]
   Float_t         pho_seed_severity[10];   //[pho_n]
   Float_t         pho_ecalsumetconedr04[10];   //[pho_n]
   Float_t         pho_hcalsumetconedr04[10];   //[pho_n]
   Float_t         pho_trksumptsolidconedr04[10];   //[pho_n]
   Float_t         pho_trksumpthollowconedr04[10];   //[pho_n]
   Float_t         pho_ntrksolidconedr04[10];   //[pho_n]
   Float_t         pho_ntrkhollowconedr04[10];   //[pho_n]
   Float_t         pho_ecalsumetconedr03[10];   //[pho_n]
   Float_t         pho_hcalsumetconedr03[10];   //[pho_n]
   Float_t         pho_trksumptsolidconedr03[10];   //[pho_n]
   Float_t         pho_trksumpthollowconedr03[10];   //[pho_n]
   Float_t         pho_ntrksolidconedr03[10];   //[pho_n]
   Float_t         pho_ntrkhollowconedr03[10];   //[pho_n]
   Int_t           pho_haspixseed[10];   //[pho_n]
   Int_t           pho_hasconvtks[10];   //[pho_n]
   Int_t           pho_nconv[10];   //[pho_n]
   Int_t           pho_conv_ntracks[10];   //[pho_n]
   Float_t         pho_conv_pairinvmass[10];   //[pho_n]
   Float_t         pho_conv_paircotthetasep[10];   //[pho_n]
   Float_t         pho_conv_eoverp[10];   //[pho_n]
   Float_t         pho_conv_zofprimvtxfromtrks[10];   //[pho_n]
   Float_t         pho_conv_distofminapproach[10];   //[pho_n]
   Float_t         pho_conv_dphitrksatvtx[10];   //[pho_n]
   Float_t         pho_conv_dphitrksatecal[10];   //[pho_n]
   Float_t         pho_conv_detatrksatecal[10];   //[pho_n]
   Float_t         pho_conv_tk1_d0[10];   //[pho_n]
   Float_t         pho_conv_tk1_pout[10];   //[pho_n]
   Float_t         pho_conv_tk1_pin[10];   //[pho_n]
   Float_t         pho_conv_tk2_d0[10];   //[pho_n]
   Float_t         pho_conv_tk2_pout[10];   //[pho_n]
   Float_t         pho_conv_tk2_pin[10];   //[pho_n]
   Float_t         pho_conv_tk1_dz[10];   //[pho_n]
   Float_t         pho_conv_tk1_dzerr[10];   //[pho_n]
   Int_t           pho_conv_tk1_nh[10];   //[pho_n]
   Float_t         pho_conv_tk2_dz[10];   //[pho_n]
   Float_t         pho_conv_tk2_dzerr[10];   //[pho_n]
   Int_t           pho_conv_tk2_nh[10];   //[pho_n]
   Int_t           pho_conv_ch1ch2[10];   //[pho_n]
   Float_t         pho_conv_chi2[10];   //[pho_n]
   Float_t         pho_conv_chi2_probability[10];   //[pho_n]
   Int_t           pho_conv_validvtx[10];   //[pho_n]
   Int_t           pho_conv_MVALikelihood[10];   //[pho_n]
   TClonesArray    *pho_p4;
   TClonesArray    *pho_calopos;
   TClonesArray    *pho_conv_vtx;
   TClonesArray    *pho_conv_pair_momentum;
   TClonesArray    *pho_conv_refitted_momentum;
   TClonesArray    *pho_conv_vertexcorrected_p4;
   Int_t           pho_scind[10];   //[pho_n]
   Float_t         pho_residCorrEnergy[10];   //[pho_n]
   Float_t         pho_residCorrResn[10];   //[pho_n]
   Float_t         pho_regr_energy[10];   //[pho_n]
   Float_t         pho_regr_energyerr[10];   //[pho_n]
   Int_t           pho_isconv[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral01[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral02[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral03[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral04[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral05[10];   //[pho_n]
   Float_t         pho_pfiso_myneutral06[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton01[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton02[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton03[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton04[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton05[10];   //[pho_n]
   Float_t         pho_pfiso_myphoton06[10];   //[pho_n]
   std::vector<std::vector<float> > *pho_pfiso_mycharged01;
   std::vector<std::vector<float> > *pho_pfiso_mycharged02;
   std::vector<std::vector<float> > *pho_pfiso_mycharged03;
   std::vector<std::vector<float> > *pho_pfiso_mycharged04;
   std::vector<std::vector<float> > *pho_pfiso_mycharged05;
   std::vector<std::vector<float> > *pho_pfiso_mycharged06;
   Int_t           pho_isPFPhoton[10];   //[pho_n]
   Int_t           pho_isPFElectron[10];   //[pho_n]
   Float_t         pho_must[10];   //[pho_n]
   Int_t           pho_mustnc[10];   //[pho_n]
   Float_t         pho_pfpresh1[10];   //[pho_n]
   Float_t         pho_pfpresh2[10];   //[pho_n]
   Float_t         pho_mustenergy[10];   //[pho_n]
   Float_t         pho_mustenergyout[10];   //[pho_n]
   Float_t         pho_pflowE[10];   //[pho_n]
   Float_t         pho_pfdeta[10];   //[pho_n]
   Float_t         pho_pfdphi[10];   //[pho_n]
   Float_t         pho_pfclusrms[10];   //[pho_n]
   Float_t         pho_pfclusrmsmust[10];   //[pho_n]
   Float_t         pho_eseffsixix[10];   //[pho_n]
   Float_t         pho_eseffsiyiy[10];   //[pho_n]
   Float_t         pho_hoe_bc[10];   //[pho_n]
   Int_t           pho_biphi[10];   //[pho_n]
   Int_t           pho_bieta[10];   //[pho_n]
   Float_t         pho_betacry[10];   //[pho_n]
   Float_t         pho_phicry[10];   //[pho_n]
   Int_t           conv_n;
   TClonesArray    *conv_p4;
   Int_t           conv_ntracks[281];   //[conv_n]
   Float_t         conv_pairinvmass[281];   //[conv_n]
   Float_t         conv_paircotthetasep[281];   //[conv_n]
   Float_t         conv_eoverp[281];   //[conv_n]
   Float_t         conv_distofminapproach[281];   //[conv_n]
   Float_t         conv_dphitrksatvtx[281];   //[conv_n]
   Float_t         conv_dphitrksatecal[281];   //[conv_n]
   Float_t         conv_detatrksatecal[281];   //[conv_n]
   Float_t         conv_dxy[281];   //[conv_n]
   Float_t         conv_dz[281];   //[conv_n]
   Float_t         conv_lxy[281];   //[conv_n]
   Float_t         conv_lz[281];   //[conv_n]
   Float_t         conv_zofprimvtxfromtrks[281];   //[conv_n]
   std::vector<std::vector<unsigned short> > *conv_nHitsBeforeVtx;
   Int_t           conv_nSharedHits[281];   //[conv_n]
   Int_t           conv_validvtx[281];   //[conv_n]
   Int_t           conv_MVALikelihood[281];   //[conv_n]
   Float_t         conv_chi2[281];   //[conv_n]
   Float_t         conv_chi2_probability[281];   //[conv_n]
   Float_t         conv_vtx_xErr[281];   //[conv_n]
   Float_t         conv_vtx_yErr[281];   //[conv_n]
   Float_t         conv_vtx_zErr[281];   //[conv_n]
   Float_t         conv_tk1_dz[281];   //[conv_n]
   Float_t         conv_tk2_dz[281];   //[conv_n]
   Float_t         conv_tk1_dzerr[281];   //[conv_n]
   Float_t         conv_tk2_dzerr[281];   //[conv_n]
   Int_t           conv_tk1_nh[281];   //[conv_n]
   Int_t           conv_tk2_nh[281];   //[conv_n]
   Int_t           conv_ch1ch2[281];   //[conv_n]
   Float_t         conv_tk1_d0[281];   //[conv_n]
   Float_t         conv_tk1_pout[281];   //[conv_n]
   Float_t         conv_tk1_pin[281];   //[conv_n]
   Float_t         conv_tk2_d0[281];   //[conv_n]
   Float_t         conv_tk2_pout[281];   //[conv_n]
   Float_t         conv_tk2_pin[281];   //[conv_n]
   TClonesArray    *conv_vtx;
   Float_t         conv_tk1_pterr[281];   //[conv_n]
   Float_t         conv_tk2_pterr[281];   //[conv_n]
   Float_t         conv_tk1_etaerr[281];   //[conv_n]
   Float_t         conv_tk2_etaerr[281];   //[conv_n]
   Float_t         conv_tk1_thetaerr[281];   //[conv_n]
   Float_t         conv_tk2_thetaerr[281];   //[conv_n]
   Float_t         conv_tk1_phierr[281];   //[conv_n]
   Float_t         conv_tk2_phierr[281];   //[conv_n]
   Float_t         conv_tk1_lambdaerr[281];   //[conv_n]
   Float_t         conv_tk2_lambdaerr[281];   //[conv_n]
   TClonesArray    *conv_pair_momentum;
   TClonesArray    *conv_refitted_momentum;
   TClonesArray    *conv_singleleg_momentum;
   Int_t           process_id;
   Float_t         weight;
   Float_t         pthat;
   Int_t           gv_n;
   TClonesArray    *gv_pos;
   Int_t           pu_n;
   std::vector<float>   *pu_zpos;
   std::vector<float>   *pu_sumpt_lowpt;
   std::vector<float>   *pu_sumpt_highpt;
   std::vector<int>     *pu_ntrks_lowpt;
   std::vector<int>     *pu_ntrks_highpt;
   Int_t           vtx_std_n;
   Int_t           vtx_std_ntks[44];   //[vtx_std_n]
   Float_t         vtx_std_x2dof[44];   //[vtx_std_n]
   TClonesArray    *vtx_std_xyz;
   TClonesArray    *vtx_std_dxdydz;
   Float_t         vtx_std_ndof[44];   //[vtx_std_n]
   std::vector<std::vector<unsigned short> > *vtx_std_tkind;
   Float_t         rho_algo1;
   Float_t         rho_algo2;
   Float_t         rho_algo3;
   TClonesArray    *bs_xyz;
   Float_t         bs_sigmaZ;
   Float_t         bs_x0Error;
   Float_t         bs_y0Error;
   Float_t         bs_z0Error;
   Float_t         bs_sigmaZ0Error;
   Float_t         met_tcmet;
   Float_t         met_phi_tcmet;
   std::vector<unsigned short> *hlt_bit;
   Int_t           hlt_n;
   std::vector<std::vector<unsigned short> > *hlt_candpath;
   std::vector<std::string>  *hlt_path_names_HLT;
   TClonesArray    *hlt_p4;
   Int_t           jet_algoPF1_n;
   Float_t         jet_algoPF1_erescale[114];   //[jet_algoPF1_n]
   TClonesArray    *jet_algoPF1_p4;
   Float_t         jet_algoPF1_beta[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_betaStar[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_betaStarClassic[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_dR2Mean[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_dRMean[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_dZ[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_frac01[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_frac02[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_frac03[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_frac04[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_frac05[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_full_mva[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_simple_mva[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_nCharged[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_nNeutrals[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_full_wp_level[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_simple_wp_level[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_cutbased_wp_level[114];   //[jet_algoPF1_n]
   Bool_t          jet_algoPF1_pfloose[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_area[114];   //[jet_algoPF1_n]
   Int_t           jet_algoPF1_nvtx;
   std::vector<std::vector<float> > *jet_algoPF1_beta_ext;
   std::vector<std::vector<float> > *jet_algoPF1_betaStar_ext;
   std::vector<std::vector<float> > *jet_algoPF1_betaStarClassic_ext;
   std::vector<std::vector<int> > *jet_algoPF1_full_wp_level_ext;
   std::vector<std::vector<int> > *jet_algoPF1_simple_wp_level_ext;
   std::vector<std::vector<int> > *jet_algoPF1_cutbased_wp_level_ext;
   std::vector<std::vector<float> > *jet_algoPF1_full_mva_ext;
   std::vector<std::vector<float> > *jet_algoPF1_simple_mva_ext;
   Int_t           jet_algoPF3_n;
   Float_t         jet_algoPF3_erescale[108];   //[jet_algoPF3_n]
   TClonesArray    *jet_algoPF3_p4;
   Float_t         jet_algoPF3_beta[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_betaStar[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_betaStarClassic[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_dR2Mean[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_dRMean[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_dZ[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_frac01[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_frac02[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_frac03[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_frac04[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_frac05[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_full_mva[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_simple_mva[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_nCharged[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_nNeutrals[108];   //[jet_algoPF3_n]
   Int_t           jet_algoPF3_full_wp_level[108];   //[jet_algoPF3_n]
   Int_t           jet_algoPF3_simple_wp_level[108];   //[jet_algoPF3_n]
   Int_t           jet_algoPF3_cutbased_wp_level[108];   //[jet_algoPF3_n]
   Bool_t          jet_algoPF3_pfloose[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_area[108];   //[jet_algoPF3_n]
   Int_t           jet_algoPF3_nvtx;
   std::vector<std::vector<float> > *jet_algoPF3_beta_ext;
   std::vector<std::vector<float> > *jet_algoPF3_betaStar_ext;
   std::vector<std::vector<float> > *jet_algoPF3_betaStarClassic_ext;
   std::vector<std::vector<int> > *jet_algoPF3_full_wp_level_ext;
   std::vector<std::vector<int> > *jet_algoPF3_simple_wp_level_ext;
   std::vector<std::vector<int> > *jet_algoPF3_cutbased_wp_level_ext;
   std::vector<std::vector<float> > *jet_algoPF3_full_mva_ext;
   std::vector<std::vector<float> > *jet_algoPF3_simple_mva_ext;
   Float_t         pho_pfRawEnergy[10];   //[pho_n]
   Float_t         pho_pfe2x2[10];   //[pho_n]
   Float_t         pho_pfe3x3[10];   //[pho_n]
   Float_t         pho_pfe5x5[10];   //[pho_n]
   Float_t         pho_pfsieie[10];   //[pho_n]
   Float_t         pho_pfsieip[10];   //[pho_n]
   Float_t         pho_pfsipip[10];   //[pho_n]
   Float_t         pho_pfemaxxtal[10];   //[pho_n]
   Float_t         pho_pfe2nd[10];   //[pho_n]
   std::vector<std::vector<float> > *vtx_std_mva;
   std::vector<float>   *vtx_std_vertexz;
   std::vector<float>   *vtx_std_nconv;
   std::vector<float>   *vtx_std_nlegs;
   std::vector<std::vector<float> > *vtx_std_pulltoconv;
   std::vector<std::vector<float> > *vtx_std_limpulltoconv;
   std::vector<std::vector<float> > *vtx_std_diphopt;
   std::vector<std::vector<float> > *vtx_std_diphopx;
   std::vector<std::vector<float> > *vtx_std_diphopy;
   std::vector<std::vector<float> > *vtx_std_nch;
   std::vector<std::vector<float> > *vtx_std_ptmax;
   std::vector<std::vector<float> > *vtx_std_sumpt;
   std::vector<std::vector<float> > *vtx_std_ptvtx;
   std::vector<std::vector<float> > *vtx_std_pxvtx;
   std::vector<std::vector<float> > *vtx_std_pyvtx;
   std::vector<std::vector<float> > *vtx_std_acosA;
   std::vector<std::vector<float> > *vtx_std_ptasym;
   std::vector<std::vector<float> > *vtx_std_ptbal;
   std::vector<std::vector<float> > *vtx_std_nchthr;
   std::vector<std::vector<float> > *vtx_std_ptmax3;
   std::vector<std::vector<float> > *vtx_std_thrust;
   std::vector<std::vector<float> > *vtx_std_sumweight;
   std::vector<std::vector<float> > *vtx_std_sumpt2;
   std::vector<std::vector<float> > *vtx_std_ptratio;
   std::vector<std::vector<float> > *vtx_std_pzasym;
   std::vector<std::vector<float> > *vtx_std_spher;
   std::vector<std::vector<float> > *vtx_std_aplan;
   std::vector<std::vector<float> > *vtx_std_sumpr;
   std::vector<std::vector<float> > *vtx_std_sumawy;
   std::vector<std::vector<float> > *vtx_std_sumtrv;
   std::vector<std::vector<float> > *vtx_std_sumtwd;
   std::vector<std::vector<float> > *vtx_std_awytwdasym;
   Int_t           vtx_std_pho1;
   Int_t           vtx_std_pho2;
   std::vector<int>     *pho_matchingConv;
   std::vector<float>   *vtx_std_evt_mva;
   std::vector<std::vector<int> > *vtx_std_ranked_list;
   Int_t           vtx_std_sel;
   std::vector<std::vector<float> > *pho_tkiso_recvtx_030_002_0000_10_01;
   Float_t         pho_tkiso_badvtx_040_002_0000_10_01[10];   //[pho_n]
   Int_t           pho_tkiso_badvtx_id[10];   //[pho_n]
   Float_t         pho_pfiso_charged_badvtx_04[10];   //[pho_n]
   Int_t           pho_pfiso_charged_badvtx_id[10];   //[pho_n]
   std::vector<std::vector<float> > *pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01;
   Float_t         pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01[10];   //[pho_n]
   Int_t           pho_ZeeVal_tkiso_badvtx_id[10];   //[pho_n]
   std::vector<std::vector<float> > *pho_mitmva;
   Float_t         pho_drtotk_25_99[10];   //[pho_n]
   Int_t           dipho_n;
   Int_t           dipho_leadind[10];   //[dipho_n]
   Int_t           dipho_subleadind[10];   //[dipho_n]
   Int_t           dipho_vtxind[10];   //[dipho_n]
   Float_t         dipho_sumpt[10];   //[dipho_n]
   std::vector<std::vector<short> > *pho_cic6cutlevel_lead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic6passcuts_lead;
   std::vector<std::vector<short> > *pho_cic6cutlevel_sublead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic6passcuts_sublead;
   std::vector<std::vector<short> > *pho_cic4cutlevel_lead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4passcuts_lead;
   std::vector<std::vector<short> > *pho_cic4cutlevel_sublead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4passcuts_sublead;
   std::vector<std::vector<short> > *pho_cic4pfcutlevel_lead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4pfpasscuts_lead;
   std::vector<std::vector<short> > *pho_cic4pfcutlevel_sublead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4pfpasscuts_sublead;
   Bool_t          pho_genmatched[10];   //[pho_n]
   Float_t         pho_regr_energy_otf[10];   //[pho_n]
   Float_t         pho_regr_energyerr_otf[10];   //[pho_n]
   Bool_t          jet_algoPF1_genMatched[114];   //[jet_algoPF1_n]
   Bool_t          jet_algoPF1_vbfMatched[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_genPt[114];   //[jet_algoPF1_n]
   Float_t         jet_algoPF1_genDr[114];   //[jet_algoPF1_n]
   Bool_t          jet_algoPF3_genMatched[108];   //[jet_algoPF3_n]
   Bool_t          jet_algoPF3_vbfMatched[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_genPt[108];   //[jet_algoPF3_n]
   Float_t         jet_algoPF3_genDr[108];   //[jet_algoPF3_n]
   Float_t         shiftMET_pt[10];   //[dipho_n]
   Float_t         shiftMET_phi[10];   //[dipho_n]
   Float_t         smearMET_pt[10];   //[dipho_n]
   Float_t         smearMET_phi[10];   //[dipho_n]
   Float_t         shiftsmearMET_pt[10];   //[dipho_n]
   Float_t         shiftsmearMET_phi[10];   //[dipho_n]
   Float_t         shiftscaleMET_pt[10];   //[dipho_n]
   Float_t         shiftscaleMET_phi[10];   //[dipho_n]
   Float_t         shiftMET_eta[10];   //[dipho_n]
   Float_t         shiftMET_e[10];   //[dipho_n]
   Float_t         shiftscaleMET_eta[10];   //[dipho_n]
   Float_t         shiftscaleMET_e[10];   //[dipho_n]
   Int_t           gh_gen2reco1;
   Int_t           gh_gen2reco2;
   Int_t           gh_vbfq1_pdgid;
   Int_t           gh_vbfq2_pdgid;
   Int_t           gh_vh_pdgid;
   Int_t           gh_vh1_pdgid;
   Int_t           gh_vh2_pdgid;
   TClonesArray    *gh_higgs_p4;
   TClonesArray    *gh_pho1_p4;
   TClonesArray    *gh_pho2_p4;
   TClonesArray    *gh_vbfq1_p4;
   TClonesArray    *gh_vbfq2_p4;
   TClonesArray    *gh_vh1_p4;
   TClonesArray    *gh_vh2_p4;
   Int_t           mu_glo_hasgsftrack[39];   //[mu_glo_n]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_run;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_sc_p4;   //!
   TBranch        *b_sc_xyz;   //!
   TBranch        *b_sc_n;   //!
   TBranch        *b_sc_sphi;   //!
   TBranch        *b_sc_seta;   //!
   TBranch        *b_sc_raw;   //!
   TBranch        *b_sc_pre;   //!
   TBranch        *b_sc_bcseedind;   //!
   TBranch        *b_sc_nbc;   //!
   TBranch        *b_ecalhit_n;   //!
   TBranch        *b_ecalhit_detid;   //!
   TBranch        *b_ecalhit_p4;   //!
   TBranch        *b_bc_n;   //!
   TBranch        *b_bc_p4;   //!
   TBranch        *b_bc_s25;   //!
   TBranch        *b_tk_n;   //!
   TBranch        *b_met_pfmet;   //!
   TBranch        *b_met_phi_pfmet;   //!
   TBranch        *b_met_sumet_pfmet;   //!
   TBranch        *b_met_mEtSig_pfmet;   //!
   TBranch        *b_met_significance_pfmet;   //!
   TBranch        *b_met_pfmetType1;   //!
   TBranch        *b_met_phi_pfmetType1;   //!
   TBranch        *b_met_sumet_pfmetType1;   //!
   TBranch        *b_met_mEtSig_pfmetType1;   //!
   TBranch        *b_met_significance_pfmetType1;   //!
   TBranch        *b_el_std_n;   //!
   TBranch        *b_el_std_p4;   //!
   TBranch        *b_el_std_sc;   //!
   TBranch        *b_el_std_scind;   //!
   TBranch        *b_el_std_sieie;   //!
   TBranch        *b_el_std_sieiesc;   //!
   TBranch        *b_el_std_hp_expin;   //!
   TBranch        *b_el_std_detain;   //!
   TBranch        *b_el_std_dphiin;   //!
   TBranch        *b_el_std_eopin;   //!
   TBranch        *b_el_std_fbrem;   //!
   TBranch        *b_el_std_ip_gsf;   //!
   TBranch        *b_el_std_hoe;   //!
   TBranch        *b_el_std_eseedopin;   //!
   TBranch        *b_el_std_tkiso03;   //!
   TBranch        *b_el_std_ecaliso03;   //!
   TBranch        *b_el_std_ecaliso04;   //!
   TBranch        *b_el_std_hcaliso03;   //!
   TBranch        *b_el_std_hcaliso04;   //!
   TBranch        *b_el_std_dcot;   //!
   TBranch        *b_el_std_dist;   //!
   TBranch        *b_el_std_tkind;   //!
   TBranch        *b_el_std_d0;   //!
   TBranch        *b_el_std_z0;   //!
   TBranch        *b_el_std_posvtx;   //!
   TBranch        *b_el_std_D0Vtx;   //!
   TBranch        *b_el_std_DZVtx;   //!
   TBranch        *b_el_std_ip_ctf;   //!
   TBranch        *b_el_std_tkdrv;   //!
   TBranch        *b_el_std_ecaldrv;   //!
   TBranch        *b_el_std_pfiso_charged;   //!
   TBranch        *b_el_std_pfiso_neutral;   //!
   TBranch        *b_el_std_pfiso_photon;   //!
   TBranch        *b_el_std_conv_vtxProb;   //!
   TBranch        *b_el_std_conv;   //!
   TBranch        *b_el_std_pin;   //!
   TBranch        *b_el_std_mva_nontrig;   //!
   TBranch        *b_el_std_mva_trig;   //!
   TBranch        *b_el_std_regr_energy;   //!
   TBranch        *b_el_std_regr_energyerr;   //!
   TBranch        *b_mu_glo_n;   //!
   TBranch        *b_mu_glo_p4;   //!
   TBranch        *b_mu_glo_dof;   //!
   TBranch        *b_mu_glo_chi2;   //!
   TBranch        *b_mu_glo_type;   //!
   TBranch        *b_mu_glo_charge;   //!
   TBranch        *b_mu_glo_nmatches;   //!
   TBranch        *b_mu_glo_posvtx;   //!
   TBranch        *b_mu_glo_losthits;   //!
   TBranch        *b_mu_glo_validhits;   //!
   TBranch        *b_mu_glo_innerhits;   //!
   TBranch        *b_mu_glo_pixelhits;   //!
   TBranch        *b_mu_glo_validChmbhits;   //!
   TBranch        *b_mu_glo_tkpterr;   //!
   TBranch        *b_mu_glo_ecaliso03;   //!
   TBranch        *b_mu_glo_hcaliso03;   //!
   TBranch        *b_mu_glo_tkiso03;   //!
   TBranch        *b_mu_glo_dz;   //!
   TBranch        *b_mu_glo_D0Vtx;   //!
   TBranch        *b_mu_glo_DZVtx;   //!
   TBranch        *b_mu_tkLayers;   //!
   TBranch        *b_mu_glo_chhadiso04;   //!
   TBranch        *b_mu_glo_nehadiso04;   //!
   TBranch        *b_mu_glo_photiso04;   //!
   TBranch        *b_mu_dbCorr;   //!
   TBranch        *b_pho_n;   //!
   TBranch        *b_pho_feta;   //!
   TBranch        *b_pho_crackcorr;   //!
   TBranch        *b_pho_localcorr;   //!
   TBranch        *b_pho_isEB;   //!
   TBranch        *b_pho_isEE;   //!
   TBranch        *b_pho_see;   //!
   TBranch        *b_pho_sieie;   //!
   TBranch        *b_pho_sipip;   //!
   TBranch        *b_pho_sieip;   //!
   TBranch        *b_pho_e1x5;   //!
   TBranch        *b_pho_e3x3;   //!
   TBranch        *b_pho_e5x5;   //!
   TBranch        *b_pho_emaxxtal;   //!
   TBranch        *b_pho_hoe;   //!
   TBranch        *b_pho_r1x5;   //!
   TBranch        *b_pho_r2x5;   //!
   TBranch        *b_pho_r9;   //!
   TBranch        *b_pho_isEBGap;   //!
   TBranch        *b_pho_isEEGap;   //!
   TBranch        *b_pho_isEBEEGap;   //!
   TBranch        *b_pho_zernike20;   //!
   TBranch        *b_pho_zernike42;   //!
   TBranch        *b_pho_e2nd;   //!
   TBranch        *b_pho_e2x5max;   //!
   TBranch        *b_pho_e2x5right;   //!
   TBranch        *b_pho_e2x5left;   //!
   TBranch        *b_pho_e2x5top;   //!
   TBranch        *b_pho_e2x5bottom;   //!
   TBranch        *b_pho_eright;   //!
   TBranch        *b_pho_eleft;   //!
   TBranch        *b_pho_etop;   //!
   TBranch        *b_pho_ebottom;   //!
   TBranch        *b_pho_r19;   //!
   TBranch        *b_pho_maxoraw;   //!
   TBranch        *b_pho_cep;   //!
   TBranch        *b_pho_lambdaratio;   //!
   TBranch        *b_pho_lambdadivcov;   //!
   TBranch        *b_pho_etawidth;   //!
   TBranch        *b_pho_brem;   //!
   TBranch        *b_pho_smaj;   //!
   TBranch        *b_pho_e2x2;   //!
   TBranch        *b_pho_seed_time;   //!
   TBranch        *b_pho_seed_outoftimechi2;   //!
   TBranch        *b_pho_seed_chi2;   //!
   TBranch        *b_pho_seed_recoflag;   //!
   TBranch        *b_pho_seed_severity;   //!
   TBranch        *b_pho_ecalsumetconedr04;   //!
   TBranch        *b_pho_hcalsumetconedr04;   //!
   TBranch        *b_pho_trksumptsolidconedr04;   //!
   TBranch        *b_pho_trksumpthollowconedr04;   //!
   TBranch        *b_pho_ntrksolidconedr04;   //!
   TBranch        *b_pho_ntrkhollowconedr04;   //!
   TBranch        *b_pho_ecalsumetconedr03;   //!
   TBranch        *b_pho_hcalsumetconedr03;   //!
   TBranch        *b_pho_trksumptsolidconedr03;   //!
   TBranch        *b_pho_trksumpthollowconedr03;   //!
   TBranch        *b_pho_ntrksolidconedr03;   //!
   TBranch        *b_pho_ntrkhollowconedr03;   //!
   TBranch        *b_pho_haspixseed;   //!
   TBranch        *b_pho_hasconvtks;   //!
   TBranch        *b_pho_nconv;   //!
   TBranch        *b_pho_conv_ntracks;   //!
   TBranch        *b_pho_conv_pairinvmass;   //!
   TBranch        *b_pho_conv_paircotthetasep;   //!
   TBranch        *b_pho_conv_eoverp;   //!
   TBranch        *b_pho_conv_zofprimvtxfromtrks;   //!
   TBranch        *b_pho_conv_distofminapproach;   //!
   TBranch        *b_pho_conv_dphitrksatvtx;   //!
   TBranch        *b_pho_conv_dphitrksatecal;   //!
   TBranch        *b_pho_conv_detatrksatecal;   //!
   TBranch        *b_pho_conv_tk1_d0;   //!
   TBranch        *b_pho_conv_tk1_pout;   //!
   TBranch        *b_pho_conv_tk1_pin;   //!
   TBranch        *b_pho_conv_tk2_d0;   //!
   TBranch        *b_pho_conv_tk2_pout;   //!
   TBranch        *b_pho_conv_tk2_pin;   //!
   TBranch        *b_pho_conv_tk1_dz;   //!
   TBranch        *b_pho_conv_tk1_dzerr;   //!
   TBranch        *b_pho_conv_tk1_nh;   //!
   TBranch        *b_pho_conv_tk2_dz;   //!
   TBranch        *b_pho_conv_tk2_dzerr;   //!
   TBranch        *b_pho_conv_tk2_nh;   //!
   TBranch        *b_pho_conv_ch1ch2;   //!
   TBranch        *b_pho_conv_chi2;   //!
   TBranch        *b_pho_conv_chi2_probability;   //!
   TBranch        *b_pho_conv_validvtx;   //!
   TBranch        *b_pho_conv_MVALikelihood;   //!
   TBranch        *b_pho_p4;   //!
   TBranch        *b_pho_calopos;   //!
   TBranch        *b_pho_conv_vtx;   //!
   TBranch        *b_pho_conv_pair_momentum;   //!
   TBranch        *b_pho_conv_refitted_momentum;   //!
   TBranch        *b_pho_conv_vertexcorrected_p4;   //!
   TBranch        *b_pho_scind;   //!
   TBranch        *b_pho_residCorrEnergy;   //!
   TBranch        *b_pho_residCorrResn;   //!
   TBranch        *b_pho_regr_energy;   //!
   TBranch        *b_pho_regr_energyerr;   //!
   TBranch        *b_pho_isconv;   //!
   TBranch        *b_pho_pfiso_myneutral01;   //!
   TBranch        *b_pho_pfiso_myneutral02;   //!
   TBranch        *b_pho_pfiso_myneutral03;   //!
   TBranch        *b_pho_pfiso_myneutral04;   //!
   TBranch        *b_pho_pfiso_myneutral05;   //!
   TBranch        *b_pho_pfiso_myneutral06;   //!
   TBranch        *b_pho_pfiso_myphoton01;   //!
   TBranch        *b_pho_pfiso_myphoton02;   //!
   TBranch        *b_pho_pfiso_myphoton03;   //!
   TBranch        *b_pho_pfiso_myphoton04;   //!
   TBranch        *b_pho_pfiso_myphoton05;   //!
   TBranch        *b_pho_pfiso_myphoton06;   //!
   TBranch        *b_pho_pfiso_mycharged01;   //!
   TBranch        *b_pho_pfiso_mycharged02;   //!
   TBranch        *b_pho_pfiso_mycharged03;   //!
   TBranch        *b_pho_pfiso_mycharged04;   //!
   TBranch        *b_pho_pfiso_mycharged05;   //!
   TBranch        *b_pho_pfiso_mycharged06;   //!
   TBranch        *b_pho_isPFPhoton;   //!
   TBranch        *b_pho_isPFElectron;   //!
   TBranch        *b_pho_must;   //!
   TBranch        *b_pho_mustnc;   //!
   TBranch        *b_pho_pfpresh1;   //!
   TBranch        *b_pho_pfpresh2;   //!
   TBranch        *b_pho_mustenergy;   //!
   TBranch        *b_pho_mustenergyout;   //!
   TBranch        *b_pho_pflowE;   //!
   TBranch        *b_pho_pfdeta;   //!
   TBranch        *b_pho_pfdphi;   //!
   TBranch        *b_pho_pfclusrms;   //!
   TBranch        *b_pho_pfclusrmsmust;   //!
   TBranch        *b_pho_eseffsixix;   //!
   TBranch        *b_pho_eseffsiyiy;   //!
   TBranch        *b_pho_hoe_bc;   //!
   TBranch        *b_pho_biphi;   //!
   TBranch        *b_pho_bieta;   //!
   TBranch        *b_pho_betacry;   //!
   TBranch        *b_pho_phicry;   //!
   TBranch        *b_conv_n;   //!
   TBranch        *b_conv_p4;   //!
   TBranch        *b_conv_ntracks;   //!
   TBranch        *b_conv_pairinvmass;   //!
   TBranch        *b_conv_paircotthetasep;   //!
   TBranch        *b_conv_eoverp;   //!
   TBranch        *b_conv_distofminapproach;   //!
   TBranch        *b_conv_dphitrksatvtx;   //!
   TBranch        *b_conv_dphitrksatecal;   //!
   TBranch        *b_conv_detatrksatecal;   //!
   TBranch        *b_conv_dxy;   //!
   TBranch        *b_conv_dz;   //!
   TBranch        *b_conv_lxy;   //!
   TBranch        *b_conv_lz;   //!
   TBranch        *b_conv_zofprimvtxfromtrks;   //!
   TBranch        *b_conv_nHitsBeforeVtx;   //!
   TBranch        *b_conv_nSharedHits;   //!
   TBranch        *b_conv_validvtx;   //!
   TBranch        *b_conv_MVALikelihood;   //!
   TBranch        *b_conv_chi2;   //!
   TBranch        *b_conv_chi2_probability;   //!
   TBranch        *b_conv_vtx_xErr;   //!
   TBranch        *b_conv_vtx_yErr;   //!
   TBranch        *b_conv_vtx_zErr;   //!
   TBranch        *b_conv_tk1_dz;   //!
   TBranch        *b_conv_tk2_dz;   //!
   TBranch        *b_conv_tk1_dzerr;   //!
   TBranch        *b_conv_tk2_dzerr;   //!
   TBranch        *b_conv_tk1_nh;   //!
   TBranch        *b_conv_tk2_nh;   //!
   TBranch        *b_conv_ch1ch2;   //!
   TBranch        *b_conv_tk1_d0;   //!
   TBranch        *b_conv_tk1_pout;   //!
   TBranch        *b_conv_tk1_pin;   //!
   TBranch        *b_conv_tk2_d0;   //!
   TBranch        *b_conv_tk2_pout;   //!
   TBranch        *b_conv_tk2_pin;   //!
   TBranch        *b_conv_vtx;   //!
   TBranch        *b_conv_tk1_pterr;   //!
   TBranch        *b_conv_tk2_pterr;   //!
   TBranch        *b_conv_tk1_etaerr;   //!
   TBranch        *b_conv_tk2_etaerr;   //!
   TBranch        *b_conv_tk1_thetaerr;   //!
   TBranch        *b_conv_tk2_thetaerr;   //!
   TBranch        *b_conv_tk1_phierr;   //!
   TBranch        *b_conv_tk2_phierr;   //!
   TBranch        *b_conv_tk1_lambdaerr;   //!
   TBranch        *b_conv_tk2_lambdaerr;   //!
   TBranch        *b_conv_pair_momentum;   //!
   TBranch        *b_conv_refitted_momentum;   //!
   TBranch        *b_conv_singleleg_momentum;   //!
   TBranch        *b_process_id;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_gv_n;   //!
   TBranch        *b_gv_pos;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_zpos;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_vtx_std_n;   //!
   TBranch        *b_vtx_std_ntks;   //!
   TBranch        *b_vtx_std_x2dof;   //!
   TBranch        *b_vtx_std_xyz;   //!
   TBranch        *b_vtx_std_dxdydz;   //!
   TBranch        *b_vtx_std_ndof;   //!
   TBranch        *b_vtx_std_tkind;   //!
   TBranch        *b_rho_algo1;   //!
   TBranch        *b_rho_algo2;   //!
   TBranch        *b_rho_algo3;   //!
   TBranch        *b_bs_xyz;   //!
   TBranch        *b_bs_sigmaZ;   //!
   TBranch        *b_bs_x0Error;   //!
   TBranch        *b_bs_y0Error;   //!
   TBranch        *b_bs_z0Error;   //!
   TBranch        *b_bs_sigmaZ0Error;   //!
   TBranch        *b_met_tcmet;   //!
   TBranch        *b_met_phi_tcmet;   //!
   TBranch        *b_hlt_bit;   //!
   TBranch        *b_hlt_n;   //!
   TBranch        *b_hlt_candpath;   //!
   TBranch        *b_hlt_path_names_HLT;   //!
   TBranch        *b_hlt_p4;   //!
   TBranch        *b_jet_algoPF1_n;   //!
   TBranch        *b_jet_algoPF1_erescale;   //!
   TBranch        *b_jet_algoPF1_p4;   //!
   TBranch        *b_jet_algoPF1_beta;   //!
   TBranch        *b_jet_algoPF1_betaStar;   //!
   TBranch        *b_jet_algoPF1_betaStarClassic;   //!
   TBranch        *b_jet_algoPF1_dR2Mean;   //!
   TBranch        *b_jet_algoPF1_dRMean;   //!
   TBranch        *b_jet_algoPF1_dZ;   //!
   TBranch        *b_jet_algoPF1_frac01;   //!
   TBranch        *b_jet_algoPF1_frac02;   //!
   TBranch        *b_jet_algoPF1_frac03;   //!
   TBranch        *b_jet_algoPF1_frac04;   //!
   TBranch        *b_jet_algoPF1_frac05;   //!
   TBranch        *b_jet_algoPF1_full_mva;   //!
   TBranch        *b_jet_algoPF1_simple_mva;   //!
   TBranch        *b_jet_algoPF1_nCharged;   //!
   TBranch        *b_jet_algoPF1_nNeutrals;   //!
   TBranch        *b_jet_algoPF1_full_wp_level;   //!
   TBranch        *b_jet_algoPF1_simple_wp_level;   //!
   TBranch        *b_jet_algoPF1_cutbased_wp_level;   //!
   TBranch        *b_jet_algoPF1_pfloose;   //!
   TBranch        *b_jet_algoPF1_area;   //!
   TBranch        *b_jet_algoPF1_nvtx;   //!
   TBranch        *b_jet_algoPF1_beta_ext;   //!
   TBranch        *b_jet_algoPF1_betaStar_ext;   //!
   TBranch        *b_jet_algoPF1_betaStarClassic_ext;   //!
   TBranch        *b_jet_algoPF1_full_wp_level_ext;   //!
   TBranch        *b_jet_algoPF1_simple_wp_level_ext;   //!
   TBranch        *b_jet_algoPF1_cutbased_wp_level_ext;   //!
   TBranch        *b_jet_algoPF1_full_mva_ext;   //!
   TBranch        *b_jet_algoPF1_simple_mva_ext;   //!
   TBranch        *b_jet_algoPF3_n;   //!
   TBranch        *b_jet_algoPF3_erescale;   //!
   TBranch        *b_jet_algoPF3_p4;   //!
   TBranch        *b_jet_algoPF3_beta;   //!
   TBranch        *b_jet_algoPF3_betaStar;   //!
   TBranch        *b_jet_algoPF3_betaStarClassic;   //!
   TBranch        *b_jet_algoPF3_dR2Mean;   //!
   TBranch        *b_jet_algoPF3_dRMean;   //!
   TBranch        *b_jet_algoPF3_dZ;   //!
   TBranch        *b_jet_algoPF3_frac01;   //!
   TBranch        *b_jet_algoPF3_frac02;   //!
   TBranch        *b_jet_algoPF3_frac03;   //!
   TBranch        *b_jet_algoPF3_frac04;   //!
   TBranch        *b_jet_algoPF3_frac05;   //!
   TBranch        *b_jet_algoPF3_full_mva;   //!
   TBranch        *b_jet_algoPF3_simple_mva;   //!
   TBranch        *b_jet_algoPF3_nCharged;   //!
   TBranch        *b_jet_algoPF3_nNeutrals;   //!
   TBranch        *b_jet_algoPF3_full_wp_level;   //!
   TBranch        *b_jet_algoPF3_simple_wp_level;   //!
   TBranch        *b_jet_algoPF3_cutbased_wp_level;   //!
   TBranch        *b_jet_algoPF3_pfloose;   //!
   TBranch        *b_jet_algoPF3_area;   //!
   TBranch        *b_jet_algoPF3_nvtx;   //!
   TBranch        *b_jet_algoPF3_beta_ext;   //!
   TBranch        *b_jet_algoPF3_betaStar_ext;   //!
   TBranch        *b_jet_algoPF3_betaStarClassic_ext;   //!
   TBranch        *b_jet_algoPF3_full_wp_level_ext;   //!
   TBranch        *b_jet_algoPF3_simple_wp_level_ext;   //!
   TBranch        *b_jet_algoPF3_cutbased_wp_level_ext;   //!
   TBranch        *b_jet_algoPF3_full_mva_ext;   //!
   TBranch        *b_jet_algoPF3_simple_mva_ext;   //!
   TBranch        *b_pho_pfRawEnergy;   //!
   TBranch        *b_pho_pfe2x2;   //!
   TBranch        *b_pho_pfe3x3;   //!
   TBranch        *b_pho_pfe5x5;   //!
   TBranch        *b_pho_pfsieie;   //!
   TBranch        *b_pho_pfsieip;   //!
   TBranch        *b_pho_pfsipip;   //!
   TBranch        *b_pho_pfemaxxtal;   //!
   TBranch        *b_pho_pfe2nd;   //!
   TBranch        *b_vtx_std_mva;   //!
   TBranch        *b_vtx_std_vertexz;   //!
   TBranch        *b_vtx_std_nconv;   //!
   TBranch        *b_vtx_std_nlegs;   //!
   TBranch        *b_vtx_std_pulltoconv;   //!
   TBranch        *b_vtx_std_limpulltoconv;   //!
   TBranch        *b_vtx_std_diphopt;   //!
   TBranch        *b_vtx_std_diphopx;   //!
   TBranch        *b_vtx_std_diphopy;   //!
   TBranch        *b_vtx_std_nch;   //!
   TBranch        *b_vtx_std_ptmax;   //!
   TBranch        *b_vtx_std_sumpt;   //!
   TBranch        *b_vtx_std_ptvtx;   //!
   TBranch        *b_vtx_std_pxvtx;   //!
   TBranch        *b_vtx_std_pyvtx;   //!
   TBranch        *b_vtx_std_acosA;   //!
   TBranch        *b_vtx_std_ptasym;   //!
   TBranch        *b_vtx_std_ptbal;   //!
   TBranch        *b_vtx_std_nchthr;   //!
   TBranch        *b_vtx_std_ptmax3;   //!
   TBranch        *b_vtx_std_thrust;   //!
   TBranch        *b_vtx_std_sumweight;   //!
   TBranch        *b_vtx_std_sumpt2;   //!
   TBranch        *b_vtx_std_ptratio;   //!
   TBranch        *b_vtx_std_pzasym;   //!
   TBranch        *b_vtx_std_spher;   //!
   TBranch        *b_vtx_std_aplan;   //!
   TBranch        *b_vtx_std_sumpr;   //!
   TBranch        *b_vtx_std_sumawy;   //!
   TBranch        *b_vtx_std_sumtrv;   //!
   TBranch        *b_vtx_std_sumtwd;   //!
   TBranch        *b_vtx_std_awytwdasym;   //!
   TBranch        *b_vtx_std_pho1;   //!
   TBranch        *b_vtx_std_pho2;   //!
   TBranch        *b_pho_matchingConv;   //!
   TBranch        *b_vtx_std_evt_mva;   //!
   TBranch        *b_vtx_std_ranked_list;   //!
   TBranch        *b_vtx_std_sel;   //!
   TBranch        *b_pho_tkiso_recvtx_030_002_0000_10_01;   //!
   TBranch        *b_pho_tkiso_badvtx_040_002_0000_10_01;   //!
   TBranch        *b_pho_tkiso_badvtx_id;   //!
   TBranch        *b_pho_pfiso_charged_badvtx_04;   //!
   TBranch        *b_pho_pfiso_charged_badvtx_id;   //!
   TBranch        *b_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01;   //!
   TBranch        *b_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01;   //!
   TBranch        *b_pho_ZeeVal_tkiso_badvtx_id;   //!
   TBranch        *b_pho_mitmva;   //!
   TBranch        *b_pho_drtotk_25_99;   //!
   TBranch        *b_dipho_n;   //!
   TBranch        *b_dipho_leadind;   //!
   TBranch        *b_dipho_subleadind;   //!
   TBranch        *b_dipho_vtxind;   //!
   TBranch        *b_dipho_sumpt;   //!
   TBranch        *b_pho_cic6cutlevel_lead;   //!
   TBranch        *b_pho_cic6passcuts_lead;   //!
   TBranch        *b_pho_cic6cutlevel_sublead;   //!
   TBranch        *b_pho_cic6passcuts_sublead;   //!
   TBranch        *b_pho_cic4cutlevel_lead;   //!
   TBranch        *b_pho_cic4passcuts_lead;   //!
   TBranch        *b_pho_cic4cutlevel_sublead;   //!
   TBranch        *b_pho_cic4passcuts_sublead;   //!
   TBranch        *b_pho_cic4pfcutlevel_lead;   //!
   TBranch        *b_pho_cic4pfpasscuts_lead;   //!
   TBranch        *b_pho_cic4pfcutlevel_sublead;   //!
   TBranch        *b_pho_cic4pfpasscuts_sublead;   //!
   TBranch        *b_pho_genmatched;   //!
   TBranch        *b_pho_regr_energy_otf;   //!
   TBranch        *b_pho_regr_energyerr_otf;   //!
   TBranch        *b_jet_algoPF1_genMatched;   //!
   TBranch        *b_jet_algoPF1_vbfMatched;   //!
   TBranch        *b_jet_algoPF1_genPt;   //!
   TBranch        *b_jet_algoPF1_genDr;   //!
   TBranch        *b_jet_algoPF3_genMatched;   //!
   TBranch        *b_jet_algoPF3_vbfMatched;   //!
   TBranch        *b_jet_algoPF3_genPt;   //!
   TBranch        *b_jet_algoPF3_genDr;   //!
   TBranch        *b_shiftMET_pt;   //!
   TBranch        *b_shiftMET_phi;   //!
   TBranch        *b_smearMET_pt;   //!
   TBranch        *b_smearMET_phi;   //!
   TBranch        *b_shiftsmearMET_pt;   //!
   TBranch        *b_shiftsmearMET_phi;   //!
   TBranch        *b_shiftscaleMET_pt;   //!
   TBranch        *b_shiftscaleMET_phi;   //!
   TBranch        *b_shiftMET_eta;   //!
   TBranch        *b_shiftMET_e;   //!
   TBranch        *b_shiftscaleMET_eta;   //!
   TBranch        *b_shiftscaleMET_e;   //!
   TBranch        *b_gh_gen2reco1;   //!
   TBranch        *b_gh_gen2reco2;   //!
   TBranch        *b_gh_vbfq1_pdgid;   //!
   TBranch        *b_gh_vbfq2_pdgid;   //!
   TBranch        *b_gh_vh_pdgid;   //!
   TBranch        *b_gh_vh1_pdgid;   //!
   TBranch        *b_gh_vh2_pdgid;   //!
   TBranch        *b_gh_higgs_p4;   //!
   TBranch        *b_gh_pho1_p4;   //!
   TBranch        *b_gh_pho2_p4;   //!
   TBranch        *b_gh_vbfq1_p4;   //!
   TBranch        *b_gh_vbfq2_p4;   //!
   TBranch        *b_gh_vh1_p4;   //!
   TBranch        *b_gh_vh2_p4;   //!
   TBranch        *b_mu_glo_hasgsftrack;   //!

   GlobeReader(TTree *tree=0);
   virtual ~GlobeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GlobeReader_cxx
GlobeReader::GlobeReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/moriond2013_reduction_v1/data/DoublePhoton_Run2012B-13Jul2012-v1/DoublePhoton_Run2012B-13Jul2012-v1_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/moriond2013_reduction_v1/data/DoublePhoton_Run2012B-13Jul2012-v1/DoublePhoton_Run2012B-13Jul2012-v1_0.root");
      }
      f->GetObject("event",tree);

   }
   Init(tree);
}

GlobeReader::~GlobeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GlobeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GlobeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GlobeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   sc_p4 = 0;
   sc_xyz = 0;
   ecalhit_p4 = 0;
   bc_p4 = 0;
   el_std_p4 = 0;
   el_std_sc = 0;
   el_std_posvtx = 0;
   mu_glo_p4 = 0;
   mu_glo_posvtx = 0;
   pho_p4 = 0;
   pho_calopos = 0;
   pho_conv_vtx = 0;
   pho_conv_pair_momentum = 0;
   pho_conv_refitted_momentum = 0;
   pho_conv_vertexcorrected_p4 = 0;
   pho_pfiso_mycharged01 = 0;
   pho_pfiso_mycharged02 = 0;
   pho_pfiso_mycharged03 = 0;
   pho_pfiso_mycharged04 = 0;
   pho_pfiso_mycharged05 = 0;
   pho_pfiso_mycharged06 = 0;
   conv_p4 = 0;
   conv_nHitsBeforeVtx = 0;
   conv_vtx = 0;
   conv_pair_momentum = 0;
   conv_refitted_momentum = 0;
   conv_singleleg_momentum = 0;
   gv_pos = 0;
   pu_zpos = 0;
   pu_sumpt_lowpt = 0;
   pu_sumpt_highpt = 0;
   pu_ntrks_lowpt = 0;
   pu_ntrks_highpt = 0;
   vtx_std_xyz = 0;
   vtx_std_dxdydz = 0;
   vtx_std_tkind = 0;
   bs_xyz = 0;
   hlt_bit = 0;
   hlt_candpath = 0;
   hlt_path_names_HLT = 0;
   hlt_p4 = 0;
   jet_algoPF1_p4 = 0;
   jet_algoPF1_beta_ext = 0;
   jet_algoPF1_betaStar_ext = 0;
   jet_algoPF1_betaStarClassic_ext = 0;
   jet_algoPF1_full_wp_level_ext = 0;
   jet_algoPF1_simple_wp_level_ext = 0;
   jet_algoPF1_cutbased_wp_level_ext = 0;
   jet_algoPF1_full_mva_ext = 0;
   jet_algoPF1_simple_mva_ext = 0;
   jet_algoPF3_p4 = 0;
   jet_algoPF3_beta_ext = 0;
   jet_algoPF3_betaStar_ext = 0;
   jet_algoPF3_betaStarClassic_ext = 0;
   jet_algoPF3_full_wp_level_ext = 0;
   jet_algoPF3_simple_wp_level_ext = 0;
   jet_algoPF3_cutbased_wp_level_ext = 0;
   jet_algoPF3_full_mva_ext = 0;
   jet_algoPF3_simple_mva_ext = 0;
   vtx_std_mva = 0;
   vtx_std_vertexz = 0;
   vtx_std_nconv = 0;
   vtx_std_nlegs = 0;
   vtx_std_pulltoconv = 0;
   vtx_std_limpulltoconv = 0;
   vtx_std_diphopt = 0;
   vtx_std_diphopx = 0;
   vtx_std_diphopy = 0;
   vtx_std_nch = 0;
   vtx_std_ptmax = 0;
   vtx_std_sumpt = 0;
   vtx_std_ptvtx = 0;
   vtx_std_pxvtx = 0;
   vtx_std_pyvtx = 0;
   vtx_std_acosA = 0;
   vtx_std_ptasym = 0;
   vtx_std_ptbal = 0;
   vtx_std_nchthr = 0;
   vtx_std_ptmax3 = 0;
   vtx_std_thrust = 0;
   vtx_std_sumweight = 0;
   vtx_std_sumpt2 = 0;
   vtx_std_ptratio = 0;
   vtx_std_pzasym = 0;
   vtx_std_spher = 0;
   vtx_std_aplan = 0;
   vtx_std_sumpr = 0;
   vtx_std_sumawy = 0;
   vtx_std_sumtrv = 0;
   vtx_std_sumtwd = 0;
   vtx_std_awytwdasym = 0;
   pho_matchingConv = 0;
   vtx_std_evt_mva = 0;
   vtx_std_ranked_list = 0;
   pho_tkiso_recvtx_030_002_0000_10_01 = 0;
   pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01 = 0;
   pho_mitmva = 0;
   pho_cic6cutlevel_lead = 0;
   pho_cic6passcuts_lead = 0;
   pho_cic6cutlevel_sublead = 0;
   pho_cic6passcuts_sublead = 0;
   pho_cic4cutlevel_lead = 0;
   pho_cic4passcuts_lead = 0;
   pho_cic4cutlevel_sublead = 0;
   pho_cic4passcuts_sublead = 0;
   pho_cic4pfcutlevel_lead = 0;
   pho_cic4pfpasscuts_lead = 0;
   pho_cic4pfcutlevel_sublead = 0;
   pho_cic4pfpasscuts_sublead = 0;
   gh_higgs_p4 = 0;
   gh_pho1_p4 = 0;
   gh_pho2_p4 = 0;
   gh_vbfq1_p4 = 0;
   gh_vbfq2_p4 = 0;
   gh_vh1_p4 = 0;
   gh_vh2_p4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("sc_p4", &sc_p4, &b_sc_p4);
   fChain->SetBranchAddress("sc_xyz", &sc_xyz, &b_sc_xyz);
   fChain->SetBranchAddress("sc_n", &sc_n, &b_sc_n);
   fChain->SetBranchAddress("sc_sphi", sc_sphi, &b_sc_sphi);
   fChain->SetBranchAddress("sc_seta", sc_seta, &b_sc_seta);
   fChain->SetBranchAddress("sc_raw", sc_raw, &b_sc_raw);
   fChain->SetBranchAddress("sc_pre", sc_pre, &b_sc_pre);
   fChain->SetBranchAddress("sc_bcseedind", sc_bcseedind, &b_sc_bcseedind);
   fChain->SetBranchAddress("sc_nbc", sc_nbc, &b_sc_nbc);
   fChain->SetBranchAddress("ecalhit_n", &ecalhit_n, &b_ecalhit_n);
   fChain->SetBranchAddress("ecalhit_detid", ecalhit_detid, &b_ecalhit_detid);
   fChain->SetBranchAddress("ecalhit_p4", &ecalhit_p4, &b_ecalhit_p4);
   fChain->SetBranchAddress("bc_n", &bc_n, &b_bc_n);
   fChain->SetBranchAddress("bc_p4", &bc_p4, &b_bc_p4);
   fChain->SetBranchAddress("bc_s25", bc_s25, &b_bc_s25);
   fChain->SetBranchAddress("tk_n", &tk_n, &b_tk_n);
   fChain->SetBranchAddress("met_pfmet", &met_pfmet, &b_met_pfmet);
   fChain->SetBranchAddress("met_phi_pfmet", &met_phi_pfmet, &b_met_phi_pfmet);
   fChain->SetBranchAddress("met_sumet_pfmet", &met_sumet_pfmet, &b_met_sumet_pfmet);
   fChain->SetBranchAddress("met_mEtSig_pfmet", &met_mEtSig_pfmet, &b_met_mEtSig_pfmet);
   fChain->SetBranchAddress("met_significance_pfmet", &met_significance_pfmet, &b_met_significance_pfmet);
   fChain->SetBranchAddress("met_pfmetType1", &met_pfmetType1, &b_met_pfmetType1);
   fChain->SetBranchAddress("met_phi_pfmetType1", &met_phi_pfmetType1, &b_met_phi_pfmetType1);
   fChain->SetBranchAddress("met_sumet_pfmetType1", &met_sumet_pfmetType1, &b_met_sumet_pfmetType1);
   fChain->SetBranchAddress("met_mEtSig_pfmetType1", &met_mEtSig_pfmetType1, &b_met_mEtSig_pfmetType1);
   fChain->SetBranchAddress("met_significance_pfmetType1", &met_significance_pfmetType1, &b_met_significance_pfmetType1);
   fChain->SetBranchAddress("el_std_n", &el_std_n, &b_el_std_n);
   fChain->SetBranchAddress("el_std_p4", &el_std_p4, &b_el_std_p4);
   fChain->SetBranchAddress("el_std_sc", &el_std_sc, &b_el_std_sc);
   fChain->SetBranchAddress("el_std_scind", el_std_scind, &b_el_std_scind);
   fChain->SetBranchAddress("el_std_sieie", el_std_sieie, &b_el_std_sieie);
   fChain->SetBranchAddress("el_std_sieiesc", el_std_sieiesc, &b_el_std_sieiesc);
   fChain->SetBranchAddress("el_std_hp_expin", el_std_hp_expin, &b_el_std_hp_expin);
   fChain->SetBranchAddress("el_std_detain", el_std_detain, &b_el_std_detain);
   fChain->SetBranchAddress("el_std_dphiin", el_std_dphiin, &b_el_std_dphiin);
   fChain->SetBranchAddress("el_std_eopin", el_std_eopin, &b_el_std_eopin);
   fChain->SetBranchAddress("el_std_fbrem", el_std_fbrem, &b_el_std_fbrem);
   fChain->SetBranchAddress("el_std_ip_gsf", el_std_ip_gsf, &b_el_std_ip_gsf);
   fChain->SetBranchAddress("el_std_hoe", el_std_hoe, &b_el_std_hoe);
   fChain->SetBranchAddress("el_std_eseedopin", el_std_eseedopin, &b_el_std_eseedopin);
   fChain->SetBranchAddress("el_std_tkiso03", el_std_tkiso03, &b_el_std_tkiso03);
   fChain->SetBranchAddress("el_std_ecaliso03", el_std_ecaliso03, &b_el_std_ecaliso03);
   fChain->SetBranchAddress("el_std_ecaliso04", el_std_ecaliso04, &b_el_std_ecaliso04);
   fChain->SetBranchAddress("el_std_hcaliso03", el_std_hcaliso03, &b_el_std_hcaliso03);
   fChain->SetBranchAddress("el_std_hcaliso04", el_std_hcaliso04, &b_el_std_hcaliso04);
   fChain->SetBranchAddress("el_std_dcot", el_std_dcot, &b_el_std_dcot);
   fChain->SetBranchAddress("el_std_dist", el_std_dist, &b_el_std_dist);
   fChain->SetBranchAddress("el_std_tkind", el_std_tkind, &b_el_std_tkind);
   fChain->SetBranchAddress("el_std_d0", el_std_d0, &b_el_std_d0);
   fChain->SetBranchAddress("el_std_z0", el_std_z0, &b_el_std_z0);
   fChain->SetBranchAddress("el_std_posvtx", &el_std_posvtx, &b_el_std_posvtx);
   fChain->SetBranchAddress("el_std_D0Vtx", el_std_D0Vtx, &b_el_std_D0Vtx);
   fChain->SetBranchAddress("el_std_DZVtx", el_std_DZVtx, &b_el_std_DZVtx);
   fChain->SetBranchAddress("el_std_ip_ctf", el_std_ip_ctf, &b_el_std_ip_ctf);
   fChain->SetBranchAddress("el_std_tkdrv", el_std_tkdrv, &b_el_std_tkdrv);
   fChain->SetBranchAddress("el_std_ecaldrv", el_std_ecaldrv, &b_el_std_ecaldrv);
   fChain->SetBranchAddress("el_std_pfiso_charged", el_std_pfiso_charged, &b_el_std_pfiso_charged);
   fChain->SetBranchAddress("el_std_pfiso_neutral", el_std_pfiso_neutral, &b_el_std_pfiso_neutral);
   fChain->SetBranchAddress("el_std_pfiso_photon", el_std_pfiso_photon, &b_el_std_pfiso_photon);
   fChain->SetBranchAddress("el_std_conv_vtxProb", el_std_conv_vtxProb, &b_el_std_conv_vtxProb);
   fChain->SetBranchAddress("el_std_conv", el_std_conv, &b_el_std_conv);
   fChain->SetBranchAddress("el_std_pin", el_std_pin, &b_el_std_pin);
   fChain->SetBranchAddress("el_std_mva_nontrig", el_std_mva_nontrig, &b_el_std_mva_nontrig);
   fChain->SetBranchAddress("el_std_mva_trig", el_std_mva_trig, &b_el_std_mva_trig);
   fChain->SetBranchAddress("el_std_regr_energy", el_std_regr_energy, &b_el_std_regr_energy);
   fChain->SetBranchAddress("el_std_regr_energyerr", el_std_regr_energyerr, &b_el_std_regr_energyerr);
   fChain->SetBranchAddress("mu_glo_n", &mu_glo_n, &b_mu_glo_n);
   fChain->SetBranchAddress("mu_glo_p4", &mu_glo_p4, &b_mu_glo_p4);
   fChain->SetBranchAddress("mu_glo_dof", mu_glo_dof, &b_mu_glo_dof);
   fChain->SetBranchAddress("mu_glo_chi2", mu_glo_chi2, &b_mu_glo_chi2);
   fChain->SetBranchAddress("mu_glo_type", mu_glo_type, &b_mu_glo_type);
   fChain->SetBranchAddress("mu_glo_charge", mu_glo_charge, &b_mu_glo_charge);
   fChain->SetBranchAddress("mu_glo_nmatches", mu_glo_nmatches, &b_mu_glo_nmatches);
   fChain->SetBranchAddress("mu_glo_posvtx", &mu_glo_posvtx, &b_mu_glo_posvtx);
   fChain->SetBranchAddress("mu_glo_losthits", mu_glo_losthits, &b_mu_glo_losthits);
   fChain->SetBranchAddress("mu_glo_validhits", mu_glo_validhits, &b_mu_glo_validhits);
   fChain->SetBranchAddress("mu_glo_innerhits", mu_glo_innerhits, &b_mu_glo_innerhits);
   fChain->SetBranchAddress("mu_glo_pixelhits", mu_glo_pixelhits, &b_mu_glo_pixelhits);
   fChain->SetBranchAddress("mu_glo_validChmbhits", mu_glo_validChmbhits, &b_mu_glo_validChmbhits);
   fChain->SetBranchAddress("mu_glo_tkpterr", mu_glo_tkpterr, &b_mu_glo_tkpterr);
   fChain->SetBranchAddress("mu_glo_ecaliso03", mu_glo_ecaliso03, &b_mu_glo_ecaliso03);
   fChain->SetBranchAddress("mu_glo_hcaliso03", mu_glo_hcaliso03, &b_mu_glo_hcaliso03);
   fChain->SetBranchAddress("mu_glo_tkiso03", mu_glo_tkiso03, &b_mu_glo_tkiso03);
   fChain->SetBranchAddress("mu_glo_dz", mu_glo_dz, &b_mu_glo_dz);
   fChain->SetBranchAddress("mu_glo_D0Vtx", mu_glo_D0Vtx, &b_mu_glo_D0Vtx);
   fChain->SetBranchAddress("mu_glo_DZVtx", mu_glo_DZVtx, &b_mu_glo_DZVtx);
   fChain->SetBranchAddress("mu_tkLayers", mu_tkLayers, &b_mu_tkLayers);
   fChain->SetBranchAddress("mu_glo_chhadiso04", mu_glo_chhadiso04, &b_mu_glo_chhadiso04);
   fChain->SetBranchAddress("mu_glo_nehadiso04", mu_glo_nehadiso04, &b_mu_glo_nehadiso04);
   fChain->SetBranchAddress("mu_glo_photiso04", mu_glo_photiso04, &b_mu_glo_photiso04);
   fChain->SetBranchAddress("mu_dbCorr", mu_dbCorr, &b_mu_dbCorr);
   fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
   fChain->SetBranchAddress("pho_feta", pho_feta, &b_pho_feta);
   fChain->SetBranchAddress("pho_crackcorr", pho_crackcorr, &b_pho_crackcorr);
   fChain->SetBranchAddress("pho_localcorr", pho_localcorr, &b_pho_localcorr);
   fChain->SetBranchAddress("pho_isEB", pho_isEB, &b_pho_isEB);
   fChain->SetBranchAddress("pho_isEE", pho_isEE, &b_pho_isEE);
   fChain->SetBranchAddress("pho_see", pho_see, &b_pho_see);
   fChain->SetBranchAddress("pho_sieie", pho_sieie, &b_pho_sieie);
   fChain->SetBranchAddress("pho_sipip", pho_sipip, &b_pho_sipip);
   fChain->SetBranchAddress("pho_sieip", pho_sieip, &b_pho_sieip);
   fChain->SetBranchAddress("pho_e1x5", pho_e1x5, &b_pho_e1x5);
   fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
   fChain->SetBranchAddress("pho_e5x5", pho_e5x5, &b_pho_e5x5);
   fChain->SetBranchAddress("pho_emaxxtal", pho_emaxxtal, &b_pho_emaxxtal);
   fChain->SetBranchAddress("pho_hoe", pho_hoe, &b_pho_hoe);
   fChain->SetBranchAddress("pho_r1x5", pho_r1x5, &b_pho_r1x5);
   fChain->SetBranchAddress("pho_r2x5", pho_r2x5, &b_pho_r2x5);
   fChain->SetBranchAddress("pho_r9", pho_r9, &b_pho_r9);
   fChain->SetBranchAddress("pho_isEBGap", pho_isEBGap, &b_pho_isEBGap);
   fChain->SetBranchAddress("pho_isEEGap", pho_isEEGap, &b_pho_isEEGap);
   fChain->SetBranchAddress("pho_isEBEEGap", pho_isEBEEGap, &b_pho_isEBEEGap);
   fChain->SetBranchAddress("pho_zernike20", pho_zernike20, &b_pho_zernike20);
   fChain->SetBranchAddress("pho_zernike42", pho_zernike42, &b_pho_zernike42);
   fChain->SetBranchAddress("pho_e2nd", pho_e2nd, &b_pho_e2nd);
   fChain->SetBranchAddress("pho_e2x5max", pho_e2x5max, &b_pho_e2x5max);
   fChain->SetBranchAddress("pho_e2x5right", pho_e2x5right, &b_pho_e2x5right);
   fChain->SetBranchAddress("pho_e2x5left", pho_e2x5left, &b_pho_e2x5left);
   fChain->SetBranchAddress("pho_e2x5top", pho_e2x5top, &b_pho_e2x5top);
   fChain->SetBranchAddress("pho_e2x5bottom", pho_e2x5bottom, &b_pho_e2x5bottom);
   fChain->SetBranchAddress("pho_eright", pho_eright, &b_pho_eright);
   fChain->SetBranchAddress("pho_eleft", pho_eleft, &b_pho_eleft);
   fChain->SetBranchAddress("pho_etop", pho_etop, &b_pho_etop);
   fChain->SetBranchAddress("pho_ebottom", pho_ebottom, &b_pho_ebottom);
   fChain->SetBranchAddress("pho_r19", pho_r19, &b_pho_r19);
   fChain->SetBranchAddress("pho_maxoraw", pho_maxoraw, &b_pho_maxoraw);
   fChain->SetBranchAddress("pho_cep", pho_cep, &b_pho_cep);
   fChain->SetBranchAddress("pho_lambdaratio", pho_lambdaratio, &b_pho_lambdaratio);
   fChain->SetBranchAddress("pho_lambdadivcov", pho_lambdadivcov, &b_pho_lambdadivcov);
   fChain->SetBranchAddress("pho_etawidth", pho_etawidth, &b_pho_etawidth);
   fChain->SetBranchAddress("pho_brem", pho_brem, &b_pho_brem);
   fChain->SetBranchAddress("pho_smaj", pho_smaj, &b_pho_smaj);
   fChain->SetBranchAddress("pho_e2x2", pho_e2x2, &b_pho_e2x2);
   fChain->SetBranchAddress("pho_seed_time", pho_seed_time, &b_pho_seed_time);
   fChain->SetBranchAddress("pho_seed_outoftimechi2", pho_seed_outoftimechi2, &b_pho_seed_outoftimechi2);
   fChain->SetBranchAddress("pho_seed_chi2", pho_seed_chi2, &b_pho_seed_chi2);
   fChain->SetBranchAddress("pho_seed_recoflag", pho_seed_recoflag, &b_pho_seed_recoflag);
   fChain->SetBranchAddress("pho_seed_severity", pho_seed_severity, &b_pho_seed_severity);
   fChain->SetBranchAddress("pho_ecalsumetconedr04", pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
   fChain->SetBranchAddress("pho_hcalsumetconedr04", pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
   fChain->SetBranchAddress("pho_trksumptsolidconedr04", pho_trksumptsolidconedr04, &b_pho_trksumptsolidconedr04);
   fChain->SetBranchAddress("pho_trksumpthollowconedr04", pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
   fChain->SetBranchAddress("pho_ntrksolidconedr04", pho_ntrksolidconedr04, &b_pho_ntrksolidconedr04);
   fChain->SetBranchAddress("pho_ntrkhollowconedr04", pho_ntrkhollowconedr04, &b_pho_ntrkhollowconedr04);
   fChain->SetBranchAddress("pho_ecalsumetconedr03", pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
   fChain->SetBranchAddress("pho_hcalsumetconedr03", pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
   fChain->SetBranchAddress("pho_trksumptsolidconedr03", pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
   fChain->SetBranchAddress("pho_trksumpthollowconedr03", pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
   fChain->SetBranchAddress("pho_ntrksolidconedr03", pho_ntrksolidconedr03, &b_pho_ntrksolidconedr03);
   fChain->SetBranchAddress("pho_ntrkhollowconedr03", pho_ntrkhollowconedr03, &b_pho_ntrkhollowconedr03);
   fChain->SetBranchAddress("pho_haspixseed", pho_haspixseed, &b_pho_haspixseed);
   fChain->SetBranchAddress("pho_hasconvtks", pho_hasconvtks, &b_pho_hasconvtks);
   fChain->SetBranchAddress("pho_nconv", pho_nconv, &b_pho_nconv);
   fChain->SetBranchAddress("pho_conv_ntracks", pho_conv_ntracks, &b_pho_conv_ntracks);
   fChain->SetBranchAddress("pho_conv_pairinvmass", pho_conv_pairinvmass, &b_pho_conv_pairinvmass);
   fChain->SetBranchAddress("pho_conv_paircotthetasep", pho_conv_paircotthetasep, &b_pho_conv_paircotthetasep);
   fChain->SetBranchAddress("pho_conv_eoverp", pho_conv_eoverp, &b_pho_conv_eoverp);
   fChain->SetBranchAddress("pho_conv_zofprimvtxfromtrks", pho_conv_zofprimvtxfromtrks, &b_pho_conv_zofprimvtxfromtrks);
   fChain->SetBranchAddress("pho_conv_distofminapproach", pho_conv_distofminapproach, &b_pho_conv_distofminapproach);
   fChain->SetBranchAddress("pho_conv_dphitrksatvtx", pho_conv_dphitrksatvtx, &b_pho_conv_dphitrksatvtx);
   fChain->SetBranchAddress("pho_conv_dphitrksatecal", pho_conv_dphitrksatecal, &b_pho_conv_dphitrksatecal);
   fChain->SetBranchAddress("pho_conv_detatrksatecal", pho_conv_detatrksatecal, &b_pho_conv_detatrksatecal);
   fChain->SetBranchAddress("pho_conv_tk1_d0", pho_conv_tk1_d0, &b_pho_conv_tk1_d0);
   fChain->SetBranchAddress("pho_conv_tk1_pout", pho_conv_tk1_pout, &b_pho_conv_tk1_pout);
   fChain->SetBranchAddress("pho_conv_tk1_pin", pho_conv_tk1_pin, &b_pho_conv_tk1_pin);
   fChain->SetBranchAddress("pho_conv_tk2_d0", pho_conv_tk2_d0, &b_pho_conv_tk2_d0);
   fChain->SetBranchAddress("pho_conv_tk2_pout", pho_conv_tk2_pout, &b_pho_conv_tk2_pout);
   fChain->SetBranchAddress("pho_conv_tk2_pin", pho_conv_tk2_pin, &b_pho_conv_tk2_pin);
   fChain->SetBranchAddress("pho_conv_tk1_dz", pho_conv_tk1_dz, &b_pho_conv_tk1_dz);
   fChain->SetBranchAddress("pho_conv_tk1_dzerr", pho_conv_tk1_dzerr, &b_pho_conv_tk1_dzerr);
   fChain->SetBranchAddress("pho_conv_tk1_nh", pho_conv_tk1_nh, &b_pho_conv_tk1_nh);
   fChain->SetBranchAddress("pho_conv_tk2_dz", pho_conv_tk2_dz, &b_pho_conv_tk2_dz);
   fChain->SetBranchAddress("pho_conv_tk2_dzerr", pho_conv_tk2_dzerr, &b_pho_conv_tk2_dzerr);
   fChain->SetBranchAddress("pho_conv_tk2_nh", pho_conv_tk2_nh, &b_pho_conv_tk2_nh);
   fChain->SetBranchAddress("pho_conv_ch1ch2", pho_conv_ch1ch2, &b_pho_conv_ch1ch2);
   fChain->SetBranchAddress("pho_conv_chi2", pho_conv_chi2, &b_pho_conv_chi2);
   fChain->SetBranchAddress("pho_conv_chi2_probability", pho_conv_chi2_probability, &b_pho_conv_chi2_probability);
   fChain->SetBranchAddress("pho_conv_validvtx", pho_conv_validvtx, &b_pho_conv_validvtx);
   fChain->SetBranchAddress("pho_conv_MVALikelihood", pho_conv_MVALikelihood, &b_pho_conv_MVALikelihood);
   fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
   fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
   fChain->SetBranchAddress("pho_conv_vtx", &pho_conv_vtx, &b_pho_conv_vtx);
   fChain->SetBranchAddress("pho_conv_pair_momentum", &pho_conv_pair_momentum, &b_pho_conv_pair_momentum);
   fChain->SetBranchAddress("pho_conv_refitted_momentum", &pho_conv_refitted_momentum, &b_pho_conv_refitted_momentum);
   fChain->SetBranchAddress("pho_conv_vertexcorrected_p4", &pho_conv_vertexcorrected_p4, &b_pho_conv_vertexcorrected_p4);
   fChain->SetBranchAddress("pho_scind", pho_scind, &b_pho_scind);
   fChain->SetBranchAddress("pho_residCorrEnergy", pho_residCorrEnergy, &b_pho_residCorrEnergy);
   fChain->SetBranchAddress("pho_residCorrResn", pho_residCorrResn, &b_pho_residCorrResn);
   fChain->SetBranchAddress("pho_regr_energy", pho_regr_energy, &b_pho_regr_energy);
   fChain->SetBranchAddress("pho_regr_energyerr", pho_regr_energyerr, &b_pho_regr_energyerr);
   fChain->SetBranchAddress("pho_isconv", pho_isconv, &b_pho_isconv);
   fChain->SetBranchAddress("pho_pfiso_myneutral01", pho_pfiso_myneutral01, &b_pho_pfiso_myneutral01);
   fChain->SetBranchAddress("pho_pfiso_myneutral02", pho_pfiso_myneutral02, &b_pho_pfiso_myneutral02);
   fChain->SetBranchAddress("pho_pfiso_myneutral03", pho_pfiso_myneutral03, &b_pho_pfiso_myneutral03);
   fChain->SetBranchAddress("pho_pfiso_myneutral04", pho_pfiso_myneutral04, &b_pho_pfiso_myneutral04);
   fChain->SetBranchAddress("pho_pfiso_myneutral05", pho_pfiso_myneutral05, &b_pho_pfiso_myneutral05);
   fChain->SetBranchAddress("pho_pfiso_myneutral06", pho_pfiso_myneutral06, &b_pho_pfiso_myneutral06);
   fChain->SetBranchAddress("pho_pfiso_myphoton01", pho_pfiso_myphoton01, &b_pho_pfiso_myphoton01);
   fChain->SetBranchAddress("pho_pfiso_myphoton02", pho_pfiso_myphoton02, &b_pho_pfiso_myphoton02);
   fChain->SetBranchAddress("pho_pfiso_myphoton03", pho_pfiso_myphoton03, &b_pho_pfiso_myphoton03);
   fChain->SetBranchAddress("pho_pfiso_myphoton04", pho_pfiso_myphoton04, &b_pho_pfiso_myphoton04);
   fChain->SetBranchAddress("pho_pfiso_myphoton05", pho_pfiso_myphoton05, &b_pho_pfiso_myphoton05);
   fChain->SetBranchAddress("pho_pfiso_myphoton06", pho_pfiso_myphoton06, &b_pho_pfiso_myphoton06);
   fChain->SetBranchAddress("pho_pfiso_mycharged01", &pho_pfiso_mycharged01, &b_pho_pfiso_mycharged01);
   fChain->SetBranchAddress("pho_pfiso_mycharged02", &pho_pfiso_mycharged02, &b_pho_pfiso_mycharged02);
   fChain->SetBranchAddress("pho_pfiso_mycharged03", &pho_pfiso_mycharged03, &b_pho_pfiso_mycharged03);
   fChain->SetBranchAddress("pho_pfiso_mycharged04", &pho_pfiso_mycharged04, &b_pho_pfiso_mycharged04);
   fChain->SetBranchAddress("pho_pfiso_mycharged05", &pho_pfiso_mycharged05, &b_pho_pfiso_mycharged05);
   fChain->SetBranchAddress("pho_pfiso_mycharged06", &pho_pfiso_mycharged06, &b_pho_pfiso_mycharged06);
   fChain->SetBranchAddress("pho_isPFPhoton", pho_isPFPhoton, &b_pho_isPFPhoton);
   fChain->SetBranchAddress("pho_isPFElectron", pho_isPFElectron, &b_pho_isPFElectron);
   fChain->SetBranchAddress("pho_must", pho_must, &b_pho_must);
   fChain->SetBranchAddress("pho_mustnc", pho_mustnc, &b_pho_mustnc);
   fChain->SetBranchAddress("pho_pfpresh1", pho_pfpresh1, &b_pho_pfpresh1);
   fChain->SetBranchAddress("pho_pfpresh2", pho_pfpresh2, &b_pho_pfpresh2);
   fChain->SetBranchAddress("pho_mustenergy", pho_mustenergy, &b_pho_mustenergy);
   fChain->SetBranchAddress("pho_mustenergyout", pho_mustenergyout, &b_pho_mustenergyout);
   fChain->SetBranchAddress("pho_pflowE", pho_pflowE, &b_pho_pflowE);
   fChain->SetBranchAddress("pho_pfdeta", pho_pfdeta, &b_pho_pfdeta);
   fChain->SetBranchAddress("pho_pfdphi", pho_pfdphi, &b_pho_pfdphi);
   fChain->SetBranchAddress("pho_pfclusrms", pho_pfclusrms, &b_pho_pfclusrms);
   fChain->SetBranchAddress("pho_pfclusrmsmust", pho_pfclusrmsmust, &b_pho_pfclusrmsmust);
   fChain->SetBranchAddress("pho_eseffsixix", pho_eseffsixix, &b_pho_eseffsixix);
   fChain->SetBranchAddress("pho_eseffsiyiy", pho_eseffsiyiy, &b_pho_eseffsiyiy);
   fChain->SetBranchAddress("pho_hoe_bc", pho_hoe_bc, &b_pho_hoe_bc);
   fChain->SetBranchAddress("pho_biphi", pho_biphi, &b_pho_biphi);
   fChain->SetBranchAddress("pho_bieta", pho_bieta, &b_pho_bieta);
   fChain->SetBranchAddress("pho_betacry", pho_betacry, &b_pho_betacry);
   fChain->SetBranchAddress("pho_phicry", pho_phicry, &b_pho_phicry);
   fChain->SetBranchAddress("conv_n", &conv_n, &b_conv_n);
   fChain->SetBranchAddress("conv_p4", &conv_p4, &b_conv_p4);
   fChain->SetBranchAddress("conv_ntracks", conv_ntracks, &b_conv_ntracks);
   fChain->SetBranchAddress("conv_pairinvmass", conv_pairinvmass, &b_conv_pairinvmass);
   fChain->SetBranchAddress("conv_paircotthetasep", conv_paircotthetasep, &b_conv_paircotthetasep);
   fChain->SetBranchAddress("conv_eoverp", conv_eoverp, &b_conv_eoverp);
   fChain->SetBranchAddress("conv_distofminapproach", conv_distofminapproach, &b_conv_distofminapproach);
   fChain->SetBranchAddress("conv_dphitrksatvtx", conv_dphitrksatvtx, &b_conv_dphitrksatvtx);
   fChain->SetBranchAddress("conv_dphitrksatecal", conv_dphitrksatecal, &b_conv_dphitrksatecal);
   fChain->SetBranchAddress("conv_detatrksatecal", conv_detatrksatecal, &b_conv_detatrksatecal);
   fChain->SetBranchAddress("conv_dxy", conv_dxy, &b_conv_dxy);
   fChain->SetBranchAddress("conv_dz", conv_dz, &b_conv_dz);
   fChain->SetBranchAddress("conv_lxy", conv_lxy, &b_conv_lxy);
   fChain->SetBranchAddress("conv_lz", conv_lz, &b_conv_lz);
   fChain->SetBranchAddress("conv_zofprimvtxfromtrks", conv_zofprimvtxfromtrks, &b_conv_zofprimvtxfromtrks);
   fChain->SetBranchAddress("conv_nHitsBeforeVtx", &conv_nHitsBeforeVtx, &b_conv_nHitsBeforeVtx);
   fChain->SetBranchAddress("conv_nSharedHits", conv_nSharedHits, &b_conv_nSharedHits);
   fChain->SetBranchAddress("conv_validvtx", conv_validvtx, &b_conv_validvtx);
   fChain->SetBranchAddress("conv_MVALikelihood", conv_MVALikelihood, &b_conv_MVALikelihood);
   fChain->SetBranchAddress("conv_chi2", conv_chi2, &b_conv_chi2);
   fChain->SetBranchAddress("conv_chi2_probability", conv_chi2_probability, &b_conv_chi2_probability);
   fChain->SetBranchAddress("conv_vtx_xErr", conv_vtx_xErr, &b_conv_vtx_xErr);
   fChain->SetBranchAddress("conv_vtx_yErr", conv_vtx_yErr, &b_conv_vtx_yErr);
   fChain->SetBranchAddress("conv_vtx_zErr", conv_vtx_zErr, &b_conv_vtx_zErr);
   fChain->SetBranchAddress("conv_tk1_dz", conv_tk1_dz, &b_conv_tk1_dz);
   fChain->SetBranchAddress("conv_tk2_dz", conv_tk2_dz, &b_conv_tk2_dz);
   fChain->SetBranchAddress("conv_tk1_dzerr", conv_tk1_dzerr, &b_conv_tk1_dzerr);
   fChain->SetBranchAddress("conv_tk2_dzerr", conv_tk2_dzerr, &b_conv_tk2_dzerr);
   fChain->SetBranchAddress("conv_tk1_nh", conv_tk1_nh, &b_conv_tk1_nh);
   fChain->SetBranchAddress("conv_tk2_nh", conv_tk2_nh, &b_conv_tk2_nh);
   fChain->SetBranchAddress("conv_ch1ch2", conv_ch1ch2, &b_conv_ch1ch2);
   fChain->SetBranchAddress("conv_tk1_d0", conv_tk1_d0, &b_conv_tk1_d0);
   fChain->SetBranchAddress("conv_tk1_pout", conv_tk1_pout, &b_conv_tk1_pout);
   fChain->SetBranchAddress("conv_tk1_pin", conv_tk1_pin, &b_conv_tk1_pin);
   fChain->SetBranchAddress("conv_tk2_d0", conv_tk2_d0, &b_conv_tk2_d0);
   fChain->SetBranchAddress("conv_tk2_pout", conv_tk2_pout, &b_conv_tk2_pout);
   fChain->SetBranchAddress("conv_tk2_pin", conv_tk2_pin, &b_conv_tk2_pin);
   fChain->SetBranchAddress("conv_vtx", &conv_vtx, &b_conv_vtx);
   fChain->SetBranchAddress("conv_tk1_pterr", conv_tk1_pterr, &b_conv_tk1_pterr);
   fChain->SetBranchAddress("conv_tk2_pterr", conv_tk2_pterr, &b_conv_tk2_pterr);
   fChain->SetBranchAddress("conv_tk1_etaerr", conv_tk1_etaerr, &b_conv_tk1_etaerr);
   fChain->SetBranchAddress("conv_tk2_etaerr", conv_tk2_etaerr, &b_conv_tk2_etaerr);
   fChain->SetBranchAddress("conv_tk1_thetaerr", conv_tk1_thetaerr, &b_conv_tk1_thetaerr);
   fChain->SetBranchAddress("conv_tk2_thetaerr", conv_tk2_thetaerr, &b_conv_tk2_thetaerr);
   fChain->SetBranchAddress("conv_tk1_phierr", conv_tk1_phierr, &b_conv_tk1_phierr);
   fChain->SetBranchAddress("conv_tk2_phierr", conv_tk2_phierr, &b_conv_tk2_phierr);
   fChain->SetBranchAddress("conv_tk1_lambdaerr", conv_tk1_lambdaerr, &b_conv_tk1_lambdaerr);
   fChain->SetBranchAddress("conv_tk2_lambdaerr", conv_tk2_lambdaerr, &b_conv_tk2_lambdaerr);
   fChain->SetBranchAddress("conv_pair_momentum", &conv_pair_momentum, &b_conv_pair_momentum);
   fChain->SetBranchAddress("conv_refitted_momentum", &conv_refitted_momentum, &b_conv_refitted_momentum);
   fChain->SetBranchAddress("conv_singleleg_momentum", &conv_singleleg_momentum, &b_conv_singleleg_momentum);
   fChain->SetBranchAddress("process_id", &process_id, &b_process_id);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("gv_n", &gv_n, &b_gv_n);
   fChain->SetBranchAddress("gv_pos", &gv_pos, &b_gv_pos);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("pu_zpos", &pu_zpos, &b_pu_zpos);
   fChain->SetBranchAddress("pu_sumpt_lowpt", &pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", &pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_ntrks_lowpt", &pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_ntrks_highpt", &pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("vtx_std_n", &vtx_std_n, &b_vtx_std_n);
   fChain->SetBranchAddress("vtx_std_ntks", vtx_std_ntks, &b_vtx_std_ntks);
   fChain->SetBranchAddress("vtx_std_x2dof", vtx_std_x2dof, &b_vtx_std_x2dof);
   fChain->SetBranchAddress("vtx_std_xyz", &vtx_std_xyz, &b_vtx_std_xyz);
   fChain->SetBranchAddress("vtx_std_dxdydz", &vtx_std_dxdydz, &b_vtx_std_dxdydz);
   fChain->SetBranchAddress("vtx_std_ndof", vtx_std_ndof, &b_vtx_std_ndof);
   fChain->SetBranchAddress("vtx_std_tkind", &vtx_std_tkind, &b_vtx_std_tkind);
   fChain->SetBranchAddress("rho_algo1", &rho_algo1, &b_rho_algo1);
   fChain->SetBranchAddress("rho_algo2", &rho_algo2, &b_rho_algo2);
   fChain->SetBranchAddress("rho_algo3", &rho_algo3, &b_rho_algo3);
   fChain->SetBranchAddress("bs_xyz", &bs_xyz, &b_bs_xyz);
   fChain->SetBranchAddress("bs_sigmaZ", &bs_sigmaZ, &b_bs_sigmaZ);
   fChain->SetBranchAddress("bs_x0Error", &bs_x0Error, &b_bs_x0Error);
   fChain->SetBranchAddress("bs_y0Error", &bs_y0Error, &b_bs_y0Error);
   fChain->SetBranchAddress("bs_z0Error", &bs_z0Error, &b_bs_z0Error);
   fChain->SetBranchAddress("bs_sigmaZ0Error", &bs_sigmaZ0Error, &b_bs_sigmaZ0Error);
   fChain->SetBranchAddress("met_tcmet", &met_tcmet, &b_met_tcmet);
   fChain->SetBranchAddress("met_phi_tcmet", &met_phi_tcmet, &b_met_phi_tcmet);
   fChain->SetBranchAddress("hlt_bit", &hlt_bit, &b_hlt_bit);
   fChain->SetBranchAddress("hlt_n", &hlt_n, &b_hlt_n);
   fChain->SetBranchAddress("hlt_candpath", &hlt_candpath, &b_hlt_candpath);
   fChain->SetBranchAddress("hlt_path_names_HLT", &hlt_path_names_HLT, &b_hlt_path_names_HLT);
   fChain->SetBranchAddress("hlt_p4", &hlt_p4, &b_hlt_p4);
   fChain->SetBranchAddress("jet_algoPF1_n", &jet_algoPF1_n, &b_jet_algoPF1_n);
   fChain->SetBranchAddress("jet_algoPF1_erescale", jet_algoPF1_erescale, &b_jet_algoPF1_erescale);
   fChain->SetBranchAddress("jet_algoPF1_p4", &jet_algoPF1_p4, &b_jet_algoPF1_p4);
   fChain->SetBranchAddress("jet_algoPF1_beta", jet_algoPF1_beta, &b_jet_algoPF1_beta);
   fChain->SetBranchAddress("jet_algoPF1_betaStar", jet_algoPF1_betaStar, &b_jet_algoPF1_betaStar);
   fChain->SetBranchAddress("jet_algoPF1_betaStarClassic", jet_algoPF1_betaStarClassic, &b_jet_algoPF1_betaStarClassic);
   fChain->SetBranchAddress("jet_algoPF1_dR2Mean", jet_algoPF1_dR2Mean, &b_jet_algoPF1_dR2Mean);
   fChain->SetBranchAddress("jet_algoPF1_dRMean", jet_algoPF1_dRMean, &b_jet_algoPF1_dRMean);
   fChain->SetBranchAddress("jet_algoPF1_dZ", jet_algoPF1_dZ, &b_jet_algoPF1_dZ);
   fChain->SetBranchAddress("jet_algoPF1_frac01", jet_algoPF1_frac01, &b_jet_algoPF1_frac01);
   fChain->SetBranchAddress("jet_algoPF1_frac02", jet_algoPF1_frac02, &b_jet_algoPF1_frac02);
   fChain->SetBranchAddress("jet_algoPF1_frac03", jet_algoPF1_frac03, &b_jet_algoPF1_frac03);
   fChain->SetBranchAddress("jet_algoPF1_frac04", jet_algoPF1_frac04, &b_jet_algoPF1_frac04);
   fChain->SetBranchAddress("jet_algoPF1_frac05", jet_algoPF1_frac05, &b_jet_algoPF1_frac05);
   fChain->SetBranchAddress("jet_algoPF1_full_mva", jet_algoPF1_full_mva, &b_jet_algoPF1_full_mva);
   fChain->SetBranchAddress("jet_algoPF1_simple_mva", jet_algoPF1_simple_mva, &b_jet_algoPF1_simple_mva);
   fChain->SetBranchAddress("jet_algoPF1_nCharged", jet_algoPF1_nCharged, &b_jet_algoPF1_nCharged);
   fChain->SetBranchAddress("jet_algoPF1_nNeutrals", jet_algoPF1_nNeutrals, &b_jet_algoPF1_nNeutrals);
   fChain->SetBranchAddress("jet_algoPF1_full_wp_level", jet_algoPF1_full_wp_level, &b_jet_algoPF1_full_wp_level);
   fChain->SetBranchAddress("jet_algoPF1_simple_wp_level", jet_algoPF1_simple_wp_level, &b_jet_algoPF1_simple_wp_level);
   fChain->SetBranchAddress("jet_algoPF1_cutbased_wp_level", jet_algoPF1_cutbased_wp_level, &b_jet_algoPF1_cutbased_wp_level);
   fChain->SetBranchAddress("jet_algoPF1_pfloose", jet_algoPF1_pfloose, &b_jet_algoPF1_pfloose);
   fChain->SetBranchAddress("jet_algoPF1_area", jet_algoPF1_area, &b_jet_algoPF1_area);
   fChain->SetBranchAddress("jet_algoPF1_nvtx", &jet_algoPF1_nvtx, &b_jet_algoPF1_nvtx);
   fChain->SetBranchAddress("jet_algoPF1_beta_ext", &jet_algoPF1_beta_ext, &b_jet_algoPF1_beta_ext);
   fChain->SetBranchAddress("jet_algoPF1_betaStar_ext", &jet_algoPF1_betaStar_ext, &b_jet_algoPF1_betaStar_ext);
   fChain->SetBranchAddress("jet_algoPF1_betaStarClassic_ext", &jet_algoPF1_betaStarClassic_ext, &b_jet_algoPF1_betaStarClassic_ext);
   fChain->SetBranchAddress("jet_algoPF1_full_wp_level_ext", &jet_algoPF1_full_wp_level_ext, &b_jet_algoPF1_full_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF1_simple_wp_level_ext", &jet_algoPF1_simple_wp_level_ext, &b_jet_algoPF1_simple_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF1_cutbased_wp_level_ext", &jet_algoPF1_cutbased_wp_level_ext, &b_jet_algoPF1_cutbased_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF1_full_mva_ext", &jet_algoPF1_full_mva_ext, &b_jet_algoPF1_full_mva_ext);
   fChain->SetBranchAddress("jet_algoPF1_simple_mva_ext", &jet_algoPF1_simple_mva_ext, &b_jet_algoPF1_simple_mva_ext);
   fChain->SetBranchAddress("jet_algoPF3_n", &jet_algoPF3_n, &b_jet_algoPF3_n);
   fChain->SetBranchAddress("jet_algoPF3_erescale", jet_algoPF3_erescale, &b_jet_algoPF3_erescale);
   fChain->SetBranchAddress("jet_algoPF3_p4", &jet_algoPF3_p4, &b_jet_algoPF3_p4);
   fChain->SetBranchAddress("jet_algoPF3_beta", jet_algoPF3_beta, &b_jet_algoPF3_beta);
   fChain->SetBranchAddress("jet_algoPF3_betaStar", jet_algoPF3_betaStar, &b_jet_algoPF3_betaStar);
   fChain->SetBranchAddress("jet_algoPF3_betaStarClassic", jet_algoPF3_betaStarClassic, &b_jet_algoPF3_betaStarClassic);
   fChain->SetBranchAddress("jet_algoPF3_dR2Mean", jet_algoPF3_dR2Mean, &b_jet_algoPF3_dR2Mean);
   fChain->SetBranchAddress("jet_algoPF3_dRMean", jet_algoPF3_dRMean, &b_jet_algoPF3_dRMean);
   fChain->SetBranchAddress("jet_algoPF3_dZ", jet_algoPF3_dZ, &b_jet_algoPF3_dZ);
   fChain->SetBranchAddress("jet_algoPF3_frac01", jet_algoPF3_frac01, &b_jet_algoPF3_frac01);
   fChain->SetBranchAddress("jet_algoPF3_frac02", jet_algoPF3_frac02, &b_jet_algoPF3_frac02);
   fChain->SetBranchAddress("jet_algoPF3_frac03", jet_algoPF3_frac03, &b_jet_algoPF3_frac03);
   fChain->SetBranchAddress("jet_algoPF3_frac04", jet_algoPF3_frac04, &b_jet_algoPF3_frac04);
   fChain->SetBranchAddress("jet_algoPF3_frac05", jet_algoPF3_frac05, &b_jet_algoPF3_frac05);
   fChain->SetBranchAddress("jet_algoPF3_full_mva", jet_algoPF3_full_mva, &b_jet_algoPF3_full_mva);
   fChain->SetBranchAddress("jet_algoPF3_simple_mva", jet_algoPF3_simple_mva, &b_jet_algoPF3_simple_mva);
   fChain->SetBranchAddress("jet_algoPF3_nCharged", jet_algoPF3_nCharged, &b_jet_algoPF3_nCharged);
   fChain->SetBranchAddress("jet_algoPF3_nNeutrals", jet_algoPF3_nNeutrals, &b_jet_algoPF3_nNeutrals);
   fChain->SetBranchAddress("jet_algoPF3_full_wp_level", jet_algoPF3_full_wp_level, &b_jet_algoPF3_full_wp_level);
   fChain->SetBranchAddress("jet_algoPF3_simple_wp_level", jet_algoPF3_simple_wp_level, &b_jet_algoPF3_simple_wp_level);
   fChain->SetBranchAddress("jet_algoPF3_cutbased_wp_level", jet_algoPF3_cutbased_wp_level, &b_jet_algoPF3_cutbased_wp_level);
   fChain->SetBranchAddress("jet_algoPF3_pfloose", jet_algoPF3_pfloose, &b_jet_algoPF3_pfloose);
   fChain->SetBranchAddress("jet_algoPF3_area", jet_algoPF3_area, &b_jet_algoPF3_area);
   fChain->SetBranchAddress("jet_algoPF3_nvtx", &jet_algoPF3_nvtx, &b_jet_algoPF3_nvtx);
   fChain->SetBranchAddress("jet_algoPF3_beta_ext", &jet_algoPF3_beta_ext, &b_jet_algoPF3_beta_ext);
   fChain->SetBranchAddress("jet_algoPF3_betaStar_ext", &jet_algoPF3_betaStar_ext, &b_jet_algoPF3_betaStar_ext);
   fChain->SetBranchAddress("jet_algoPF3_betaStarClassic_ext", &jet_algoPF3_betaStarClassic_ext, &b_jet_algoPF3_betaStarClassic_ext);
   fChain->SetBranchAddress("jet_algoPF3_full_wp_level_ext", &jet_algoPF3_full_wp_level_ext, &b_jet_algoPF3_full_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF3_simple_wp_level_ext", &jet_algoPF3_simple_wp_level_ext, &b_jet_algoPF3_simple_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF3_cutbased_wp_level_ext", &jet_algoPF3_cutbased_wp_level_ext, &b_jet_algoPF3_cutbased_wp_level_ext);
   fChain->SetBranchAddress("jet_algoPF3_full_mva_ext", &jet_algoPF3_full_mva_ext, &b_jet_algoPF3_full_mva_ext);
   fChain->SetBranchAddress("jet_algoPF3_simple_mva_ext", &jet_algoPF3_simple_mva_ext, &b_jet_algoPF3_simple_mva_ext);
   fChain->SetBranchAddress("pho_pfRawEnergy", pho_pfRawEnergy, &b_pho_pfRawEnergy);
   fChain->SetBranchAddress("pho_pfe2x2", pho_pfe2x2, &b_pho_pfe2x2);
   fChain->SetBranchAddress("pho_pfe3x3", pho_pfe3x3, &b_pho_pfe3x3);
   fChain->SetBranchAddress("pho_pfe5x5", pho_pfe5x5, &b_pho_pfe5x5);
   fChain->SetBranchAddress("pho_pfsieie", pho_pfsieie, &b_pho_pfsieie);
   fChain->SetBranchAddress("pho_pfsieip", pho_pfsieip, &b_pho_pfsieip);
   fChain->SetBranchAddress("pho_pfsipip", pho_pfsipip, &b_pho_pfsipip);
   fChain->SetBranchAddress("pho_pfemaxxtal", pho_pfemaxxtal, &b_pho_pfemaxxtal);
   fChain->SetBranchAddress("pho_pfe2nd", pho_pfe2nd, &b_pho_pfe2nd);
   fChain->SetBranchAddress("vtx_std_mva", &vtx_std_mva, &b_vtx_std_mva);
   fChain->SetBranchAddress("vtx_std_vertexz", &vtx_std_vertexz, &b_vtx_std_vertexz);
   fChain->SetBranchAddress("vtx_std_nconv", &vtx_std_nconv, &b_vtx_std_nconv);
   fChain->SetBranchAddress("vtx_std_nlegs", &vtx_std_nlegs, &b_vtx_std_nlegs);
   fChain->SetBranchAddress("vtx_std_pulltoconv", &vtx_std_pulltoconv, &b_vtx_std_pulltoconv);
   fChain->SetBranchAddress("vtx_std_limpulltoconv", &vtx_std_limpulltoconv, &b_vtx_std_limpulltoconv);
   fChain->SetBranchAddress("vtx_std_diphopt", &vtx_std_diphopt, &b_vtx_std_diphopt);
   fChain->SetBranchAddress("vtx_std_diphopx", &vtx_std_diphopx, &b_vtx_std_diphopx);
   fChain->SetBranchAddress("vtx_std_diphopy", &vtx_std_diphopy, &b_vtx_std_diphopy);
   fChain->SetBranchAddress("vtx_std_nch", &vtx_std_nch, &b_vtx_std_nch);
   fChain->SetBranchAddress("vtx_std_ptmax", &vtx_std_ptmax, &b_vtx_std_ptmax);
   fChain->SetBranchAddress("vtx_std_sumpt", &vtx_std_sumpt, &b_vtx_std_sumpt);
   fChain->SetBranchAddress("vtx_std_ptvtx", &vtx_std_ptvtx, &b_vtx_std_ptvtx);
   fChain->SetBranchAddress("vtx_std_pxvtx", &vtx_std_pxvtx, &b_vtx_std_pxvtx);
   fChain->SetBranchAddress("vtx_std_pyvtx", &vtx_std_pyvtx, &b_vtx_std_pyvtx);
   fChain->SetBranchAddress("vtx_std_acosA", &vtx_std_acosA, &b_vtx_std_acosA);
   fChain->SetBranchAddress("vtx_std_ptasym", &vtx_std_ptasym, &b_vtx_std_ptasym);
   fChain->SetBranchAddress("vtx_std_ptbal", &vtx_std_ptbal, &b_vtx_std_ptbal);
   fChain->SetBranchAddress("vtx_std_nchthr", &vtx_std_nchthr, &b_vtx_std_nchthr);
   fChain->SetBranchAddress("vtx_std_ptmax3", &vtx_std_ptmax3, &b_vtx_std_ptmax3);
   fChain->SetBranchAddress("vtx_std_thrust", &vtx_std_thrust, &b_vtx_std_thrust);
   fChain->SetBranchAddress("vtx_std_sumweight", &vtx_std_sumweight, &b_vtx_std_sumweight);
   fChain->SetBranchAddress("vtx_std_sumpt2", &vtx_std_sumpt2, &b_vtx_std_sumpt2);
   fChain->SetBranchAddress("vtx_std_ptratio", &vtx_std_ptratio, &b_vtx_std_ptratio);
   fChain->SetBranchAddress("vtx_std_pzasym", &vtx_std_pzasym, &b_vtx_std_pzasym);
   fChain->SetBranchAddress("vtx_std_spher", &vtx_std_spher, &b_vtx_std_spher);
   fChain->SetBranchAddress("vtx_std_aplan", &vtx_std_aplan, &b_vtx_std_aplan);
   fChain->SetBranchAddress("vtx_std_sumpr", &vtx_std_sumpr, &b_vtx_std_sumpr);
   fChain->SetBranchAddress("vtx_std_sumawy", &vtx_std_sumawy, &b_vtx_std_sumawy);
   fChain->SetBranchAddress("vtx_std_sumtrv", &vtx_std_sumtrv, &b_vtx_std_sumtrv);
   fChain->SetBranchAddress("vtx_std_sumtwd", &vtx_std_sumtwd, &b_vtx_std_sumtwd);
   fChain->SetBranchAddress("vtx_std_awytwdasym", &vtx_std_awytwdasym, &b_vtx_std_awytwdasym);
   fChain->SetBranchAddress("vtx_std_pho1", &vtx_std_pho1, &b_vtx_std_pho1);
   fChain->SetBranchAddress("vtx_std_pho2", &vtx_std_pho2, &b_vtx_std_pho2);
   fChain->SetBranchAddress("pho_matchingConv", &pho_matchingConv, &b_pho_matchingConv);
   fChain->SetBranchAddress("vtx_std_evt_mva", &vtx_std_evt_mva, &b_vtx_std_evt_mva);
   fChain->SetBranchAddress("vtx_std_ranked_list", &vtx_std_ranked_list, &b_vtx_std_ranked_list);
   fChain->SetBranchAddress("vtx_std_sel", &vtx_std_sel, &b_vtx_std_sel);
   fChain->SetBranchAddress("pho_tkiso_recvtx_030_002_0000_10_01", &pho_tkiso_recvtx_030_002_0000_10_01, &b_pho_tkiso_recvtx_030_002_0000_10_01);
   fChain->SetBranchAddress("pho_tkiso_badvtx_040_002_0000_10_01", pho_tkiso_badvtx_040_002_0000_10_01, &b_pho_tkiso_badvtx_040_002_0000_10_01);
   fChain->SetBranchAddress("pho_tkiso_badvtx_id", pho_tkiso_badvtx_id, &b_pho_tkiso_badvtx_id);
   fChain->SetBranchAddress("pho_pfiso_charged_badvtx_04", pho_pfiso_charged_badvtx_04, &b_pho_pfiso_charged_badvtx_04);
   fChain->SetBranchAddress("pho_pfiso_charged_badvtx_id", pho_pfiso_charged_badvtx_id, &b_pho_pfiso_charged_badvtx_id);
   fChain->SetBranchAddress("pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01", &pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01, &b_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01);
   fChain->SetBranchAddress("pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01", pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01, &b_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01);
   fChain->SetBranchAddress("pho_ZeeVal_tkiso_badvtx_id", pho_ZeeVal_tkiso_badvtx_id, &b_pho_ZeeVal_tkiso_badvtx_id);
   fChain->SetBranchAddress("pho_mitmva", &pho_mitmva, &b_pho_mitmva);
   fChain->SetBranchAddress("pho_drtotk_25_99", pho_drtotk_25_99, &b_pho_drtotk_25_99);
   fChain->SetBranchAddress("dipho_n", &dipho_n, &b_dipho_n);
   fChain->SetBranchAddress("dipho_leadind", dipho_leadind, &b_dipho_leadind);
   fChain->SetBranchAddress("dipho_subleadind", dipho_subleadind, &b_dipho_subleadind);
   fChain->SetBranchAddress("dipho_vtxind", dipho_vtxind, &b_dipho_vtxind);
   fChain->SetBranchAddress("dipho_sumpt", dipho_sumpt, &b_dipho_sumpt);
   fChain->SetBranchAddress("pho_cic6cutlevel_lead", &pho_cic6cutlevel_lead, &b_pho_cic6cutlevel_lead);
   fChain->SetBranchAddress("pho_cic6passcuts_lead", &pho_cic6passcuts_lead, &b_pho_cic6passcuts_lead);
   fChain->SetBranchAddress("pho_cic6cutlevel_sublead", &pho_cic6cutlevel_sublead, &b_pho_cic6cutlevel_sublead);
   fChain->SetBranchAddress("pho_cic6passcuts_sublead", &pho_cic6passcuts_sublead, &b_pho_cic6passcuts_sublead);
   fChain->SetBranchAddress("pho_cic4cutlevel_lead", &pho_cic4cutlevel_lead, &b_pho_cic4cutlevel_lead);
   fChain->SetBranchAddress("pho_cic4passcuts_lead", &pho_cic4passcuts_lead, &b_pho_cic4passcuts_lead);
   fChain->SetBranchAddress("pho_cic4cutlevel_sublead", &pho_cic4cutlevel_sublead, &b_pho_cic4cutlevel_sublead);
   fChain->SetBranchAddress("pho_cic4passcuts_sublead", &pho_cic4passcuts_sublead, &b_pho_cic4passcuts_sublead);
   fChain->SetBranchAddress("pho_cic4pfcutlevel_lead", &pho_cic4pfcutlevel_lead, &b_pho_cic4pfcutlevel_lead);
   fChain->SetBranchAddress("pho_cic4pfpasscuts_lead", &pho_cic4pfpasscuts_lead, &b_pho_cic4pfpasscuts_lead);
   fChain->SetBranchAddress("pho_cic4pfcutlevel_sublead", &pho_cic4pfcutlevel_sublead, &b_pho_cic4pfcutlevel_sublead);
   fChain->SetBranchAddress("pho_cic4pfpasscuts_sublead", &pho_cic4pfpasscuts_sublead, &b_pho_cic4pfpasscuts_sublead);
   fChain->SetBranchAddress("pho_genmatched", pho_genmatched, &b_pho_genmatched);
   fChain->SetBranchAddress("pho_regr_energy_otf", pho_regr_energy_otf, &b_pho_regr_energy_otf);
   fChain->SetBranchAddress("pho_regr_energyerr_otf", pho_regr_energyerr_otf, &b_pho_regr_energyerr_otf);
   fChain->SetBranchAddress("jet_algoPF1_genMatched", jet_algoPF1_genMatched, &b_jet_algoPF1_genMatched);
   fChain->SetBranchAddress("jet_algoPF1_vbfMatched", jet_algoPF1_vbfMatched, &b_jet_algoPF1_vbfMatched);
   fChain->SetBranchAddress("jet_algoPF1_genPt", jet_algoPF1_genPt, &b_jet_algoPF1_genPt);
   fChain->SetBranchAddress("jet_algoPF1_genDr", jet_algoPF1_genDr, &b_jet_algoPF1_genDr);
   fChain->SetBranchAddress("jet_algoPF3_genMatched", jet_algoPF3_genMatched, &b_jet_algoPF3_genMatched);
   fChain->SetBranchAddress("jet_algoPF3_vbfMatched", jet_algoPF3_vbfMatched, &b_jet_algoPF3_vbfMatched);
   fChain->SetBranchAddress("jet_algoPF3_genPt", jet_algoPF3_genPt, &b_jet_algoPF3_genPt);
   fChain->SetBranchAddress("jet_algoPF3_genDr", jet_algoPF3_genDr, &b_jet_algoPF3_genDr);
   fChain->SetBranchAddress("shiftMET_pt", shiftMET_pt, &b_shiftMET_pt);
   fChain->SetBranchAddress("shiftMET_phi", shiftMET_phi, &b_shiftMET_phi);
   fChain->SetBranchAddress("smearMET_pt", smearMET_pt, &b_smearMET_pt);
   fChain->SetBranchAddress("smearMET_phi", smearMET_phi, &b_smearMET_phi);
   fChain->SetBranchAddress("shiftsmearMET_pt", shiftsmearMET_pt, &b_shiftsmearMET_pt);
   fChain->SetBranchAddress("shiftsmearMET_phi", shiftsmearMET_phi, &b_shiftsmearMET_phi);
   fChain->SetBranchAddress("shiftscaleMET_pt", shiftscaleMET_pt, &b_shiftscaleMET_pt);
   fChain->SetBranchAddress("shiftscaleMET_phi", shiftscaleMET_phi, &b_shiftscaleMET_phi);
   fChain->SetBranchAddress("shiftMET_eta", shiftMET_eta, &b_shiftMET_eta);
   fChain->SetBranchAddress("shiftMET_e", shiftMET_e, &b_shiftMET_e);
   fChain->SetBranchAddress("shiftscaleMET_eta", shiftscaleMET_eta, &b_shiftscaleMET_eta);
   fChain->SetBranchAddress("shiftscaleMET_e", shiftscaleMET_e, &b_shiftscaleMET_e);
   fChain->SetBranchAddress("gh_gen2reco1", &gh_gen2reco1, &b_gh_gen2reco1);
   fChain->SetBranchAddress("gh_gen2reco2", &gh_gen2reco2, &b_gh_gen2reco2);
   fChain->SetBranchAddress("gh_vbfq1_pdgid", &gh_vbfq1_pdgid, &b_gh_vbfq1_pdgid);
   fChain->SetBranchAddress("gh_vbfq2_pdgid", &gh_vbfq2_pdgid, &b_gh_vbfq2_pdgid);
   fChain->SetBranchAddress("gh_vh_pdgid", &gh_vh_pdgid, &b_gh_vh_pdgid);
   fChain->SetBranchAddress("gh_vh1_pdgid", &gh_vh1_pdgid, &b_gh_vh1_pdgid);
   fChain->SetBranchAddress("gh_vh2_pdgid", &gh_vh2_pdgid, &b_gh_vh2_pdgid);
   fChain->SetBranchAddress("gh_higgs_p4", &gh_higgs_p4, &b_gh_higgs_p4);
   fChain->SetBranchAddress("gh_pho1_p4", &gh_pho1_p4, &b_gh_pho1_p4);
   fChain->SetBranchAddress("gh_pho2_p4", &gh_pho2_p4, &b_gh_pho2_p4);
   fChain->SetBranchAddress("gh_vbfq1_p4", &gh_vbfq1_p4, &b_gh_vbfq1_p4);
   fChain->SetBranchAddress("gh_vbfq2_p4", &gh_vbfq2_p4, &b_gh_vbfq2_p4);
   fChain->SetBranchAddress("gh_vh1_p4", &gh_vh1_p4, &b_gh_vh1_p4);
   fChain->SetBranchAddress("gh_vh2_p4", &gh_vh2_p4, &b_gh_vh2_p4);
   fChain->SetBranchAddress("mu_glo_hasgsftrack", mu_glo_hasgsftrack, &b_mu_glo_hasgsftrack);
   Notify();
}

Bool_t GlobeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GlobeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GlobeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GlobeReader_cxx
