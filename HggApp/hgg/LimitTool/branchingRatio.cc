#ifndef branchingRatio_cc
#define branchingRatio_cc
  //copied from the twiki: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
class branchingRatio{
public:
  branchingRatio();
  std::map<float,float> BranchingRatios;
  std::map<float,float> Proc8TeV;
};
#endif
branchingRatio::branchingRatio(){
  BranchingRatios[90.]  = 1.23E-03;
  BranchingRatios[95.]  = 1.40E-03;
  BranchingRatios[100.] = 1.59E-03;
  BranchingRatios[105.] = 1.78E-03;
  BranchingRatios[110.] = 1.97E-03;
  BranchingRatios[115.] = 2.13E-03;
  BranchingRatios[120.] = 2.25E-03;
  BranchingRatios[123.] = 2.28E-03;
  BranchingRatios[124.] = 2.29E-03;
  BranchingRatios[125.] = 2.29E-03;
  BranchingRatios[129.] = 2.27E-03;
  BranchingRatios[130.] = 2.26E-03;
  BranchingRatios[135.] = 2.13E-03;
  BranchingRatios[140.] = 1.93E-03;
  BranchingRatios[145.] = 1.67E-03;
  BranchingRatios[150.] = 1.36E-03;
  BranchingRatios[155.] = 0.999E-03;

  Proc8TeV[90.]  = 36.772;
  Proc8TeV[95.]  = 33.153;
  Proc8TeV[100.] = 30.089;
  Proc8TeV[105.] = 27.360;
  Proc8TeV[110.] = 25.012;
  Proc8TeV[115.] = 22.937;
  Proc8TeV[120.] = 21.109;
  Proc8TeV[123.] = 20.113;
  Proc8TeV[124.] = 19.796;
  Proc8TeV[125.] = 19.487;
  Proc8TeV[129.] = 18.317;
  Proc8TeV[130.] = 18.040;
  Proc8TeV[135.] = 16.745;
  Proc8TeV[140.] = 15.578;
  Proc8TeV[145.] = 14.525;
  Proc8TeV[150.] = 13.567;
  Proc8TeV[155.] = 12.678;
}
