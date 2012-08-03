#ifndef ROOT_TREE_VECTOR_LINKDEF_H 
#define ROOT_TREE_VECTOR_LINKDEF_H 1

#ifdef __CINT__
#include <vector>
#pragma link off all classes;

#pragma link C++ class std::pair<unsigned int, int>+;
#pragma link C++ class std::vector<std::pair<unsigned int, int> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::pair<float,int> >+;
#pragma link C++ class std::vector<std::pair<int,float> >+;
#pragma link C++ class std::vector<short>+;

#pragma link C++ class VecbosPho+;
#pragma link C++ class VecbosEle+;
#pragma link C++ class VecbosSC+;
#pragma link C++ class VecbosBC+;
#pragma link C++ class VecbosConversion+;
#pragma link C++ class VecbosMu+;
#pragma link C++ class VecbosJet+;
#pragma link C++ class VecbosMET+;
#pragma link C++ class VecbosPFMET+;
#pragma link C++ class VecbosCaloMET+;
#pragma link C++ class VecbosTCMET+;
#pragma link C++ class VecbosGen+;
#pragma link C++ class VecbosVtx+;
#pragma link C++ class VecbosTrack+;
#pragma link C++ class VecbosGsfTrack+;
#pragma link C++ class VecbosMuonTrack+;
#pragma link C++ class VecbosBaseObject+;
#pragma link C++ class VecbosPhysicsObject+;

#pragma link C++ class std::vector<VecbosPho>+;
#pragma link C++ class std::vector<VecbosEle>+;
#pragma link C++ class std::vector<VecbosSC>+;
#pragma link C++ class std::vector<VecbosBC>+;
#pragma link C++ class std::vector<VecbosConversion>+;
#pragma link C++ class std::vector<VecbosMu>+;
#pragma link C++ class std::vector<VecbosJet>+;
#pragma link C++ class std::vector<VecbosMET>+;
#pragma link C++ class std::vector<VecbosCaloMET>+;
#pragma link C++ class std::vector<VecbosPFMET>+;
#pragma link C++ class std::vector<VecbosTCMET>+;
#pragma link C++ class std::vector<VecbosGen>+;
#pragma link C++ class std::vector<VecbosVtx>+;
#pragma link C++ class std::vector<VecbosTrack>+;
#pragma link C++ class std::vector<VecbosGsfTrack>+;
#pragma link C++ class std::vector<VecbosMuonTrack>+;
#pragma link C++ class std::vector<VecbosBaseObject>+;
#pragma link C++ class std::vector<VecbosPhysicsObject>+;

#endif

#endif // ROOT_TREE_VECTOR_LINKDEF_H

