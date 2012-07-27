// $Id: MitPhysicsUtilsLinkDef.h,v 1.1 2012/05/28 09:25:47 amott Exp $

#ifndef MITPHYSICS_UTILS_LINKDEF_H
#define MITPHYSICS_UTILS_LINKDEF_H

#include "include/GBRTree.h"
#include "include/GBRForest.h"
/* #include "MitPhysics/Utils/interface/EGEnergyCorrector.h" */


#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
/* #pragma link C++ namespace mithep; */

#pragma link C++ class GBRTree+; 
#pragma link C++ class GBRForest+; 


/* #pragma link C++ class mithep::EGEnergyCorrector;  */
#endif
