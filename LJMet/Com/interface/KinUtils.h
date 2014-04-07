// -*- C++ -*-
//
// Miscellaneous kinematic utilities
//
// Gena Kukartsev, June 2012
//

#ifndef LJMet_Com_interface_KinUtils_h
#define LJMet_Com_interface_KinUtils_h



#include "TLorentzVector.h"



namespace KinUtils{

  inline TLorentzVector GetTLorentzVector( ROOT::Math::XYZTVector & inVec );

}



}



inline
TLorentzVector KinUtils::GetTLorentzVector( ROOT::Math::XYZTVector & inVec ){
  TLorentzVector outVec(inVec.x(), inVec.y(), inVec.z(), inVec.t());
  return outVec;
}



#endif // LJMet_Com_interface_KinUtils_h
