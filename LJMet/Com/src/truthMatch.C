#include "TLorentzVector.h"



bool truthMatch(float reco_pT, float reco_eta, float reco_phi,
		float truth_pT, float truth_eta, float truth_phi,
		float delta_R=0.5, float delta_pT_over_pT=0.2) {

  if( sqrt( (reco_phi-truth_phi)*(reco_phi-truth_phi) + (reco_eta-truth_eta)*(reco_eta-truth_eta) ) < delta_R){
    if ( fabs((reco_pT-truth_pT)/truth_pT) < delta_pT_over_pT ){
      return true;
    }
  }


  return false;
}


bool truthMatch(TLorentzVector* p1, TLorentzVector* p2, float delta_R=0.5){
  if( p1->DeltaR(*p2) < delta_R) return true;
  return false;
}

