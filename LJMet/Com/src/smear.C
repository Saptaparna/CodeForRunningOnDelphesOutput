#include "TRandom3.h"
#include "math.h"
#include <cstdlib>

double smearJet(double pT, double eta, double nPU, bool isCMS=true){
  if(nPU==0) return pT;
  //return pT;  
  double c;
  if(fabs(eta)<1.3) c=0.136;
  else if(fabs(eta)<2) c=0.102;
  else if(fabs(eta)<3) c=0.131;
  else std::exit(-1);

  return gRandom->Gaus(pT,1.1*c*nPU);
}

double smearMET(double MET, double nPU, bool isCMS=true){
  if(nPU==0) return MET;
  //return MET;

  TRandom3 _random(0); // random seed

  return MET+_random.Gaus(0,3.5*sqrt(nPU));
}
