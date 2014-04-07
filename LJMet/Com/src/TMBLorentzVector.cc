#include "LJMet/Com/interface/TMBLorentzVector.h"
#include "TLorentzRotation.h"

TMBLorentzVector::~TMBLorentzVector(){
}



Double_t TMBLorentzVector::E() const {
  /// get energy
  return TMath::Sqrt(fM*fM+Mag32());
}



const TMBVector3& TMBLorentzVector::Vect() const {
   /// Get spatial component (const).
   return *this; 
 }



TMBVector3& TMBLorentzVector::Vect() {
  /// Get spatial component (non-const).
  return *this; 
}



Double_t TMBLorentzVector::operator () (int i) const {
// return i'th component of 4 vector
   switch (i)
   {
   case kX:
   case kY:
   case kZ:
      return Vect ()(i);
      case kE:return E ();
      default:Error ("operator()()", "bad index (%d) returning 0", i);
   }
   return 0.;
}

void TMBLorentzVector::Boost(const TVector3 & b) {
// apply Lorentz boost by vector b
  Double_t b2 = b.Mag2 ();
  Double_t gamma = 1.0 / TMath::Sqrt (1.0 - b2);
  Double_t bp = b.Dot (Vect ());
  Double_t gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;
  Vect () += (gamma2 * bp + gamma * E ()) * b;
  // SetE (gamma * (E () + bp));
}

TMBLorentzVector & TMBLorentzVector::operator *= (const TLorentzRotation & m) {
// apply lorentz rotation, see TLorentzRotation::VectorMultiplication()
  return *this = m.VectorMultiplication((TLorentzVector)*this);
}

TMBLorentzVector & TMBLorentzVector::Transform (const TLorentzRotation & m){
// apply lorentz transformation, see TLorentzRotation::VectorMultiplication()
  return *this = m.VectorMultiplication((TLorentzVector)*this);
}

/// Equivalence operator. True if all components are
/// equivalent within machine precision.
Bool_t TMBLorentzVector::is_equal (const TMBLorentzVector &lv) const
{
  return TMBVector3::is_equal (lv.E(), E()) && Vect().is_equal (lv.Vect());
}

//ClassImp(TMBLorentzVector)
