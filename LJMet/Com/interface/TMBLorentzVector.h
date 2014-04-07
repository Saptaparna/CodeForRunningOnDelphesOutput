#ifndef INCLUDE_TMBLORENTZVECTOR
#define INCLUDE_TMBLORENTZVECTOR


// See copyright statements at end of file.

#include "LJMet/Com/interface/TMBVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"

class TLorentzRotation;

/** 
// TMBLorentzVector is based on ROOT's TLorentzVector, implementing all of its
// interfaces. Changes compared to TLorentzVector:
// * derives from TMBVector3 (instead of having a TVector3 member)
// * mass is stored as member, not energy
// * constructor allows several parameter sets, not just xyzE
// * Rapidity() is rapidity, not pseudo-rapidity
// * Some protections against divisions by zero
// * some virtual methods
//
// NOTE on mass: 
//   by setting an energy < |momentum|, Mag2() is negative.
//   In that case, M() will return -sqrt(-Mag2()).
// \ingroup  reco
*/

class TMBLorentzVector: public TMBVector3 {

public:

    /// Safe indexing of the coordinates when using with matrices, arrays, etc.
   enum ECoordinates { 
       kX = 0, 
       kY = 1, 
       kZ = 2, 
       kE = 3, 
       kNUM_COORDINATES = 4, 
       kSIZE = kNUM_COORDINATES 
   };

    /// Initializer sets for the constructor
   enum EInitializer { 
       kXYZE, ///< the four values are given in Cartesian coords (XYZ) and energy
       kXYZM, ///< the four values are given in Cartesian coords (XYZ) and mass
       kPEtaPhiE, ///< the four values are momentum, eta, phi, and energy
       kPEtaPhiM, ///< the four values are momentum, eta, phi, and mass 
       kPtEtaPhiE, ///< the four values are transverse momentum (see Pt()), eta, phi, and energy
       kPtEtaPhiM ///< the four values are transverse momentum (see Pt()), eta, phi, and mass
   };

   TMBLorentzVector(Double_t a = 0.0, Double_t b = 0.0,
		    Double_t c = 0.0, Double_t d = 0.0, 
		    EInitializer init = kXYZE) 
       : TMBVector3 ()
   {
       /// Initialize a TMBLorentzVector.
       /// The meaning of a,b,c,d depends on the value of init,
       /// see EInitializer
       Set(a,b,c,d,init);
   }

   TMBLorentzVector(const TVector3& vect, Double_t v4, EInitializer init=kXYZE):
      TMBVector3(vect) {
       /// Initialize a TMBLorentzVector with a 3 vector and energy or mass,
       /// depending on the value of init (e.g. it's the mass for init==kPtEtaPhiM)
      switch (init) {
         case kXYZM: case kPEtaPhiM: case kPtEtaPhiM: fM=v4;
         default: SetE(v4);
      }
   }
   TMBLorentzVector(const Double_t * a, EInitializer init=kXYZE) 
    : TMBVector3()
   {
       /// Initialize a TMBLorentzVector with an array (array size not checked!)
       /// The meaning of the array values depends on the value of init,
       /// see EInitializer
       Set(a[0],a[1],a[2],a[3],init);
   }
   TMBLorentzVector(const Float_t * a, EInitializer init=kXYZE) {
       /// Initialize a TMBLorentzVector with an array (array size not checked!)
       /// The meaning of the array values depends on the value of init,
       /// see EInitializer
       Set(a[0],a[1],a[2],a[3],init);
   }

   TMBLorentzVector(const TMBLorentzVector& lv): 
     TMBVector3(lv.Vect()), fM(lv.fM) {
       /// Copy constructor.
     }

   TMBLorentzVector(const TLorentzVector& lv): TMBVector3(lv.Vect()) {
       /// Conversion constructor from a TLorentzVector
      SetE(lv.E()); }

     virtual ~TMBLorentzVector();
    /// Destructor


    /// @name Getters
    //@{

   Double_t M() const { 
      /// get mass, see remark on mass in TMBLorentzVector class description
      return fM; }

   virtual Double_t E()  const;

   Double_t Energy() const { /// get energy
      return E(); }
   Double_t T() const { /// get energy / time component
      return E(); }

   virtual const TMBVector3& Vect() const;  // Get spatial component (non-const).
     
   virtual TMBVector3& Vect();// Get spatial component (non-const).

   void GetXYZT(Double_t *carray) const { 
      /// Getters into an arry (no checking of array bounds!)
      GetXYZ(carray); carray[3]=E(); }
   void GetXYZT(Float_t *carray) const {
      /// Getters into an arry (no checking of array bounds!)
      GetXYZ(carray); carray[3]=E(); }

   Double_t operator () (int i) const;
   Double_t operator [] (int i) const { 
      /// Get components by index (as if it was an array).
      /// You can use ECoordinates to index a dimension.
      return operator ()(i); } 

    //@}

    /// @name Setters
    //@{

   void Set(Double_t a = 0.0, Double_t b = 0.0,
       Double_t c = 0.0, Double_t d = 0.0, 
       EInitializer init = kXYZE) {
       /// Initialize a TMBLorentzVector.
       /// The meaning of a,b,c,d depends on the value of init,
       /// see EInitializer
       switch (init) {
       case kXYZM: SetXYZM(a,b,c,d); break;
       case kPEtaPhiE: SetPEtaPhiE(a,b,c,d); break;
       case kPEtaPhiM: SetPEtaPhiM(a,b,c,d); break;
       case kPtEtaPhiE: SetPtEtaPhiE(a,b,c,d); break;
       case kPtEtaPhiM: SetPtEtaPhiM(a,b,c,d); break;
       default: SetXYZT(a,b,c,d);
       }
   }

   void SetM(Double_t a) { /// set mass
      fM=a; }

   void SetT(Double_t a) { /// set E component. See remark on mass in TMBLorentzVector class description
      Double_t m2=a*a-Mag32(); 
      fM=m2>0?TMath::Sqrt(m2) : -TMath::Sqrt(-m2); 
   }
   void SetE(Double_t a) { /// set E component. See remark on mass in TMBLorentzVector class description
      SetT(a); }

   void SetVect(const TVector3 & vect3) { /// Set spatial component.
      Vect()=vect3; }

   void SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e) {
      /// set all members using Cartesian coordinates and energy
      SetXYZT(px,py,pz,e); }
   void SetXYZT(Double_t  x, Double_t  y, Double_t  z, Double_t t) {
      /// set all members using Cartesian coordinates and energy
      SetXYZ(x,y,z); SetE(t); }
   void SetXYZM(Double_t  x, Double_t  y, Double_t  z, Double_t m) {
      /// set all members using Cartesian coordinates and mass
      SetXYZ(x,y,z); SetM(m); }
   void SetPtEtaPhiE(Double_t pt, Double_t eta, Double_t phi, Double_t e) {
      /// set all members using transverse momentum, eta, phi and energy
      SetPtEtaPhi(pt,eta,phi); SetE(e); }
   void SetPtEtaPhiM(Double_t pt, Double_t eta, Double_t phi, Double_t m) {
      /// set all members using transverse momentum, eta, phi and mass
      SetPtEtaPhi(pt,eta,phi); fM=m; }
   void SetPEtaPhiE(Double_t p, Double_t eta, Double_t phi, Double_t e) {
      /// set all members using |momentum|, eta, phi and energy
      SetPtEtaPhi(p*TMath::Sin(2.*TMath::ATan(TMath::Exp(-eta))),eta,phi); SetE(e); }
   void SetPEtaPhiM(Double_t p, Double_t eta, Double_t phi, Double_t m) {
      /// set all members using |momentum|, eta, phi and mass
      SetPtEtaPhi(p*TMath::Sin(2.*TMath::ATan(TMath::Exp(-eta))),eta,phi); fM=m; }

    //@}

   TMBLorentzVector & operator = (const TMBLorentzVector & lv) {
      /// assignment
      Vect()=lv.Vect(); fM = lv.fM; return *this; }
   

   TMBLorentzVector   operator +  (const TMBLorentzVector &lv) const {
      /// 4vector addition
      return TMBLorentzVector(Vect()+lv.Vect(), E()+lv.E()); }
   TMBLorentzVector & operator += (const TMBLorentzVector &lv) {
      /// 4vector addition
      Double_t e=E(); Vect()+=lv.Vect(); SetE(e+lv.E()); return *this; }

   TMBLorentzVector   operator -  (const TMBLorentzVector &lv) const {
      /// 4vector substraction
      return TMBLorentzVector(Vect()-lv.Vect(), E()-lv.E()); }
   TMBLorentzVector & operator -= (const TMBLorentzVector &lv) {
      /// 4vector substraction
      Double_t e=E(); Vect()-=lv.Vect(); SetE(e-lv.E()); return *this; }

   TMBLorentzVector operator - () const {
      /// unary minus
      return TMBLorentzVector(-Vect(),-E()); }

   TMBLorentzVector operator * (Double_t a) const {
      /// scaling with real numbers (const)
      return TMBLorentzVector(Vect()*a, E()*a); }
   TMBLorentzVector & operator *= (Double_t a) {
      /// scaling with real numbers (non-const)
      Vect()*=a; fM*=a; return *this; }

   Bool_t operator == (const TMBLorentzVector &lv) const {
      /// Comparison, equality operator. 
      /// Two lorentz vectors are equal if all their componens are equal.
      return (fM==lv.fM && Vect()==lv.Vect()); }
   Bool_t operator != (const TMBLorentzVector &lv) const {
      /// Comparison, in-equality operator. 
      /// Two lorentz vectors are not equal if any of their componens are not equal.
      return !(operator == (lv)); }

   /// Equivalence operator. True if all components are
   /// equivalent within machine precision.
   Bool_t is_equal (const TMBLorentzVector &lv) const;

   Double_t Mag2() const { /// get invariant mass squared. See remark on mass in TMBLorentzVector class description
      return fM*fM; }
   Double_t M2() const { /// get invariant mass squared. See remark on mass in TMBLorentzVector class description
      return Mag2(); }

   Double_t Mag() const { /// get invariant mass. See remark on mass in TMBLorentzVector class description
      return fM;
   }

   Double_t Mt2() const { /// transverse mass squared
      return E()*E() - Z()*Z(); }
   

   Double_t Mt() const { /// Transverse mass. If Mt2() is negative then -sqrt(-Mt2()) is returned.
      Double_t mm=Mt2(); 
      return mm < 0.0 ? -TMath::Sqrt(-mm) : TMath::Sqrt(mm); }

   Double_t Beta() const { 
      /// get beta
      if (E()!=0) return Mag3() / E(); 
      else Warning("Beta()", "E()==0.!"); return 0.; }
   Double_t Gamma() const {
      /// get gamma
      Double_t b = Beta();
      if (b!=1.) return 1.0/TMath::Sqrt(1- b*b);
      else Warning("Gamma()", "Beta()==1.!"); return 0.; }

   Double_t Dot(const TMBLorentzVector &lv) const {
      /// scalar (dot) product
      return E()*lv.E() - Vect().Dot(lv.Vect()); }
   Double_t operator * (const TMBLorentzVector &lv) const {
      /// scalar (dot) product
      return Dot(lv);}

   void SetVectMag(const TVector3 & spatial, Double_t magnitude) {
      /// set vector and mass
      Vect()=spatial; SetE(TMath::Sqrt(magnitude * magnitude + spatial * spatial)); }
   void SetVectM(const TVector3 & spatial, Double_t mass) {
      /// set vector and mass
      SetVectMag(spatial, mass); }
   void SetVectE(const TVector3 & spatial, Double_t e) {
      /// set vector and energy
      Vect()=spatial; SetE(e); }

   Double_t Plus() const { 
      /// Returns t + z.
      /// Related to the positive/negative light-cone component,
      /// which some define this way and others define as (t +/- z)/sqrt(2)
      /// CAVEAT: The values returned are T{+,-}Z. It is known that some authors
      /// find it easier to define these components as (T{+,-}Z)/sqrt(2). Thus
      /// check what definition is used in the physics you're working in and adapt
      /// your code accordingly.
      return E() + Vect().Z(); }
   Double_t Minus() const { 
      /// Returns t - z.
      /// Related to the positive/negative light-cone component,
      /// which some define this way and others define as (t +/- z)/sqrt(2)
      /// CAVEAT: The values returned are T{+,-}Z. It is known that some authors
      /// find it easier to define these components as (T{+,-}Z)/sqrt(2). Thus
      /// check what definition is used in the physics you're working in and adapt
      /// your code accordingly.
      return E() - Vect().Z(); }

   TVector3 BoostVector() const {
      /// Returns the spatial components divided by the time component.
      if (E()!=0.) return 1./E()*Vect();
      else Warning("BoostVector()","E()==0.!"); return TVector3(); }

   void Boost(Double_t x, Double_t y, Double_t z) { 
      /// Lorentz boost by vector xyz
      Boost(TVector3(x,y,z)); }
   void Boost(const TVector3 & b);

   Double_t Rapidity() const { 
      /// get rapidity (i.e. NOT pseudo-rapidity!)
      if (E()!=Pz()) return .5*log( (E()+Pz()) / (E()-Pz()) );
      else Warning("Rapidity", "E()==Pz()!"); return 0.; }


   TMBLorentzVector & operator *= (const TRotation &m) {
      /// rotate using a TRotation
      Vect() *= m; return *this; }
   TMBLorentzVector & Transform(const TRotation &m) {
      /// transform using a TRotation, see TVector3::Transform()
      Vect().Transform(m); return *this; }
   TMBLorentzVector & operator *= (const TLorentzRotation &);
   TMBLorentzVector & Transform(const TLorentzRotation &);
   /// Transformation with HepLorenzRotation.

   operator TLorentzVector() const { 
      /// cast to a TLorentzVector
      return TLorentzVector(X(),Y(),Z(),E()); }

private:
   Double32_t fM; ///< mass
     //ClassDef(TMBLorentzVector,1) // A four vector with (-,-,-,+) metric
};

inline
TMBLorentzVector operator * (Double_t a, const TMBLorentzVector & p) {
    /// global scope scale-by-float operator
   return TMBLorentzVector(a*p.Vect(), a*p.E());
}

/*************************************************************************
* Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

//------------------------------------------------------------------------------
// Copyright(c) 1995-1997, P.Murat (CDF collaboration, FNAL)
//
// Permission to use, copy, modify and distribute this software and its
// documentation for non-commercial purposes is hereby granted without fee,
// provided that the above copyright notice appears in all copies and
// that both the copyright notice and this permission notice appear in
// the supporting documentation. The authors make no claims about the
// suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.
// *0001 Mar 29 1999 P.Murat: add forgotten scalar product (dot operator)
//------------------------------------------------------------------------------

#endif // INCLUDE_TMBLORENTZVECTOR
