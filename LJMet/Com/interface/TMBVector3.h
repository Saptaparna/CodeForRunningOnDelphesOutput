#ifndef INCLUDE_TMBVECTOR3
#define INCLUDE_TMBVECTOR3

#include "TMath.h"
#include "TError.h"
#include "TVector2.h"
#include "TMatrix.h"
#include "TRotation.h"

/*! \brief 
The Physics Vector package 

The Physics Vector package consists of five classes:                
  - TVector2 
  - TVector3 
  - TRotation 
  - TLorentzVector 
  - TLorentzRotation 
It is a combination of CLHEPs Vector package written by            
Leif Lonnblad, Andreas Nilsson and Evgueni Tcherniaev              
and a ROOT package written by Pasha Murat.                        
for CLHEP see:  http://wwwinfo.cern.ch/asd/lhc++/clhep/           

\ingroup reco

<H2>
TVector3</H2>
<TT>TVector3</TT> is a general three vector class, which can be used for
the description of different vectors in 3D.
<H3>
Declaration / Access to the components</H3>
<TT>TVector3</TT> has been implemented as a vector of three <TT>Double_t</TT>
variables, representing the cartesian coordinates. By default all components
are initialized to zero:

<P><TT>&nbsp; TVector3 v1;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; //
v1 = (0,0,0)</TT>
<BR><TT>&nbsp; TVector3 v2(1);&nbsp;&nbsp;&nbsp;&nbsp; // v2 = (1,0,0)</TT>
<BR><TT>&nbsp; TVector3 v3(1,2,3); // v3 = (1,2,3)</TT>
<BR><TT>&nbsp; TVector3 v4(v2);&nbsp;&nbsp;&nbsp; // v4 = v2</TT>

<P>It is also possible (but not recommended) to initialize a <TT>TVector3</TT>
with a <TT>Double_t</TT> or <TT>Float_t</TT> C array.

<P>You can get the basic components either by name or by index using <TT>operator()</TT>:

<P><TT>&nbsp; xx = v1.X();&nbsp;&nbsp;&nbsp; or&nbsp;&nbsp;&nbsp; xx =
v1(0);</TT>
<BR><TT>&nbsp; yy = v1.Y();&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
yy = v1(1);</TT>
<BR><TT>&nbsp; zz = v1.Z();&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
zz = v1(2);</TT>

<P>The memberfunctions <TT>SetX()</TT>, <TT>SetY()</TT>, <TT>SetZ()</TT>
and<TT> SetXYZ()</TT> allow to set the components:

<P><TT>&nbsp; v1.SetX(1.); v1.SetY(2.); v1.SetZ(3.);</TT>
<BR><TT>&nbsp; v1.SetXYZ(1.,2.,3.);</TT>
<BR>&nbsp;
<H3>
Noncartesian coordinates</H3>
To get information on the <TT>TVector3</TT> in spherical (rho,phi,theta)
or cylindrical (z,r,theta) coordinates, the
<BR>the member functions <TT>Mag3()</TT> (=magnitude=rho in spherical coordinates),
<TT>Mag32()</TT>, <TT>Theta()</TT>, <TT>CosTheta()</TT>, <TT>Phi()</TT>,
<TT>Perp()</TT> (the transverse component=r in cylindrical coordinates),
<TT>Perp2()</TT> can be used:

<P><TT>&nbsp; Double_t m&nbsp; = v.Mag3();&nbsp;&nbsp;&nbsp; // get magnitude
(=rho=Sqrt(x*x+y*y+z*z)))</TT>
<BR><TT>&nbsp; Double_t m2 = v.Mag32();&nbsp;&nbsp; // get magnitude squared</TT>
<BR><TT>&nbsp; Double_t t&nbsp; = v.Theta();&nbsp; // get polar angle</TT>
<BR><TT>&nbsp; Double_t ct = v.CosTheta();// get cos of theta</TT>
<BR><TT>&nbsp; Double_t p&nbsp; = v.Phi();&nbsp;&nbsp;&nbsp; // get azimuth angle</TT>
<BR><TT>&nbsp; Double_t pp = v.Perp();&nbsp;&nbsp; // get transverse component</TT>
<BR><TT>&nbsp; Double_t pp2= v.Perp2();&nbsp; // get transvers component
squared</TT>

<P>It is also possible to get the transverse component with respect to
another vector:

<P><TT>&nbsp; Double_t ppv1 = v.Perp(v1);</TT>
<BR><TT>&nbsp; Double_t pp2v1 = v.Perp2(v1);</TT>

<P>The pseudorapiditiy ( eta=-ln (tan (phi/2)) ) can be get by <TT>Eta()</TT>
or <TT>PseudoRapidity()</TT>:
<BR>&nbsp;
<BR><TT>&nbsp; Double_t eta = v.PseudoRapidity();</TT>

<P>There are set functions to change one of the noncartesian coordinates:

<P><TT>&nbsp; v.SetTheta(.5); // keeping rho and phi</TT>
<BR><TT>&nbsp; v.SetPhi(.8);&nbsp;&nbsp; // keeping rho and theta</TT>
<BR><TT>&nbsp; v.SetMag3(10.);&nbsp; // keeping theta and phi</TT>
<BR><TT>&nbsp; v.SetPerp(3.);&nbsp; // keeping z and phi</TT>
<BR>&nbsp;
<H3>
Arithmetic / Comparison</H3>
The <TT>TVector3</TT> class provides the operators to add, subtract, scale and compare
vectors:

<P><TT>&nbsp; v3&nbsp; = -v1;</TT>
<BR><TT>&nbsp; v1&nbsp; = v2+v3;</TT>
<BR><TT>&nbsp; v1 += v3;</TT>
<BR><TT>&nbsp; v1&nbsp; = v1 - v3</TT>
<BR><TT>&nbsp; v1 -= v3;</TT>
<BR><TT>&nbsp; v1 *= 10;</TT>
<BR><TT>&nbsp; v1&nbsp; = 5*v2;</TT>

<P><TT>&nbsp; if(v1==v2) {...}</TT>
<BR><TT>&nbsp; if(v1!=v2) {...}</TT>
<BR>&nbsp;
<H3>
Related Vectors</H3>
<TT>&nbsp; v2 = v1.Unit();&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; // get unit
vector parallel to v1</TT>
<BR><TT>&nbsp; v2 = v1.Orthogonal(); // get vector orthogonal to v1</TT>
<H3>
Scalar and vector products</H3>
<TT>&nbsp; s = v1.Dot(v2);&nbsp;&nbsp; // scalar product</TT>
<BR><TT>&nbsp; s = v1 * v2;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; // scalar product</TT>
<BR><TT>&nbsp; v = v1.Cross(v2); // vector product</TT>
<H3>
&nbsp;Angle between two vectors</H3>
<TT>&nbsp; Double_t a = v1.Angle(v2);</TT>
<H3>
Rotations</H3>

<H5>
Rotation around axes</H5>
<TT>&nbsp; v.RotateX(.5);</TT>
<BR><TT>&nbsp; v.RotateY(TMath::Pi());</TT>
<BR><TT>&nbsp; v.RotateZ(angle);</TT>
<H5>
Rotation around a vector</H5>
<TT>&nbsp; v1.Rotate(TMath::Pi()/4, v2); // rotation around v2</TT>
<H5>
Rotation by TRotation</H5>
<TT>TVector3</TT> objects can be rotated by objects of the <TT>TRotation</TT>
class using the <TT>Transform()</TT> member functions,
<BR>the <TT>operator *=</TT> or the <TT>operator *</TT> of the TRotation
class:

<P><TT>&nbsp; TRotation m;</TT>
<BR><TT>&nbsp; ...</TT>
<BR><TT>&nbsp; v1.transform(m);</TT>
<BR><TT>&nbsp; v1 = m*v1;</TT>
<BR><TT>&nbsp; v1 *= m; // Attention v1 = m*v1</TT>
<H5>
Transformation from rotated frame</H5>
<TT>&nbsp; TVector3 direction = v.Unit()</TT>
<BR><TT>&nbsp; v1.RotateUz(direction); // direction must be TVector3 of
unit length</TT>

<P>transforms v1 from the rotated frame (z' parallel to direction, x' in
the theta plane and y' in the xy plane as well as perpendicular to the
theta plane) to the (x,y,z) frame.'

 See copyright statements at end of file.

*/

class TMBVector3 : public TObject {
public:

   TMBVector3(Double_t x = 0.0, Double_t y = 0.0, Double_t z = 0.0):
       fX(x), fY(y), fZ(z) {  /// The constructor.
   }

   TMBVector3(const Double_t *a): fX(a[0]), fY(a[1]), fZ(a[2]) 
   { /// Constructor from an array
   }

   TMBVector3(const Float_t *a): fX(a[0]), fY(a[1]), fZ(a[2]) 
   { /// Constructor from an array
   }

   TMBVector3(const TMBVector3 &p): TObject(p),
      fX(p.fX), fY(p.fY), fZ(p.fZ)
   {  /// Copy constructor.
   }

   TMBVector3(const TVector3& v): fX(v.X()), fY(v.Y()), fZ(v.Z())
   {
      /// conversion constructor from TVector3
   }

     virtual ~TMBVector3();
  // Destructor

  Double_t operator () (int i) const {
     /// Get components by index.
     switch(i) {
     case 0: return fX;
     case 1: return fY;
     case 2: return fZ;
     default: Error("operator()(i)", "bad index (%d) returning &fX",i);
     }
     return 0.;
  }
     
  inline Double_t operator [] (int i) const { /// Get components by index.
     return operator()(i); }


  Double_t & operator () (int i) {
     /// Get components by index.

     switch(i) {
     case 0: return fX;
     case 1: return fY;
     case 2: return fZ;
     default: Error("operator()(i)", "bad index (%d) returning &fX",i);
     }
     return fX;
  }
  inline Double_t & operator [] (int i) { /// Set components by index.
     return operator()(i); }

  Double_t x() const { /// get X component 
  return fX; }
  
   Double_t X() const { /// get X component 
      return x(); }
   Double_t Px() const { /// get X component
      return x(); }

   Double_t y() const { /// get Y component
      return fY; }
   Double_t Y() const { /// get Y component
      return y(); }
   Double_t Py() const { /// get X component
      return y(); }

   Double_t z() const { /// get Y component
      return fZ; }
   Double_t Z() const { /// get Z component
      return z(); }
   Double_t Pz() const { /// get X component
      return z(); }

   void SetX(Double_t a) { /// set X component, see TVector3::SetX()
      fX=a; }
   void SetPx(Double_t a) { /// set X component, see TVector3::SetX()
      SetX(a); }
   void SetY(Double_t a) { /// set X component, see TVector3::SetY()
      fY=a; }
   void SetPy(Double_t a) { /// set X component, see TVector3::SetX()
      SetY(a); }
   void SetZ(Double_t a) { /// set X component, see TVector3::SetZ()
      fZ=a; }
   void SetPz(Double_t a) { /// set X component, see TVector3::SetX()
      SetZ(a); }
   inline void SetXYZ(Double_t x, Double_t y, Double_t z) {
      /// set all members using Cartesian coordinates
       fX=x; fY=y; fZ=z; }
    
   void SetPtEtaPhi(Double_t pt, Double_t eta, Double_t phi) {
      /// set all members using |momentum|, eta, and phi
      Double_t apt = TMath::Abs(pt);
      Double_t d=TMath::Tan(2.0*TMath::ATan(TMath::Exp(-eta)));
      if (d==0.) 
         Warning("SetPtEtaPhi", "tan(2*atan(exp(-eta)))==0.! Setting z=0.");
      SetXYZ(apt*TMath::Cos(phi), apt*TMath::Sin(phi), d==0.?0.:apt/d );
   }
   void SetPtThetaPhi(Double_t pt, Double_t theta, Double_t phi) {
      /// set all members using |momentum|, theta, and phi
      Double_t apt = TMath::Abs(pt);
      fX = apt * TMath::Cos(phi);
      fY = apt * TMath::Sin(phi); 
      Double_t tanTheta = TMath::Tan(theta);
      if (tanTheta==0.) 
         Warning("SetPtThetaPhi", "tan(theta)==0.! Setting z=0.");
      fZ = tanTheta ? apt / tanTheta : 0;
   }
   
   inline void GetXYZ(Double_t *carray) const {
      /// Getters into an arry (no checking of array bounds!)
      carray[0]=fX; carray[1]=fY; carray[2]=fZ; }
   inline void GetXYZ(Float_t *carray) const {
      /// Getters into an arry (no checking of array bounds!)
      carray[0]=fX; carray[1]=fY; carray[2]=fZ; }

   Double_t Eta() const { /// get eta (aka pseudo-rapidity)
      return PseudoRapidity(); }
   Double_t Theta() const { /// get theta
      return fX==.0 && fY==.0 && fZ==.0 ? .0 : TMath::ATan2(Perp(),fZ); }
   Double_t CosTheta() const { /// get cos(Theta)
      Double_t ptot = Mag3();
      return ptot == 0.0 ? 1.0 : fZ/ptot;
   }
   Double_t Phi() const { /// get phi from 0..2Pi
      return (fX==.0 && fY==.0 ? .0 : TVector2::Phi_0_2pi(TMath::ATan2(fY,fX))); }

   Double_t Rho() const { /// get magnitude of momentum ("radius" in spherical coords)
      return Mag3(); }

   virtual Double_t Mag32() const {
      /// The magnitude squared (rho^2 in spherical coordinate system).
      return fX*fX + fY*fY + fZ*fZ; }
   virtual Double_t Mag3() const {
      /// The magnitude (rho in spherical coordinate system).
      return TMath::Sqrt(Mag32()); }
   Double_t P()  const { /// get |P|
      return Mag3(); }

   Double_t DeltaPhi(const TMBVector3 &v) const { 
      /// Get difference in phi, between -pi and pi
      return TVector2::Phi_mpi_pi(Phi()-v.Phi()); }

   Double_t DeltaR(const TMBVector3 &v) const {
      /// Get "distance" dR of v to *this, dR=sqrt(dPhi*dPhi+dEta*dEta)
      Double_t deta = Eta()-v.Eta();
      Double_t dphi = TVector2::Phi_mpi_pi(Phi()-v.Phi());
      return TMath::Sqrt( deta*deta+dphi*dphi );
   }
   Double_t DrEtaPhi(const TMBVector3 &lv) const {
      /// Get "distance" dR of lv to *this, dR=sqrt(dPhi*dPhi+dEta*dEta)
      return DeltaR(lv); }

    //   TVector2 EtaPhiVector() { /// return TVector2(eta,phi)
    //      return TVector2 (Eta(),Phi()); }

   Double_t Angle(const TVector3 & q) const { /// Angle wrt. another vector.
      Double_t ptot2 = Mag32()*q.Mag2();
      if(ptot2 <= 0.) return .0;
      Double_t arg = Dot(q)/TMath::Sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
      return TMath::ACos(arg);
   }

   void SetPhi(Double_t ph) { /// set phi keeping mag and theta constant
      Double_t xy   = Perp();
      SetX(xy*TMath::Cos(ph));
      SetY(xy*TMath::Sin(ph));
   }

  inline void SetTheta(Double_t th) { /// Set theta keeping mag and phi constant
     Double_t ma   = Mag3();
     Double_t ph   = Phi();
     SetX(ma*TMath::Sin(th)*TMath::Cos(ph));
     SetY(ma*TMath::Sin(th)*TMath::Sin(ph));
     SetZ(ma*TMath::Cos(th));
  }

  inline void SetMag3(Double_t ma) { /// Set magnitude keeping theta and phi constant
     Double_t factor = Mag3();
     if (factor == 0) {
        Warning("SetMag3","zero vector can't be stretched");
     }else{
        factor = ma/factor;
        SetX(fX*factor);
        SetY(fY*factor);
        SetZ(fZ*factor);
     }
  }
  inline Double_t Perp2() const { 
     /// The transverse component squared (R^2 in cylindrical coordinate system)
     return fX*fX + fY*fY; }
  inline Double_t Pt() const { 
     /// transverse component; projection of 3 vector onto XY plane (R in cylindrical coords)
     return Perp(); }
  inline Double_t Perp() const { 
     /// The transverse component (R in cylindrical coordinate system)
     return TMath::Sqrt(Perp2()); }
     
  inline void SetPerp(Double_t r) { 
     /// Set the transverse component keeping phi and z constant.
     Double_t p = Perp();
     if (p != 0.0) {
        fX *= r/p;
        fY *= r/p;
     }
  }

  inline Double_t Perp2(const TMBVector3 &p) const { 
     /// The transverse component w.r.t. given axis squared.
     Double_t tot = p.Mag32();
     Double_t ss  = Dot(p);
     return tot > 0.0 ? Mag32()-ss*ss/tot : Mag32();
  }

  inline Double_t Pt(const TMBVector3 &p) const {
     /// The transverse component w.r.t. given axis.
     return Perp(p); }
  inline Double_t Perp(const TMBVector3 &p) const {
     /// The transverse component w.r.t. given axis.
     return TMath::Sqrt(Perp2(p)); }

  inline void SetMag3ThetaPhi(Double_t mag, Double_t theta, Double_t phi) {
     /// Set all components as magnitude, theta, and phi
     Double_t amag = TMath::Abs(mag);
     fX = amag * TMath::Sin(theta) * TMath::Cos(phi);
     fY = amag * TMath::Sin(theta) * TMath::Sin(phi);
     fZ = amag * TMath::Cos(theta);
  }

  inline TMBVector3 & operator = (const TMBVector3 &p) {
     /// Assignment.
     fX = p.fX; fY = p.fY; fZ = p.fZ;
     return *this;
  }

  inline Bool_t operator == (const TMBVector3 &v) const {
     /// Equality operator. Equal if all components X,Y,Z are equal.
     return (v.fX==fX && v.fY==fY && v.fZ==fZ); }
  inline Bool_t operator != (const TMBVector3 &v) const {
     /// Inequality operator. Not equal if any of X,Y,Z are not equal.
     return (v.fX!=fX || v.fY!=fY || v.fZ!=fZ); }

  /// Equivalence operator. True if all components are
  /// equivalent within machine precision.
  Bool_t is_equal (const TMBVector3 &v) const;

  inline TMBVector3 & operator += (const TMBVector3 &p) { /// Addition.
     fX += p.fX; fY += p.fY; fZ += p.fZ;
     return *this;
  }

  inline TMBVector3 & operator -= (const TMBVector3 &p) { /// Subtraction.
     fX -= p.fX; fY -= p.fY; fZ -= p.fZ;
     return *this;
  }

  inline TMBVector3 operator - () const { /// Unary minus.
     return TMBVector3(-fX, -fY, -fZ); }

  inline TMBVector3 & operator *= (Double_t a) { 
     /// Scaling with real numbers.
     fX *= a; fY *= a; fZ *= a;
     return *this;
  }

  inline TMBVector3 Unit() const { /// Unit vector parallel to this.
     Double_t  tot = Mag32();
     TMBVector3 p(fX,fY,fZ);
     return tot>.0 ? p *= (1.0/TMath::Sqrt(tot)) : p;
  }

  inline TMBVector3 Orthogonal() const { /// Vector orthogonal to this.
     Double_t x = fX < 0.0 ? -fX : fX;
     Double_t y = fY < 0.0 ? -fY : fY;
     Double_t z = fZ < 0.0 ? -fZ : fZ;
     if (x < y)
        return x<z ? TMBVector3(0,fZ,-fY) : TMBVector3(fY,-fX,0);
     return y<z ? TMBVector3(-fZ,0,fX) : TMBVector3(fY,-fX,0);
  }

  Double_t Dot(const TMBVector3 &p) const { /// Scalar product.
     return fX*p.fX + fY*p.fY + fZ*p.fZ; }

  inline TMBVector3 Cross(const TMBVector3 &p) const { /// Cross product
     return TMBVector3(fY*p.fZ-p.fY*fZ, fZ*p.fX-p.fZ*fX, fX*p.fY-p.fX*fY);
  }

  Double_t PseudoRapidity() const {
     /// Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))
     double cosTheta = CosTheta();
     if (cosTheta*cosTheta < 1) 
        return (-0.5* TMath::Log( (1.0-cosTheta)/(1.0+cosTheta) ));
     Warning("PseudoRapidity","transvers momentum = 0! return +/- 10e10");
     if (fZ > 0) return (10e10);
     else        return (-10e10);
  }
  void RotateX(Double_t angle) {
     /// Rotates around the x-axis.
     Double_t s = TMath::Sin(angle);
     Double_t c = TMath::Cos(angle);
     Double_t yy = fY;
     fY = c*yy - s*fZ;
     fZ = s*yy + c*fZ;
  }

  void RotateY(Double_t angle) {
     /// Rotates around the y-axis.
     Double_t s = TMath::Sin(angle);
     Double_t c = TMath::Cos(angle);
     Double_t zz = fZ;
     fZ = c*zz - s*fX;
     fX = s*zz + c*fX;
  }

  void RotateZ(Double_t angle) {
     /// Rotates around the z-axis.
     Double_t s = TMath::Sin(angle);
     Double_t c = TMath::Cos(angle);
     Double_t xx = fX;
     fX = c*xx - s*fY;
     fY = s*xx + c*fY;
  }

  void RotateUz(const TMBVector3& NewUzVector) {
     /// Rotates reference frame from Uz to newUz (must be unit vector!)
     Double_t u1 = NewUzVector.fX;
     Double_t u2 = NewUzVector.fY;
     Double_t u3 = NewUzVector.fZ;
     Double_t up = u1*u1 + u2*u2;

     if (up>0.) {
        up = TMath::Sqrt(up);
        Double_t px = fX,  py = fY,  pz = fZ;
        fX = (u1*u3*px - u2*py + u1*up*pz)/up;
        fY = (u2*u3*px + u1*py + u2*up*pz)/up;
        fZ = (u3*u3*px -    px + u3*up*pz)/up;
     }
     else if (u3 < 0.) { fX = -fX; fZ = -fZ; }      // phi=0  teta=pi
     else {};
  }

  void Rotate(Double_t angle, const TMBVector3 &axis) {
     /// Rotates around the axis specified by another vector
     TRotation trans;
     trans.Rotate(angle, (TVector3)axis);
     operator*=(trans);
  }

  TMBVector3 & operator *= (const TRotation &m) {
     /// Transform with a rotation matrix
     return *this = m * (*this);
  }
  TMBVector3 & Transform(const TRotation &m) {
     /// Transform with a rotation matrix
     return *this = m * (*this);
  }

   TVector2 XYvector() const {
     /// TVector2 containing x and y components
     return TVector2(fX,fY); }

   operator TVector3() const { 
      /// cast to a TVector3
      return TVector3(X(),Y(),Z()); }

   void Print(Option_t* option="") const {
      /// print components to stdout
      Printf("%s %s (x,y,z)=(%f,%f,%f) (rho,theta,phi)=(%f,%f,%f)",
             GetName(),GetTitle(),X(),Y(),Z(),
             Mag3(),Theta()*TMath::RadToDeg(),Phi()*TMath::RadToDeg());
   }

   // Helper to test if two FP numbers are equivalent within machine precision.
   static bool is_equal (double x1, double x2);

private:
   Double32_t fX; ///< X component (left-right)
   Double32_t fY; ///< Y component (up-down)
   Double32_t fZ; ///< Z component (along the beam)

     //ClassDef(TMBVector3,1) // A 3D physics vector

};

inline
TMBVector3 operator + (const TMBVector3 &a, const TMBVector3 &b) {
   /// Addition of 3-vectors.
  return TVector3(a.X() + b.X(), a.Y() + b.Y(), a.Z() + b.Z()); }

inline
TMBVector3 operator - (const TMBVector3 &a, const TMBVector3 &b) {
   /// Subtraction of 3-vectors.
   return TVector3(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z()); }

inline
Double_t operator * (const TMBVector3 &a, const TMBVector3 &b) {
   /// Dot product of two vectors
   return a.Dot(b); }

inline
TMBVector3 operator * (const TMBVector3 &p, Double_t a) {
   /// Scalar product of 3-vectors.
   return TVector3(a*p.X(), a*p.Y(), a*p.Z()); }

inline
TMBVector3 operator * (Double_t a, const TMBVector3 &p) {
   /// Scaling of 3-vectors with a real number
   return TVector3(a*p.X(), a*p.Y(), a*p.Z()); }

inline
TMBVector3 operator * (const TMatrix &m, const TMBVector3 &v) {
   /// Transform with matrix
  return TVector3( m(0,0)*v.X()+m(0,1)*v.Y()+m(0,2)*v.Z(),
                   m(1,0)*v.X()+m(1,1)*v.Y()+m(1,2)*v.Z(),
                   m(2,0)*v.X()+m(2,1)*v.Y()+m(2,2)*v.Z());
}

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Author: Pasha Murat, Peter Malzacher   12/02/99

#endif // INCLUDE_TMBVECTOR3
