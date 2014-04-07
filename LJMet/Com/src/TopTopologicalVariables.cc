/* File TopTopologicalVariables.hpp
 *
 * Created       : Tue Oct  4 16:58:28 CDT 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 *
 * Purpose      : Class containing generic methods for 
 *                       calculating topological variables 
 *                       (Instantiate using as argument: vector<TMBLorentzVector> 
 *                        of all objects for which you want any variable 
 *                        to be calculated)  
 *
 * Last modified : Amnon Harel, 06 Mar 2006 (const correctness, cache eigenvalues)
 * Comments      : 
 */

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"

#include "LJMet/Com/interface/TopTopologicalVariables.h"
#include "LJMet/Com/interface/AnglesUtil.h"
#include <iostream>

using namespace std;

  /// Copy constructor
  TopTopologicalVariables::TopTopologicalVariables(const TopTopologicalVariables& that)
    : _myobjects(that._myobjects)
  {
    _pv = 0;
    if (that._pv) _pv = new TVectorD (*(that._pv));
  }

  /// Asignment operator
  TopTopologicalVariables& TopTopologicalVariables::operator= (const TopTopologicalVariables& that)
  {
    if (this != &that)
      {
	_myobjects = that._myobjects;
	if (_pv) delete _pv;
	if (that._pv) _pv = new TVectorD (*(that._pv));
      }
    return *this;
  }

  /// destructor
  TopTopologicalVariables::~TopTopologicalVariables()
  {
    if (_pv) delete _pv;
  }

  /// this method returns the centrality of a group of objects
  double TopTopologicalVariables::Centrality() const 
  {
    return Ht()/H();  
  } // Centrality()

    /// this method returns the aplanarity of a group of objects
    double TopTopologicalVariables::Aplanarity() const 
    {
      ensurePV();
      return 1.5 * (*_pv)[2]; // alternative syntax: ... * _pv->operator[](2)
    } // Aplanarity()
  
  /// this method returns the sphericity of a group of objects
  double TopTopologicalVariables::Sphericity() const
  {
    ensurePV();
    return 1.5 * ( (*_pv)[2] + (*_pv)[1] );
  } // Sphericity()
    
  /// this method returns the C of a group of objects (~momentum elipsiod surface area)
  double TopTopologicalVariables::C() const
  {
    ensurePV();
    return 3 * ( (*_pv)[0] * (*_pv)[1] + (*_pv)[0] * (*_pv)[2] + (*_pv)[1] * (*_pv)[2]);
  } // C()
    
  /// this method returns the D of a group of objects (~momentum elipsiod volume)
  double TopTopologicalVariables::D() const
  {
    ensurePV();
    return 27 * (*_pv)[0] * (*_pv)[1] * (*_pv)[2];
  } // D()
    
    /// this method returns the Momentum Tensor Eigenvalues of a group of objects
    TVectorD TopTopologicalVariables::GetMomentumTensorEigenvalues() const
    {
      ensurePV ();
      return *_pv; // implicit copy and cast
    } // GetMomentumTensorEigenvalues()

  /// internal method to prepare and cache the eigenvectors
  void TopTopologicalVariables::ensurePV () const
  {
    if (_pv) return;

    TMatrixD MomentumTensor(3,3);
    
    Double_t p2_sum=0.0;
    
    // for each _myobjects:
    for ( unsigned int k=0; k<_myobjects.size(); k++ ) {
      
      for ( int i=0; i<3; i++ ) // px, py, pz
	for ( int j=0; j<3; j++ ) // px, py, pz
	  MomentumTensor(i,j) += _myobjects[k][i]*_myobjects[k][j];
      
      // add the 3-momentum squared to the sum
      p2_sum += _myobjects[k].Mag32();
    } // Loop over _myobjects
    
      // Divide the sums with the p2 sum
    if ( p2_sum != 0 )
      for ( int i=0; i<3; i++ ) // px, py, pz
	for ( int j=0; j<3; j++ ) // px, py, pz
	  MomentumTensor(i,j) /= p2_sum;
    
    _pv = new TVectorD (3);
    MomentumTensor.EigenVectors(*_pv);
  }

  /// this method returns the  Pt of a group of objects
  double TopTopologicalVariables::Pt() const {
    // initialize
    TMBLorentzVector Sum(0.0,0.0,0.0,0.0);
    double pt;
    for (unsigned int i=0; i<_myobjects.size(); i++)
      Sum += _myobjects[i];
    pt = Sum.Pt();
    return pt;
  } // Pt

    /// this method returns the Ht: scalar sum of Pt
    double TopTopologicalVariables::Ht() const {
      // initialize
      double ht = 0.0;
      for (unsigned int i=0; i<_myobjects.size(); i++)
	ht += _myobjects[i].Pt();
      return ht;
    } // Ht

  /// this method returns the H: sum of the (scalar) energy, E
  double TopTopologicalVariables::H() const {
    // initialize
    double H = 0.0;
    for (unsigned int i=0; i<_myobjects.size(); i++)
      H+=_myobjects[i].E();
    return H;
  } // H

    /// this method returns the Transverse Mass, Mt(a, b)  
  /// = sqrt( sum{pT^2} - sum{px^2} - sum{px^2});
  double TopTopologicalVariables::TransverseMass() const {
    // Require at least two objects 
    if(_myobjects.size()<2) return -1.0; 
    // initialize
    double Mt = 0.0;
    double sum_pT = 0.0;
    double sum_px = 0.0;
    double sum_py = 0.0;
    // loop over _myobjects
    for (unsigned int i = 0; i < _myobjects.size(); i++) {
      sum_pT += _myobjects[i].Pt();
      sum_px += _myobjects[i].Px();
      sum_py += _myobjects[i].Py();
    } 
    
    Mt = sum_pT*sum_pT - sum_px*sum_px - sum_py*sum_py;
    
    if ( Mt >= 0.0 )
      Mt = sqrt(Mt); 
    else {
      if(fabs(Mt)<1.0e-6) 
	Mt = 0.0;  // Square of Mt is negative, but small
      else {
	std::cout << "In TopTopologicalVariables\n  Error: Square of transverse mass is negative!\n";
	return -1.0;
      } // Square of Mt is negative and large
    } // Square of Mt is negative
    
    return Mt;
    
  } // TransverseMass()
      
    /// this method returns the invariant mass: sqrt(E^2 - ThreeVector{P}^2 )
    double TopTopologicalVariables::M() const {
      // initialize
      TMBLorentzVector Sum(0,0,0,0);
      for (unsigned int i=0; i<_myobjects.size(); i++)
	Sum += _myobjects[i];
    
      if ( Sum.E()*Sum.E() - Sum.Mag32() < 0 ) {
	std::cout << "Error: Square of Invariant_mass is negative!" << std::endl;
	return -1.0;
      }
      return Sum.M();
    } // M()

  double TopTopologicalVariables::Mt() const {      
    double mt = -1.0;

    if(_myobjects.size() == 2 ) {
      double pt1 = _myobjects[0].Pt();
      double pt2 = _myobjects[1].Pt();
      double phi1 = _myobjects[0].Phi();
      double phi2 = _myobjects[1].Phi();

      double delta_phi = kinem::delta_phi(phi1, phi2);

      mt = TMath::Sqrt(2*pt1*pt2*(1-TMath::Cos(delta_phi)));
    }

    return( mt );
  } // Mt()

    /// this method returns the geometric mean of the Pt()s of the objects
    double TopTopologicalVariables::GeometricMeanPt() const {
      if(_myobjects.size() != 2) return -1.0;
      double meansq = 0.0;
      meansq = _myobjects[0].Pt() * _myobjects[1].Pt();
      if (meansq < 0) {
	std::cout << "Error: Square of GeometricMeanPt is negative!" << std::endl;
	return -1.0;
      }
      return sqrt(meansq);
    } // GeometricMeanPt()
  
  /// this method returns the weighted RMS of eta of the objects
  double TopTopologicalVariables::WeightedEtaRMS() const {
    double meaneta = 0.0;
    for (unsigned i=0; i<_myobjects.size(); i++) {
      meaneta += _myobjects[i].Pt() * _myobjects[i].Eta();
    }
    meaneta /= Ht();

    double weightsum = 0.0;
    double weightedEtaRMS = 0.0;
    for (unsigned i=0; i<_myobjects.size(); i++) {
      double sigmabg = 1.10 - 1.99e-3*_myobjects[i].Pt()  + 15.8/_myobjects[i].Pt();
      double sigmatt = 0.68 - 1.05e-3*_myobjects[i].Pt()  + 13.6/_myobjects[i].Pt();
      double weight = (sigmatt - sigmabg) / sigmatt;  // Usually this is negative.
      weightsum += weight;
      double deta = _myobjects[i].Eta()-meaneta;
      weightedEtaRMS += weight * deta * deta;
    }
    weightedEtaRMS /= weightsum;

    return weightedEtaRMS;
  } // WeightedEtaRMS()

    /// return minimum pair invariant mass
    double TopTopologicalVariables::MinimumPairMass() const
    {
      double minMass = 100000;
      Iterator end = _myobjects.end();
      for(Iterator itr1 = _myobjects.begin(); itr1 != end; ++itr1)
        {
	  Iterator itr2 = itr1 + 1;
	  for( ; itr2 != end; ++itr2)
            {
	      double tempmass = (*itr1 + *itr2).M();
	      if(tempmass < minMass)
                {
		  minMass = tempmass;
                }
            }
        }
      return minMass;
    } // MinimumPairMass()

  /// Note: this is KtMin, **NOT** KtMinPrime
    /// User needs to divide by his/her favorite energy variable to get
    /// KtMinPrime
    double TopTopologicalVariables::KtMin() const
    {
      double ptmin = 0;
      double drmin2 = 100000;
      Iterator end = _myobjects.end();
      for(Iterator itr1 = _myobjects.begin(); itr1 != end; ++itr1)
        {
	  Iterator itr2 = itr1 + 1;
	  for( ; itr2 != end; ++itr2)
            {
	      double dp = kinem::delta_phi(itr1->Phi(), itr2->Phi());
	      double de = fabs(itr1->Eta() - itr2->Eta());
	      double dr2 = (dp*dp + de*de);
	      if (dr2 < drmin2)
                {
		  drmin2 = dr2;
		  ptmin = itr1->Pt() < itr2->Pt() ? itr1->Pt() : itr2->Pt();
                }
            }
        }

      return sqrt(drmin2)*ptmin;
    } // KtMin()
  
  // Note: this is the minimal relative momentum between two objects, 
  //       a clean variation of KtMin. As for Ktmin, the user is
  //       encouraged to divide by an energy and get PtrelMinPrime
  double TopTopologicalVariables::PtrelMin() const
  {
    double ptmin = 999;
    for (unsigned int i=0; i<_myobjects.size(); ++i) {

      for (unsigned int j=0; j<_myobjects.size(); ++j) {

	if (i==j) continue;

	double cur = _myobjects[i].Pt(_myobjects[j]);
	if (cur < ptmin) ptmin = cur;
      }
    }
    return ptmin;
  }
    

  /// Calculate cos(theta*), where theta* is the angle between the
  /// 1st object and the beam axis in the objects' rest frame
  double TopTopologicalVariables::CosThetaStar() const
  {
    if (_myobjects.size() <= 0) return -1;
      
    // initialize
    TMBLorentzVector Sum(0,0,0,0);
    for (unsigned int i=0; i<_myobjects.size(); i++)
      Sum += _myobjects[i];
      
    TVector3 sumBoost (Sum.BoostVector());
      
    TMBLorentzVector lv (_myobjects[0]);
    lv.Boost (-1 * sumBoost); // boost to Sum's rest frame
      
    return cos (lv.Theta());
  } // CosThetaStar()

  /// Calculate cos(theta*), assuming no boost in x and y, where theta*
  /// is the angle between the
  /// 1st object and the beam axis in the objects' rest frame
  double TopTopologicalVariables::CosThetaStarJustZ() const
  {
    if (_myobjects.size() <= 0) return -1;
      
    // initialize
    TMBLorentzVector Sum(0,0,0,0);
    for (unsigned int i=0; i<_myobjects.size(); i++)
      Sum += _myobjects[i];
      
    TVector3 sumBoost (Sum.BoostVector());
      
    TMBLorentzVector lv (_myobjects[0]);
    lv.Boost (0, 0, -sumBoost.Z()); // boost to Sum's rest frame, assuming no x&y boost
      
    return cos (lv.Theta());
  } // CosThetaStar()


  /// Returns the lowest pT of all the objects
  double TopTopologicalVariables::SoftestPt() const
  {
    double minPt = 9999;
    for (unsigned int i=0; i<_myobjects.size(); ++i) {
      if (_myobjects[i].Pt() < minPt) minPt = _myobjects[i].Pt();
    }
    return minPt;
  }

  /// Find the lowest delta R between two objects
  double TopTopologicalVariables::MinDR() const
  {
    double minDr = 9999;
    for (unsigned int i=0; i<_myobjects.size()-1; ++i) {
      for (unsigned int j=i+1; j<_myobjects.size(); ++j) {
	double curDr = _myobjects[i].DeltaR (_myobjects[j]);
	if (curDr < minDr) minDr = curDr;
      }
    }
    return minDr;
  }

  /// Find the largest delta R between two objects
  double TopTopologicalVariables::MaxDR() const
  {
    double maxDr = -1;
    for (unsigned int i=0; i<_myobjects.size()-1; ++i) {
      for (unsigned int j=i+1; j<_myobjects.size(); ++j) {
	double curDr = _myobjects[i].DeltaR (_myobjects[j]);
	if (curDr > maxDr) maxDr = curDr;
      }
    }
    return maxDr;
  }




