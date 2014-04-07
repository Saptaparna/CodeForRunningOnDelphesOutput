/* File TopUtils.cpp
 *
 * Created       : Tue Oct  4 16:58:38 CDT 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 * 
 * Purpose       : This contains some general-purpose methods.
 *
 * Last modified : 
 * Comments      : 
 */
#include "LJMet/Com/interface/TopUtils.h"
//#include "cafe/Collection.hpp"
#include <iostream>
#include <algorithm>

using namespace std;
//using namespace cafe;

namespace top_cafe {
  
    // constructor
    TopUtils::TopUtils() {
    }
    // destructor
    TopUtils::~TopUtils() {
    }

    //
    /// this method adds together the elements of vector<TMBLorentzVector>
    //
    TMBLorentzVector TopUtils::Sum(const std::vector<TMBLorentzVector> &myobjects) {
	if(myobjects.size()==0) 	    
	    cout << "In TopUtils\n  Error: Size of vector being summed is zero!\n";
	
	TMBLorentzVector SummedElements;
	SummedElements.SetXYZM(0.0, 0.0, 0.0, 0.0);
	for(unsigned int i = 0; i < myobjects.size(); i++) 
	    SummedElements += myobjects[i];
	return SummedElements;
    } // Sum(const std::vector<TMBLorentzVector> &myobjects)
  
    //
    /// this method does the Lorentz Boost of an object in a new frame
    /// (Note: standard methods in ROOT and cafe/TMBLorentzVetor.hpp
    ///             do not do this quite correctly yet:( )
    //
    void TopUtils::boost_new(TVector3 boostvect, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject)
    {
	boost_new(boostvect.X(),boostvect.Y(),boostvect.Z(),myobject, myBoostedObject);
    } // boost_new(TVector3 boostvect, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject)
  
    void TopUtils::boost_new(double bx, double by, double bz, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject)
    {
	double b2 = bx*bx + by*by + bz*bz;
	register double gamma = 1.0 / sqrt(1.0 - b2);
	register double bp = bx*myobject.Px() + by*myobject.Py() + bz*myobject.Pz();
	register double coef= gamma*((gamma*bp)/(gamma+1) - myobject.T());
    
	myBoostedObject.SetX(myobject.Px() + coef*bx );
	myBoostedObject.SetY(myobject.Py() + coef*by );
	myBoostedObject.SetZ(myobject.Pz() + coef*bz );
	myBoostedObject.SetT(gamma*(myobject.T() - bp));
    } // boost_new(double bx, double by, double bz, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject)

    //
    /// this method calculates the PtRel of a muon: the Pt of a muon w.r.t. the jet axis
    //
    double TopUtils::PtRel(const TMBLorentzVector muon, const TMBLorentzVector jet) {
	// components of sum momentum of muon+jet
	TMBLorentzVector MuonJet;
	MuonJet = muon+jet;

	double muonTimesMuonJet = muon.Px()*MuonJet.Px() + muon.Py()*MuonJet.Py() + muon.Pz()*MuonJet.Pz();
	double pLrel2 = muonTimesMuonJet * muonTimesMuonJet / MuonJet.Mag32();
	double pTrel2 = muon.Mag32() - pLrel2;
	return (pTrel2 > 0)? sqrt(pTrel2): 0;
    }
  
  
} // namespace top_cafe 

