/* File TopTopologicalVariables.hpp
 *
 * Created       : Tue Oct  4 16:58:28 CDT 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 *
 * Purpose      : Class containing generic methods for 
 *                      calculating topological variables 
 *                     (Instantiate using as argument: vector<TMBLorentzVector> 
 *                      of all objects for which you want any variable 
 *                      to be calculated)  
 *
 * Last modified : Amnon Harel, 06 Mar 2006 (const correctness, cache eigenvalues)
 * Comments      : 
 */



#ifndef TopTopologicalVariables_HPP_
#define TopTopologicalVariables_HPP_

#include <vector>
#include "LJMet/Com/interface/TMBLorentzVector.h"

    class TopTopologicalVariables {
	/**
	   Contains generic methods for calculating topological variables
           Constructor takes most containers (vector, list, Collection, etc.)
           of any object derived from TMBLorentzVector
	*/
    public:
    
	// Constructor, destructor:
        template<class Container>
	TopTopologicalVariables(const Container& objects) 
	  : _pv (0)
        {
            copy(objects.begin(), objects.end(), back_inserter(_myobjects));
        }
	~TopTopologicalVariables();
        TopTopologicalVariables& operator= (const TopTopologicalVariables& that);
        TopTopologicalVariables(const TopTopologicalVariables& that);
      
      // Lambdas (eigen vectors of momentum tensor) and their combinations:
	TVectorD GetMomentumTensorEigenvalues() const;
	double Aplanarity() const;
	double Sphericity() const;
	double C() const;
	double D() const;

        double Centrality() const;
	double Pt() const;
	double Ht() const;
	double H() const;
	double TransverseMass() const;
	double M() const;
	double Mt() const;
	double GeometricMeanPt() const;
	double WeightedEtaRMS() const;
        double MinimumPairMass() const;
        double CosThetaStar() const;
        double CosThetaStarJustZ() const;
    
        double SoftestPt() const;
        double MinDR() const;
        double MaxDR() const;
      
        // Note: this is the minimal relative momentum between two objects, 
        //       a clean variation of KtMin. As for Ktmin, the user is
        //       encouraged to divide by an energy and get PtrelMinPrime
        double PtrelMin() const;

        // Note: this is KtMin, **NOT** KtMinPrime
        // User needs to divide by his/her favorite energy variable to get
        // KtMinPrime
        double KtMin() const;

    private:

        typedef std::vector<TMBLorentzVector>::const_iterator Iterator;
	std::vector<TMBLorentzVector> _myobjects;
    
        mutable TVectorD *_pv; // internal representation: can change even for a const object
        void ensurePV() const; // will calculate the eigen values (_pv) if it hasn't been done yet

     //ClassDef(TopTopologicalVariables,1) // top topological variables
    };
  
#endif













