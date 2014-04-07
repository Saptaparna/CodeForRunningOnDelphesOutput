/* File TopUtils.hpp
 *
 * Created       : Tue Oct  4 16:58:38 CDT 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 *
 * Last modified : 
 * Comments      : 
 */



#ifndef TopUtils_HPP_
#define TopUtils_HPP_


//#include "cafe/Event.hpp"


#include "LJMet/Com/interface/TMBLorentzVector.h"

namespace top_cafe {
  
    class TopUtils{
	/**
	   Contains some general-purpose methods.
	 */    
    public:
    
	// Constructor, destructor: 
	TopUtils();
	~TopUtils();

	TMBLorentzVector Sum(const std::vector<TMBLorentzVector> &myobjects);

	// implement the Lorentz Boost function ourselves!
	void boost_new(TVector3 boostvect, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject);
	void boost_new(double bx, double by, double bz, const TMBLorentzVector &myobject, TMBLorentzVector &myBoostedObject);
	double PtRel(const TMBLorentzVector muon, const TMBLorentzVector jet);
    }; // class TopUtils{
  
} // namespace top_cafe {

#endif
