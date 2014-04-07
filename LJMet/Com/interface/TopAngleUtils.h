/* File TopAngleUtils.hpp
 *
 * Created       : Sat Nov  5 19:32:39 CST 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 * Purpose       : 
 * Last modified : 
 * Comments      : 
 */



#ifndef TopAngleUtils_HPP_
#define TopAngleUtils_HPP_

//#include "cafe/Collection.hpp"
//#include "cafe/Event.hpp"
//#include "tmb_tree/TMBLorentzVector.hpp"
#include "LJMet/Com/interface/TMBLorentzVector.h"

namespace top_cafe {
  
    class TopAngleUtils {
	/**
	   Calculates angular variables
	 */
    
    public:
    
	// Constructor, destructor: 
	TopAngleUtils();
	~TopAngleUtils();
    
	inline bool isValid(const TMBLorentzVector & object) {
	    /// checks if a TMBLorentzVector is valid
	    if (object.Px()==0&&object.Py()==0&&object.Pz()==0&&object.M()==0)
		return false;
	    else
		return true;
	} // isValid

	double CosAngle(const TMBLorentzVector &object1, const TMBLorentzVector & object2, const TMBLorentzVector &frame1 = TMBLorentzVector (0.,0.,0.,0.), const TMBLorentzVector &frame2 = TMBLorentzVector (0.,0.,0.,0.));
    
    private:
    
    
    }; // TopAngleUtils
  
} // namespace top_cafe

#endif
