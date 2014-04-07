/* File TopAngleUtils.hpp
 *
 * Created       : Sat Nov  5 19:32:39 CST 2005
 * Author        : Supriya JAIN, sjain@fnal.gov
 *
 * Purpose       : To calculate angular variables
 *
 *
 * Last modified : S. Jain, 3 Nov. 2005 
 * Comments      : Included "using namespace std"
 */

#include "LJMet/Com/interface/TopAngleUtils.h"
#include "LJMet/Com/interface/TopUtils.h"

using namespace std;

namespace top_cafe {
  
    // Constructor, destructor: 
    TopAngleUtils::TopAngleUtils()
    {
    }
    TopAngleUtils::~TopAngleUtils()
    {
    }

    //
    /// calculates the Cosine between two TMBLorentzVectors,
    /// object1, object2, in frame1 and frame2, respectively.
    ///
    /// (All input objects are TMBLorentzVector).
    ///
    /// If no frame is specified, default is "LAB" frame, 
    /// If only frame1 is specified, both objects are boosted in frame1.
    //
    double TopAngleUtils::CosAngle(const TMBLorentzVector &object1, 
				   const TMBLorentzVector & object2, 
				   const TMBLorentzVector &frame1, 
				   const TMBLorentzVector &frame2) 
    {
	//
	// Check input object1 and object2
	if (!isValid(object1) || !isValid(object2))
	    return -2.0;
	//
	// instantiate TopUtils
	TopUtils utils;
	// initialize
	double cos_value;
	TMBLorentzVector new_object1;
	TMBLorentzVector new_object2;
	new_object1 = object1;
	new_object2 = object2;

	// Check if frame1 is specified
	// if, not, then it is the LAB, by default
	if (isValid(frame1)) {
	    // check, if frame is same as object1
	    // if not, boost it in the given frame
	    if (object1 != frame1)
		utils.boost_new(frame1.BoostVector(), object1, new_object1); 
	    // Check if frame2 is also specified
	    if (isValid(frame2)) {
		if (object2 != frame2)
		    utils.boost_new(frame2.BoostVector(), object2, new_object2); 	
	    } // if (isValid(frame2))
	    else {  
		if (object2 != frame1)
		    utils.boost_new(frame1.BoostVector(), object2, new_object2); 
	    } // if, not, then use frame1 for both object1 and object2, by default
	} // if (isValid(frame1))
    
	cos_value = TMath::Cos( new_object1.Angle( new_object2.Vect()) );
	return cos_value;   
    } // CosAngle
  
  
} // using namespace top_cafe










