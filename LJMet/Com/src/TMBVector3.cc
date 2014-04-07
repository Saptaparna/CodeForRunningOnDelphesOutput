
#include "LJMet/Com/interface/TMBVector3.h"
#include <limits>
#include <cmath>


//ClassImp(TMBVector3)


TMBVector3::~TMBVector3(){
}



bool TMBVector3::is_equal(double x1, double x2)
{
  if ( x1 == 0.0 ) return x2 == 0.0;
  if ( x2 == 0.0 ) return false;
  double eps = std::numeric_limits<double>::epsilon();
  // 1 tick = 0.5, 2tick = 1, ...
  // Choose an intermediate value
  double maxdif = 0.7*eps;
  double num = x1 - x2;
  double den = std::fabs(x1) + std::fabs(x2);
  double rat = std::fabs(num/den);
  return rat < maxdif;
}


/// Equivalence operator. True if all components are
/// equivalent within machine precision.
Bool_t TMBVector3::is_equal (const TMBVector3 &v) const
{
  return is_equal (v.fX, fX) && is_equal (v.fY, fY) && is_equal (v.fZ, fZ);
}

