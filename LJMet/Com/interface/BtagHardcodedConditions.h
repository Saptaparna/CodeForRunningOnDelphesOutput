#ifndef BtagHardcodedConditions_h
#define BtagHardcodedConditions_h



#include <iostream>
#include <vector>
#include <algorithm>



class BtagHardcodedConditions{

 public:
    
  BtagHardcodedConditions();
  ~BtagHardcodedConditions();

  /**
   *   Returns the discriminant for a particular algo/OP
   */
  float getDiscriminant(const std::string & op);
  /**
   *   Returns the name of the b-tag algo used in PAT
   */
  std::string getAlgoName(const std::string & op);

  /**
   *   Returns the algorithm tag for a particular algo/OP, e.g. CSVM -> CSV
   */
  inline std::string getAlgoTag(const std::string & op){
    return op.substr(0,op.length()-1);
  }
  /**
   *   Returns the letter of the algo/OP, e.g. CSVM -> M
   */
  inline char getOPTag(const std::string & op) {
    return op[op.length()-1];
  }

  double GetBtagEfficiency(double pt, double eta, std::string tagger="CSVM");
  double GetBtagScaleFactor(double pt, double eta, std::string tagger="CSVM", int year = 2012);
  double GetBtagSFUncertUp(double pt, double eta, std::string tagger="CSVM", int year = 2012);
  double GetBtagSFUncertDown(double pt, double eta, std::string tagger="CSVM", int year = 2012);

  double GetMistagRate(double pt, double eta, std::string tagger="CSVM");
  double GetMistagScaleFactor(double pt, double eta, std::string tagger="CSVM", int year = 2012);
  double GetMistagSFUncertUp(double pt, double eta, std::string tagger="CSVM", int year = 2012);
  double GetMistagSFUncertDown(double pt, double eta, std::string tagger="CSVM", int year = 2012);
    
private:
  double GetBtagScaleFactor2011(double pt, double eta, std::string tagger="CSVM");
  double GetBtagScaleFactor2012(double pt, double eta, std::string tagger="CSVM");
  double GetBtagSFUncertainty2011(double pt, double eta, std::string tagger="CSVM");
  double GetBtagSFUncertainty2012(double pt, double eta, std::string tagger="CSVM");
  double GetMistagSF2011(double pt, double eta, std::string tagger,
	std::string meanminmax);
  double GetMistagSF2012(double pt, double eta, std::string tagger,
	std::string meanminmax);
  inline void fillArray(float* a, float* b, int n) {
    for (int i=0;i<n;++i) a[i] = b[i];
  }

  float SFb_TCHPT_error11[14], SFb_CSVL_error11[14], SFb_CSVM_error11[14], SFb_CSVT_error11[14], SFb_JPL_error11[14], SFb_JPM_error11[14], SFb_JPT_error11[14];
  float SFb_TCHPT_error12[16], SFb_CSVL_error12[16], SFb_CSVM_error12[16], SFb_CSVT_error12[16], SFb_JPL_error12[16], SFb_JPM_error12[16], SFb_JPT_error12[16];
  float ptmin, ptmax;
  typedef std::vector< float > FVec;
  typedef std::vector< float >::iterator FVecI;
  FVec ptRange11, ptRange12;
  inline int findBin(float pt, FVec ptRange){
    return (std::upper_bound(ptRange.begin(), ptRange.end(), pt)-ptRange.begin())-1;
  }

};


#endif
