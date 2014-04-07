#ifndef BTAGWEIGHT_H
#define BTAGWEIGHT_H
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std; 


class BTagWeight 
{
 public:
  struct JetInfo {
  JetInfo(float mceff,float datasf,float t=0.) : eff(mceff), sf(datasf) , tag(t){}
    float eff;
    float sf;
    int tag;
  };

 BTagWeight(unsigned int nTaggers) : taggers(nTaggers) {}
   
  // virtual bool filter(vector<int> tags);
  template <class Filter> float weight(vector<vector<JetInfo> > jets);
 private:
  unsigned int taggers;

};





class BTag1MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] == 1;
  }
};



class BTag0MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] == 0;
  }
};



class BTag2MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] == 2;
  }
};



class BTagGE1MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] >= 1;
  }
};



class BTagGE2MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] >= 2;
  }
};













/*
// assume only two OP are used Tight in [1], and Custom in [0]
class BTag1TightFilter
{ 
 public:
  static bool filter(std::vector<int> t)
  {
    return t[1] >= 1;
  }
};

class BTag2TightFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[1] >= 2;
  }
};


class BTag1Tight2CustomFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[1] >= 1 && t[0] >= 2;
  }
};


class BTag1Tight1MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[1] >= 1 && t[0] >= 1;
  }
};


class BTag1MediumFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] >= 1;
  }
};

class BTag2CustomFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return  t[0] >= 2;
  }
};
class BTagAntiMax0CustomFilter
{
 public:
  static bool filter(std::vector<int> t)
  {
    return t[0] == 0;
  }
};
*/




template <class Filter> float BTagWeight::weight(vector<vector<JetInfo> >jets)
{
  unsigned int njets=jets.size();
  std::vector<unsigned int> comb(jets.size());
  for(size_t i=0;i < jets.size(); i++) comb[i]=0;
  unsigned int idx=0;
  unsigned int max=taggers+1; //force sorted tagging //1 << taggers;
  float pMC=0;
  float pData=0;
  if(jets.size()==0) return 0.;
  while(comb[jets.size()-1] < max)
    {
      //std::cout << std::endl << "New comb" << std::endl;
      //for(size_t i=0;i < jets.size(); i++) {std::cout << comb[i] << " ";}
      //std::cout << std::endl;

      std::vector<int> tags;
      for(size_t j=0;j<taggers;j++) tags.push_back(0);

      float mc=1.;
      float data=1.;
      for(size_t j=0;j<njets;j++) // loop on jets
	{
	  //std::cout << std::endl << "Jet" << j << endl;;

	  // if none tagged, take the 1-eff SF for the loosest:
	  float tagMc = 1.-jets[j][0].eff;
	  float tagData = 1.-jets[j][0].eff*jets[j][0].sf;

	  //std::cout <<"mc eff, data sf = " << jets[j][0].eff<< ", " << jets[j][0].sf << std::endl;

	  if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger
	    {
	      int k=comb[j]-1;
	      tagMc=jets[j][k].eff;
	      tagData=jets[j][k].eff*jets[j][k].sf;

	      if(comb[j]< taggers) //if at least one tagged take the SF for the tightest tagger
		{
		  int k1=comb[j];
		  tagMc*=1-jets[j][k1].eff/jets[j][k].eff;
		  tagData*=1-jets[j][k1].eff/jets[j][k].eff*jets[j][k1].sf/jets[j][k].sf;

		}
	    }

	  for(size_t k=0;k< taggers; k++ ) // loop on taggers
	    {
	      bool tagged = (comb[j] > k) ; ///((comb[j] >> k) & 0x1) == 1;
	      if(tagged) tags[k]++;
	    }

	  mc*=tagMc;       
	  data*=tagData;       
   
	}
      if(Filter::filter(tags))
	{
	  //std::cout << mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
	  pMC+=mc;
	  pData+=data;
	  //std::cout << "mc, data,ratioThis,ratioTot " <<  mc << " " << data << " " << data/mc << " " << pData/pMC << std::endl;
	}
      while (comb[idx] == max -1  && idx+1 < jets.size()) idx++; // find first jets for which we did not already test all configs 
      // next combination
      comb[idx]++;  // test a new config for that jet
      for(size_t i=0;i<idx;i++) { comb[i]=0;  } // reset the tested configs for all previous jets
      idx=0;
    }
  if(pMC==0) return 0; 
  return pData/pMC;
}


#endif
