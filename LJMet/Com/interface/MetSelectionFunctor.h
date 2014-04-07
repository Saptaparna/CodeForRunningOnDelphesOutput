#ifndef LJMet_Com_interface_MetSelectionFunctor_h
#define LJMet_Com_interface_MetSelectionFunctor_h


/**
  \class    MetSelectionFunctor MetSelectionFunctor.h "LJMet/Com/interface/MetSelectionFunctor.h"
  \brief    Met selector for pat::Met

  Selector functor for pat::Met inspired by JetIDSelectionFunctor

  Please see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATSelectors
  for a general overview of the selectors. 

  \author Gena Kukartsev
  \version  $Id: MetSelectionFunctor.h,v 1.17 2012/01/18 19:25:49 msegala Exp $
*/


#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TMath.h>
class MetSelectionFunctor : public Selector<pat::MET>  {
  //class MetSelectionFunctor : public Selector<reco::CaloMET>  {  
 public: // interface
  
  enum Version_t { DATA2010, METSTUDY, N_VERSIONS };
  enum Quality_t { OFF, pt50, MINIMAL, LOOSE, N_QUALITY, DiLepton, Wprime_Mu };
  
  
  MetSelectionFunctor( edm::ParameterSet const & parameters ) {

    // legend for identifying output messages
    legend = "[";
    legend.append("MetSelector");
    legend.append("]: ");

    std::string versionStr = parameters.getParameter<std::string>("version");
    std::string qualityStr = parameters.getParameter<std::string>("quality");

    Version_t version = N_VERSIONS;
    Quality_t quality = N_QUALITY;
    
    if ( versionStr == "DATA2010" ) {

      version = DATA2010;

      if      ( qualityStr == "OFF" )   quality = OFF;
      else if ( qualityStr == "pt50" )  quality = pt50;
      else if ( qualityStr == "DiLepton" )  quality = DiLepton;
      else if ( qualityStr == "Wprime_Mu" ) quality = Wprime_Mu; 
      else                              quality = N_QUALITY;
      
    }
    else if ( versionStr == "METSTUDY" ) {

      version = METSTUDY;

      if      ( qualityStr == "LOOSE" )   quality = LOOSE;
      else if ( qualityStr == "MINIMAL" ) quality = MINIMAL;
      else                                quality = N_QUALITY;

    } else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of DATA2010, METSTUDY" << std::endl;
    }

    initialize( version, quality );

  }
  
  MetSelectionFunctor( Version_t version, Quality_t quality ) {
    initialize(version, quality);
  }
  
  void initialize( Version_t version, Quality_t quality )
  {
    version_ = version;
    quality_ = quality;

    push_back("pt"); // register a cut
    
    // all on by default
    set("pt",true); // cut is considered by default


    std::cout << legend << "New selector initialized" << std::endl;
    std::cout << legend << "Version: " << getVersionString(version_) << std::endl;
    std::cout << legend << "Quality: " << getQualityString(quality_) << std::endl;


    if ( version_ == DATA2010){

      if ( quality_ == OFF ) {
	set("pt", false );
      } 

      else if ( quality_ == pt50 ) {
	set("pt", 50.0 );
      } 
      
      else if ( quality_ == DiLepton ) {
	set("pt", 30.0 );
      } 

      else if ( quality_ == Wprime_Mu ) {
	set("pt", 10.0 );
      } 


    }
    
    else if ( version_ == METSTUDY){
      push_back("pt");
      
      // all on by default
      set("pt", double(0.0));

      if ( quality_ == MINIMAL ) {
	set("pt", false);
      } 

      else if ( quality_ == LOOSE ) {
	set("pt", 10.0);
      }

      else if ( quality_ == N_QUALITY ) {
	set("pt", false);
      }
      
    }
    
    retInternal_ = getBitTemplate();
  }



  // 
  // Accessor from PAT MET
  // 
  bool operator()( const pat::MET & met, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( met, ret );
    else if ( version_ == METSTUDY ) return metStudyCuts( met, ret );
    else {
      return false;
    }
  }
  // accessor from PAT met without the ret
  ///// PUT THIS BACK
  using Selector<pat::MET>::operator();


    


  // 
  // Accessor from RECO MET, reco::CaloMET
  // 
  bool operator()( const reco::CaloMET & met, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( met, ret );
    else if ( version_ == METSTUDY ) return metStudyCuts( met, ret );
    else {
      return false;
    }
  }
  // accessor from RECO met type II without the ret
  //using Selector<reco::CaloMET>::operator();


  // 
  // Accessor from RECO MET, reco::MET
  // 
  bool operator()( const reco::MET & met, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( met, ret );
    else if ( version_ == METSTUDY ) return metStudyCuts( met, ret );
    else {
      return false;
    }
  }
  // accessor from RECO Pf met type I without the ret
  //using Selector<reco::MET>::operator();





  std::string getVersionString(Version_t ver){
    if (ver == METSTUDY) return "METSTUDY";
    else if (ver == DATA2010) return "DATA2010";
    else return "N_VERSIONS";
  }


  std::string getQualityString(Quality_t qual){
    if (qual == MINIMAL) return "MINIMAL";
    else if (qual == LOOSE) return "LOOSE";
    else if (qual == OFF) return "OFF";
    else if (qual == pt50) return "pt50";
    else if (qual == DiLepton) return "DiLepton";
    else if (qual == Wprime_Mu) return "Wprime_Mu";
    else return "N_QUALITY";
  }



  // 
  // cuts for 2010 early data
  // 
  bool data2010Cuts( pat::MET const & met, 
		     pat::strbitset & ret){

    ret.set(false);
    
    // cache some variables
    double pt = met.pt();

    //std::cout <<"pt = "<<pt<<std::endl;
    //double phi = met.phi();
    
    if ( ignoreCut("pt") || pt > cut("pt", double()) ){
      passCut( ret, "pt");
    }
 
    // set all ignored cuts as passed
    setIgnored( ret );    
    
    // combine all cut bits in one and return
    return (bool)ret;

  }
    


  // 
  // Cuts for MET studies.
  // Note that ret bitset is returned by reference.
  // It contains individual cut bit values
  bool metStudyCuts( pat::MET const & met, 
		     pat::strbitset & ret){
    
    ret.set(false);
    
    // cache some variables
    double pt = met.pt();
    //double phi = met.phi();
    
    if ( ignoreCut("pt") || pt > cut("pt", double()) ){
      passCut( ret, "pt");
    }
 
    // set all ignored cuts as passed
    setIgnored( ret );    
    
    // combine all cut bits in one and return
    return (bool)ret;
  }







 private: // member variables
  
  Version_t version_;
  Quality_t quality_;

  std::string legend;

};

#endif
