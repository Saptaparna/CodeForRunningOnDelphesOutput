#ifndef LJMet_Com_interface_MetRECOpfSelectionFunctor_h
#define LJMet_Com_interface_MetRECOpfSelectionFunctor_h


/**
  \class    MetRECOpfSelectionFunctor MetRECOpfSelectionFunctor.h "LJMet/Com/interface/MetRECOpfSelectionFunctor.h"
  \brief    Met selector for pat::Met

  Selector functor for pat::Met inspired by JetIDSelectionFunctor

  Please see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATSelectors
  for a general overview of the selectors. 

  \author Gena Kukartsev
  \version  $Id: MetRECOpfSelectionFunctor.h,v 1.1 2011/01/27 22:41:58 msegala Exp $
*/


#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TMath.h>

class MetRECOpfSelectionFunctor : public Selector<reco::MET>  {  
 public: // interface
  
  enum Version_t { DATA2010, METSTUDY, N_VERSIONS };
  enum Quality_t { MINIMAL, LOOSE, N_QUALITY };
  
  
  MetRECOpfSelectionFunctor( edm::ParameterSet const & parameters ) {
    std::string versionStr = parameters.getParameter<std::string>("version");
    std::string qualityStr = parameters.getParameter<std::string>("quality");
    minimalPt = parameters.getParameter<double>("minPt");
    Quality_t quality = N_QUALITY;
    
    if ( versionStr == "DATA2010" ) {
      if      ( qualityStr == "LOOSE" )   quality = LOOSE;
      else                                quality = MINIMAL;
      
      initialize( DATA2010, quality );
    }
    else if ( versionStr == "METSTUDY" ) {
      initialize( METSTUDY, quality );
    } else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of DATA2010, METSTUDY" << std::endl;
    }
  }
  
  MetRECOpfSelectionFunctor( Version_t version, Quality_t quality ) {
    initialize(version, quality);
  }
  
  void initialize( Version_t version, Quality_t quality )
  {
    version_ = version;
    quality_ = quality;

    push_back("MET pT"); // register a cut
    
    // all on by default
    set("MET pT"); // cut is considered by default
   

    // now set the return values for the ignored parts
    // this is a setup for when the cuts are settled
    if ( quality_ == MINIMAL ) {
      set("MET pT", false );
    } 
    else if ( quality_ == LOOSE ) {
      set("MET pT", 10.0 );
    }
    // this is the METSTUDY case, our work in progress
    else{
      set("MET pT", minimalPt);
    }

    retInternal_ = getBitTemplate();
  }

  // 
  // Accessor from PAT MET
  // 
  bool operator()( const pat::MET & met, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( met, ret );
    if ( version_ == METSTUDY ) return metStudyCuts( met, ret );
    else {
      return false;
    }
  }
  // accessor from PAT met without the ret
  ///// PUT THIS BACK
  //using Selector<pat::MET>::operator();


  // 
  // cuts for 2010 early data
  // 
  bool data2010Cuts( pat::MET const & met, 
		     pat::strbitset & ret){
    return false;
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
    
    if ( ignoreCut("MET pT") || pt > cut("MET pT", double()) ){
      passCut( ret, "MET pT");
    }
 
    // set all ignored cuts as passed
    setIgnored( ret );    
    
    // combine all cut bits in one and return
    return (bool)ret;
  }
  


  // 
  // Accessor from RECO MET, reco::MET
  // 
  bool operator()( const reco::MET & met, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( met, ret );
    if ( version_ == METSTUDY ) return metStudyCuts( met, ret );
    else {
      return false;
    }
  }
  // accessor from RECO pF met type I without the ret
  using Selector<reco::MET>::operator();


  // 
  // cuts for 2010 early data
  // 
  bool data2010Cuts( reco::MET const & met, 
		     pat::strbitset & ret){
    return false;
  }
    
  // 
  // Cuts for MET studies.
  // Note that ret bitset is returned by reference.
  // It contains individual cut bit values
  bool metStudyCuts( reco::MET const & met, 
		     pat::strbitset & ret){
    
    ret.set(false);
    
    // cache some variables
    double pt = met.pt();
    //double phi = met.phi();
    
    if ( ignoreCut("MET pT") || pt > cut("MET pT", double()) ){
      passCut( ret, "MET pT");
    }
 
    // set all ignored cuts as passed
    setIgnored( ret );    
    
    // combine all cut bits in one and return
    return (bool)ret;
  }







 private: // member variables
  
  Version_t version_;
  Quality_t quality_;
  double minimalPt;
};

#endif
