#ifndef PhysicsTools_PatUtils_interface_JetIDSelectionFunctor_h
#define PhysicsTools_PatUtils_interface_JetIDSelectionFunctor_h


/**
  \class    JetIDSelectionFunctor JetIDSelectionFunctor.h "PhysicsTools/Utilities/interface/JetIDSelectionFunctor.h"
  \brief    Jet selector for pat::Jets

  Selector functor for pat::Jets that implements quality cuts based on
  studies of noise patterns. 

  Please see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATSelectors
  for a general overview of the selectors. 

  \author Salvatore Rappoccio
  \author Gena Kukatsev (modified)
  \version  $Id: JetIDSelectionFunctor.h,v 1.12 2012/08/01 22:52:07 kukartse Exp $
*/




#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TMath.h>
class JetIDSelectionFunctor : public Selector<pat::Jet>  {

 public: // interface

  enum Version_t { DATA2010,  N_VERSIONS };
  enum Quality_t { Top_SelV3, N_QUALITY, MEENAKSHI, DiLepton};
  

  JetIDSelectionFunctor( edm::ParameterSet const & parameters ) {

    // legend for identifying output messages
    legend = "[";
    legend.append("JetSelector");
    legend.append("]: ");

    std::string versionStr = parameters.getParameter<std::string>("version");
    std::string qualityStr = parameters.getParameter<std::string>("quality");

    Version_t version = N_VERSIONS;
    Quality_t quality = N_QUALITY;

    if ( versionStr == "DATA2010" ) {

      version = DATA2010;

      if ( qualityStr == "Top_SelV3" ) quality = Top_SelV3;
      else if ( qualityStr == "MEENAKSHI" ) quality = MEENAKSHI;
      else if ( qualityStr == "DiLepton" ) quality = DiLepton;
      else quality = Top_SelV3;
      
    } else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of CRAFT08, N_VERSIONS" << std::endl;
    }

    initialize( version, quality );

  }

  JetIDSelectionFunctor( Version_t version, Quality_t quality ) {
    initialize(version, quality);
  }

 void initialize( Version_t version, Quality_t quality )
  {

    version_ = version;
    quality_ = quality;

    std::cout << legend << "New selector initialized" << std::endl;
    std::cout << legend << "Version: " << getVersionString(version_) << std::endl;
    std::cout << legend << "Quality: " << getQualityString(quality_) << std::endl;

    if ( version_ == DATA2010){
      push_back("pt");
      push_back("eta");
      push_back("CHF" );
      push_back("NHF" );
      push_back("CEF" );
      push_back("NEF" );
      push_back("NCH" );
      push_back("nConstituents");
      push_back("BTAG");

      // all on by default
      set("pt", true);
      set("eta", true);
      set("CHF", 0.0);
      set("NHF", 0.99);
      set("CEF", 0.99);
      set("NEF", 0.99);
      set("NCH", 0);
      set("nConstituents", 1);
      set("BTAG", false);
      
      // now set the return values for the ignored parts
      if ( quality_ == Top_SelV3 ) {
	set("pt", 35.0);   //put back to 30
	set("eta", 2.4);
      } 

       if ( quality_ == MEENAKSHI ) {
	set("pt", 20.0);
	set("eta", 2.4);
      } 

       if ( quality_ == DiLepton ) {
	 set("pt", 30.0);
	 set("eta", 2.5);
	 //set("btag", 0.679);
	 set("BTAG", false);
       } 

    } // end of version == DATA2010 cuts

  
    retInternal_ = getBitTemplate();
  }


 




   // 
  // Accessor from PAT jets
  // 
  bool operator()( const pat::Jet & jet, pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) {
      if ( jet.currentJECLevel() == "Uncorrected" ) 
	{
	  return data2010Cuts( jet, ret );	  
	}

      else{
	return data2010Cuts( jet.correctedJet("Uncorrected"), ret );
	//return data2010Cuts( jet, ret );
      }
    }
    else {
      return false;
    }
  }
  using Selector<pat::Jet>::operator();

  // 
  // Accessor from *CORRECTED* 4-vector, EMF, and Jet ID. 
  // This can be used with reco quantities. 
  // 
  bool operator()( reco::PFJet const & jet, 
		   pat::strbitset & ret )  
  {
    if ( version_ == DATA2010 ) return data2010Cuts( jet, ret );
    else {
      return false;
    }
  }





    // 
  // cuts based on craft 08 analysis. 
  // 
  bool data2010Cuts( reco::Jet const & jet,
		     pat::strbitset & ret) 
  {    

    ret.set(false);

    // cache some variables
    double chf = 0.0;
    double nhf = 0.0;
    double cef = 0.0;
    double nef = 0.0;
    int    nch = 0;
    int    nconstituents = 0;
    double pt = 0.0;
    double eta = 0.0;
    double btag = 0.0;

    // FIXME: unused variables
    (void)pt; // silence warning for now - clean up!
    (void)eta;

    // Have to do this because pat::Jet inherits from reco::Jet but not reco::PFJet
    reco::PFJet const * pfJet = dynamic_cast<reco::PFJet const *>(&jet);
    pat::Jet const * patJet = dynamic_cast<pat::Jet const *>(&jet);
    reco::BasicJet const * basicJet = dynamic_cast<reco::BasicJet const *>(&jet);

    if ( patJet != 0 ) {

      if ( patJet->isPFJet() ) {
	chf = patJet->chargedHadronEnergyFraction();
	nhf = patJet->neutralHadronEnergyFraction();
	cef = patJet->chargedEmEnergyFraction();
	nef = patJet->neutralEmEnergyFraction();
	nch = patJet->chargedMultiplicity();
	nconstituents = patJet->numberOfDaughters();
	pt = patJet->pt();
	eta = patJet->eta();	
	btag = patJet->bDiscriminator( "CombinedSecondaryVertex" );
      } 
      // Handle the special case where this is a composed jet for
      // subjet analyses
      else if ( patJet->isBasicJet() ) {
	double e_chf = 0.0;
	double e_nhf = 0.0;
	double e_cef = 0.0;
	double e_nef = 0.0;
	nch = 0;
	nconstituents = 0;

	for ( reco::Jet::const_iterator ibegin = patJet->begin(),
		iend = patJet->end(), isub = ibegin;
	      isub != iend; ++isub ) {
	  reco::PFJet const * pfsub = dynamic_cast<reco::PFJet const *>( &*isub );
	  e_chf += pfsub->chargedHadronEnergy();
	  e_nhf += (pfsub->neutralHadronEnergy() + pfsub->HFHadronEnergy());
	  e_cef += pfsub->chargedEmEnergy();
	  e_nef += pfsub->neutralEmEnergy();
	  nch += pfsub->chargedMultiplicity();
	  nconstituents += pfsub->numberOfDaughters();
	}
	double e = patJet->energy();
	if ( e > 0.000001 ) {
	  chf = e_chf / e;
	  nhf = e_nhf / e;
	  cef = e_cef / e;
	  nef = e_nef / e;
	} else {
	  chf = nhf = cef = nef = 0.0;
	}
      }
    } // end if pat jet
    
    else if ( pfJet != 0 ) {
      chf = pfJet->chargedHadronEnergyFraction();
      nhf = ( pfJet->neutralHadronEnergy() + pfJet->HFHadronEnergy() ) / pfJet->energy();
      cef = pfJet->chargedEmEnergyFraction();
      nef = pfJet->neutralEmEnergyFraction();
      nch = pfJet->chargedMultiplicity();
      nconstituents = pfJet->numberOfDaughters();
      pt = pfJet->pt();
      eta = pfJet->eta();
    } // end if PF jet
    
    // Handle the special case where this is a composed jet for
    // subjet analyses
    else if ( basicJet != 0 ) {
      double e_chf = 0.0;
      double e_nhf = 0.0;
      double e_cef = 0.0;
      double e_nef = 0.0;
      nch = 0;
      nconstituents = 0;
      
      for ( reco::Jet::const_iterator ibegin = basicJet->begin(),
	      iend = patJet->end(), isub = ibegin;
	    isub != iend; ++isub ) {
	reco::PFJet const * pfsub = dynamic_cast<reco::PFJet const *>( &*isub );
	e_chf += pfsub->chargedHadronEnergy();
	e_nhf += (pfsub->neutralHadronEnergy() + pfsub->HFHadronEnergy());
	e_cef += pfsub->chargedEmEnergy();
	e_nef += pfsub->neutralEmEnergy();
	nch += pfsub->chargedMultiplicity();
	nconstituents += pfsub->numberOfDaughters();
      }
      double e = basicJet->energy();
      if ( e > 0.000001 ) {
	chf = e_chf / e;
	nhf = e_nhf / e;
	cef = e_cef / e;
	nef = e_nef / e;
      }
    } // end if basic jet


    // Cuts for all |eta|:
    if ( ignoreCut("BTAG") || btag > cut("BTAG", double() ) ) passCut( ret, "BTAG");
    if ( ignoreCut("nConstituents") || nconstituents > cut("nConstituents", int() ) ) passCut( ret, "nConstituents");
    if ( ignoreCut("NEF")           || ( nef < cut("NEF", double()) ) ) passCut( ret, "NEF");
    if ( ignoreCut("NHF")           || ( nhf < cut("NHF", double()) ) ) passCut( ret, "NHF");    
    // Cuts for |eta| < 2.4:
    if ( ignoreCut("CEF")           || ( cef < cut("CEF", double()) && std::abs(jet.eta()) < 2.4 ) ) passCut( ret, "CEF");
    if ( ignoreCut("CHF")           || ( chf > cut("CHF", double()) && std::abs(jet.eta()) < 2.4 ) ) passCut( ret, "CHF");
    if ( ignoreCut("NCH")           || ( nch > cut("NCH", int())    && std::abs(jet.eta()) < 2.4 ) ) passCut( ret, "NCH");    
    
    if ( jet.pt()            >  cut("pt", double())     || ignoreCut("pt")      ) passCut(ret, "pt" );
    if ( std::abs(jet.eta()) <  cut("eta", double())     || ignoreCut("eta")      ) passCut(ret, "eta" );
    
	

    setIgnored( ret );
    return (bool)ret;
  }



  
  std::string getVersionString(Version_t ver){
    if (ver == DATA2010) return "DATA2010";
    else return "N_VERSIONS";
  }


  std::string getQualityString(Quality_t qual){
    if (qual == Top_SelV3) return "Top_SelV3";
    else if (qual == MEENAKSHI) return "MEENAKSHI";
    else if (qual == DiLepton) return "DiLepton";
    else return "N_QUALITY";
  }
  



 private: // member variables
  
  Version_t version_;
  Quality_t quality_;

  std::string legend;
};

#endif
