#ifndef LJMet_Com_interface_ElectronSelector_h
#define LJMet_Com_interface_ElectronSelector_h


/**
   \class    ElectronSelector ElectronSelector.h "LJMet/Com/interface/ElectronSelector.h"
   \brief    Electron selector for pat::Electron

   Selector functor for pat::Electron inspired by JetIDSelectionFunctor

   Please see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATSelectors
   for a general overview of the selectors. 

   \author Michael Segala, Gena Kukartsev
   \version  $Id: ElectronSelector.h,v 1.18 2012/09/06 01:11:39 kukartse Exp $
*/




#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>
class ElectronSelector : public Selector<pat::Electron>  {

public: // interface

    enum Version_t { DATA2010, N_VERSIONS };
    enum Quality_t { OFF, Top_SelV3, TOP_TIGHT, Top_Loose_SelV3, TOP_VERYLOOSE, EWK_Wenu, DiLepton, N_QUALITY };
  

    ElectronSelector( edm::ParameterSet const & parameters ) {

        // legend for identifying output messages
        legend = "[";
        legend.append("ElectronSelector");
        legend.append("]: ");

        std::string versionStr = parameters.getParameter<std::string>("version");
        std::string qualityStr = parameters.getParameter<std::string>("quality");
       
        Version_t version = N_VERSIONS;
        Quality_t quality = N_QUALITY;

        if ( versionStr == "DATA2010" ) {

            version = DATA2010;

            if      ( qualityStr == "OFF" )   quality = OFF;
            else if ( qualityStr == "Top_SelV3" ) quality = Top_SelV3;
            else if ( qualityStr == "Top_Loose_SelV3" ) quality = Top_Loose_SelV3;
            else if ( qualityStr == "TOP_TIGHT" ) quality = TOP_TIGHT;
            else if ( qualityStr == "TOP_VERYLOOSE" ) quality = TOP_VERYLOOSE;
            else if ( qualityStr == "EWK_Wenu" ) quality = EWK_Wenu;
            else if ( qualityStr == "DiLepton" ) quality = DiLepton;
                                 
        } else {
            throw cms::Exception("InvalidInput") << "Expect version to be one of DATA2010, ..." << std::endl;
        }

        initialize( version, quality );

    }

    ElectronSelector( Version_t version, Quality_t quality ) {
        initialize(version, quality);
    }
  
    void initialize( Version_t version, Quality_t quality )
    {

        version_ = version;
        quality_ = quality;

        std::cout << legend << "New selector initialized" << std::endl;
        std::cout << legend << "Version: " << getVersionString(version_) << std::endl;
        std::cout << legend << "Quality: " << getQualityString(quality_) << std::endl;

        push_back("Electron ID");
        push_back("ET");
        push_back("PT");
        push_back("ETA");
        push_back("NotEBEEGap");
        push_back("NotEBEEGapDiLepton");
        push_back("Ip");

        push_back("RelIso");
        push_back("MaxRelIso");

        push_back("Cic");
        push_back("NotConversion");
        push_back("NotConversionTight");
   
        push_back("MinDrEJet");
	push_back("vertex");    

        // all on by default
        electronId = "someElectronId";
        set("Electron ID", true);
        set("ET", true);
        set("PT", false);
        set("ETA", true);
        set("NotEBEEGap",true);
        set("RelIso", true);
        set("MaxRelIso", false);
        set("Ip", false);
        set("Cic", false);
        set("NotConversion",true);
        set("NotConversionTight",false);
        set("MinDrEJet",false);
        set("vertex",false);

	//Off except for dileptons
        set("NotEBEEGapDiLepton", false);

        // now set the return values for the ignored parts
        if ( quality_ == OFF ) {

            set("Electron ID", false);
            set("ET", false );
            set("ETA", false );
            set("RelIso", false );
            set("MaxRelIso", false);
	    set("NotEBEEGap", false);
	    set("NotEBEEGapDiLepton", false);
            set("Ip", false);
            set("Cic", false);
            set("NotConversion", false);
	}
    

        if ( quality_ == Top_SelV3 ) {

            set("Electron ID", false);
            set("ET", 15.0);
            set("ETA", 2.5);
            set("RelIso", 0.2);
            set("MaxRelIso", false);
	    set("NotEBEEGap",false);
	    set("MaxRelIso", false);
	    set("Ip", false);
	    set("Cic", false);
	    set("NotConversion",false);
        }

        if ( quality_ == Top_Loose_SelV3 ) {
	
            set("Electron ID", false);
            set("ET", 15.0 );
            set("ETA", 2.5 );
            set("NotEBEEGap",false);
            set("RelIso", 0.2 );
            set("MaxRelIso", false);
	    set("Ip", false);
            set("Cic", false);
            set("NotConversion",false);
            set("NotConversionTight",false);
            set("MinDrEJet",false);
            set("vertex",false);

        }



        if ( quality_ == TOP_TIGHT ) {

            set("Electron ID", false);
            set("ET", 35.0 );
            set("ETA", 2.5 );
            set("NotEBEEGap",true);
            set("RelIso", 0.2 );
            set("MaxRelIso", false);
            set("Ip", 0.02);
            set("Cic", true);
            set("NotConversion",false);
            set("NotConversionTight",true);
            set("MinDrEJet",0.3);
            set("vertex",1.0);

        }

        if ( quality_ == TOP_VERYLOOSE ) {

            set("Electron ID", false);
            set("ET", 15.0);
            set("ETA", false);
            set("RelIso", 0.2);
            set("MaxRelIso", false);

        }

        if ( quality_ == EWK_Wenu ) {
      
            set("Electron ID", false);
            set("ET", false );
            set("ETA", false );
            set("RelIso", false );
            set("MaxRelIso", false);
        }

        if ( quality_ == DiLepton ) {
	  
            set("Electron ID", true);
            set("ET", 20.0 );
            set("ETA", 2.5 );
            set("NotEBEEGap", false);
	    set("NotEBEEGapDiLepton", true);
            set("RelIso", 0.2 );
            set("MaxRelIso", false);
            set("Ip", 0.02);
            set("Cic", true);
            set("NotConversion", false);
        }

    
        retInternal_ = getBitTemplate();
    }


    // 
    // Accessor from PAT Electron
    //
 

    bool operator()( const pat::Electron & electron, pat::strbitset & ret )  
    {
      return operator()(electron,ret,0,0,0);
    }

    bool operator()( const pat::Electron & electron, pat::strbitset & ret, 
		     std::vector<edm::Ptr<pat::Jet> > * vp_jets,
		     std::vector<edm::Ptr<reco::Vertex> > * vp_PV,
		     double rho)  
    {
      if ( version_ == DATA2010 ) return data2010Cuts(electron,ret,vp_jets,vp_PV,rho);
        else {
            return false;
        }
    }
    
    // accessor from PAT electron without the ret
    //using Selector<pat::Electron>::operator();
  
    // 
    // cuts for 2010 early data
    // 
 
    bool data2010Cuts( const pat::Electron & electron,
                       pat::strbitset & ret, 
		       std::vector<edm::Ptr<pat::Jet> > * vp_jets = 0,
		       std::vector<edm::Ptr<reco::Vertex> > * vp_PV = 0,
		       double rho = 0) 
    {

        ret.set(false);

        double et      = electron.et();
        double pt      = electron.pt();
        double eta     = electron.eta();
        double sceta   = electron.superCluster()->eta();

        double dr03TkSumPt         = electron.dr03TkSumPt();
        double dr03EcalRecHitSumEt = electron.dr03EcalRecHitSumEt();
        double dr03HcalTowerSumEt  = electron.dr03HcalTowerSumEt();
        double dr03RelIso = (dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt) / et;

        double Iso = (electron.chargedHadronIso()+electron.neutralHadronIso()+electron.photonIso())/electron.pt();

        //bool CiC = ((int)electron.electronID("eidHyperTight1MC") & 1) == 1;       //turn this back on!!
        bool CiC = ((int)electron.electronID("eidTight") & 1) == 1;       //turn this back on!!
        double ip = electron.dB();
 
        //bool hadId = (((int)electron.electronID("eidHyperTight1MC") & 1)) == 1 && electron.isGsfCtfScPixChargeConsistent();
        bool hadId = (((int)electron.electronID("eidTight") & 1)) == 1 && electron.isGsfCtfScPixChargeConsistent();

        bool isEBEEGap = ( fabs(sceta)>1.4442 && fabs(sceta)<1.5660 );
        bool isNotEBEEGap = isEBEEGap ? false : true;

        bool isNotEBEEGapDiLepton = !electron.isEBEEGap();

        double convDist = electron.convDist();
        double convDcot = electron.convDcot();
        int nlosthits = electron.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
        bool isNotConversion = (nlosthits == 0 && (fabs(convDcot)>0.02 || fabs(convDist)>0.02));
        bool isNotConversionTight = (nlosthits == 0 && fabs(convDcot)>0.02 && fabs(convDist)>0.02);

	// jets
	double minDrEJet = -1.0;
	if (considerCut("MinDrEJet") && vp_jets){
	  // jet loop
	  for (std::vector<edm::Ptr<pat::Jet> >::const_iterator jet = vp_jets->begin();
	          jet != vp_jets->end(); ++jet) {
	  
	  double _dr = reco::deltaR(electron, **jet);
  
	  if (minDrEJet < 0.0 || minDrEJet > _dr) minDrEJet = _dr;

	  }
	}
	if ( minDrEJet == -1)  minDrEJet = 1;
     
	double Pv_match = 0.5;
	if(vp_PV) Pv_match = fabs( electron.vertex().z() - (*vp_PV->at(0)).z() );


        if ( hadId || ignoreCut("Electron ID") ) passCut(ret, "Electron ID");
        if ( et         >  cut("ET", double() ) || ignoreCut("ET")) passCut(ret, "ET");
        if ( pt         >  cut("PT", double() ) || ignoreCut("PT")) passCut(ret, "PT");
        if ( fabs(eta)  <  cut("ETA", double() ) || ignoreCut("ETA")) passCut(ret, "ETA");
        if ( isNotEBEEGap || ignoreCut("NotEBEEGap") ) passCut(ret, "NotEBEEGap");
        if ( isNotEBEEGapDiLepton || ignoreCut("NotEBEEGapDiLepton") ) passCut(ret, "NotEBEEGapDiLepton");
        if ( ip         < cut("Ip", double() ) || ignoreCut("Ip")) passCut(ret, "Ip");
        if ( Iso        <  cut("RelIso", double() ) || ignoreCut("RelIso")) passCut(ret, "RelIso");
        if ( dr03RelIso <  cut("MaxRelIso", double() ) || ignoreCut("MaxRelIso")) passCut(ret, "MaxRelIso");
        if ( CiC        || ignoreCut("Cic") ) passCut(ret, "Cic");
        if ( isNotConversion || ignoreCut("NotConversion") ) passCut(ret, "NotConversion");
        if ( isNotConversionTight || ignoreCut("NotConversionTight") ) passCut(ret, "NotConversionTight");
        if ( minDrEJet >  cut("MinDrEJet", double()) || ignoreCut("MinDrEJet")) passCut(ret, "MinDrEJet");
        if ( Pv_match   <  cut("vertex", double()) || ignoreCut("vertex")) passCut(ret, "vertex");
        
	setIgnored( ret );    

        return (bool)ret;
    }

    std::string getVersionString(Version_t ver){
        if (ver == DATA2010) return "DATA2010";
        else return "N_VERSIONS";
    }

    std::string getQualityString(Quality_t qual){
        if (qual == OFF) return "OFF";
        else if (qual == Top_SelV3) return "Top_SelV3";
        else if (qual == Top_Loose_SelV3) return "Top_Loose_SelV3";
        else if (qual == TOP_TIGHT) return "TOP_TIGHT";
        else if (qual == TOP_VERYLOOSE) return "TOP_VERYLOOSE";
        else if (qual == EWK_Wenu) return "EWK_Wenu";
        else if (qual == DiLepton) return "DiLepton";
        else return "N_QUALITY";
    }

  
  

private: // member variables
  
    Version_t version_;
    Quality_t quality_;

    std::string legend;

    std::string electronId;
};

#endif

