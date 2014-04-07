#ifndef LJMet_Com_interface_MuonSelectionFunctor_h
#define LJMet_Com_interface_MuonSelectionFunctor_h

/**
  \class    MuonSelectionFunctor MuonSelectionFunctor.h "LJMet/Com/interface/MuonSelectionFunctor.h"
  \brief    Muon selector for pat::Muon

  Selector functor for pat::Muon inspired by JetIDSelectionFunctor

  Please see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATSelectors
  for a general overview of the selectors. 

  \author Michael Segala
  \author Gena Kukartsev
  \version  $Id: MuonSelectionFunctor.h,v 1.47 2012/08/01 22:52:07 kukartse Exp $
*/




#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>
class MuonSelectionFunctor : public Selector<pat::Muon>  {

 public: // interface

  enum Version_t { DATA2010, N_VERSIONS };
  enum Quality_t { OFF, Top_SelV3, Top_Loose_SelV3, Top_SelV3_36X, TOP_TIGHT, TOP_LOOSE, EWK_Wmunu, EWK_Zmumu, EWK_Zmumu_loose, N_QUALITY, MEENAKSHI, TOP_BTAG, TOP_2012, DiLepton };
  

  MuonSelectionFunctor( edm::ParameterSet const & parameters ) {

    // legend for identifying output messages
    legend = "[";
    legend.append("MuonSelector");
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
      else if ( qualityStr == "Top_SelV3_36X" ) quality = Top_SelV3_36X;
      else if ( qualityStr == "TOP_LOOSE" ) quality = TOP_LOOSE;
      else if ( qualityStr == "TOP_TIGHT" ) quality = TOP_TIGHT;
      else if ( qualityStr == "EWK_Wmunu" ) quality = EWK_Wmunu;
      else if ( qualityStr == "EWK_Zmumu" ) quality = EWK_Zmumu;
      else if ( qualityStr == "EWK_Zmumu_loose" ) quality = EWK_Zmumu_loose;
      else if ( qualityStr == "MEENAKSHI" ) quality = MEENAKSHI;
      else if ( qualityStr == "TOP_BTAG" ) quality = TOP_BTAG;
      else if ( qualityStr == "TOP_2012" ) quality = TOP_2012;
      else if ( qualityStr == "DiLepton" ) quality = DiLepton;
                                 
    } else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of DATA2010" << std::endl;
    }

    initialize( version, quality );

  }

  MuonSelectionFunctor( Version_t version, Quality_t quality ) {
    initialize(version, quality);
  }
  
  void initialize( Version_t version, Quality_t quality )
  {
    version_ = version;
    quality_ = quality;

    std::cout << legend << "New selector initialized" << std::endl;
    std::cout << legend << "Version: " << getVersionString(version_) << std::endl;
    std::cout << legend << "Quality: " << getQualityString(quality_) << std::endl;

    push_back("Global muon");
    push_back("PF muon");
    push_back("Tracker muon");
    push_back("pT");
    push_back("eta");
    push_back("RelIso");
    push_back("Global muon prompt tight");
    push_back("MinInnerTrackHits");
    push_back("MinDrMuJet");
    push_back("MaxIpBs2d");

    push_back("MaxD0"); // uses old track IP interface, for 36X

    push_back("TrackIso");
    push_back("nChi2");  
    push_back("TrackerHits");
    push_back("MuonHits");
    push_back("PixelHits");
    push_back("nMatches");

    // legacy cuts, kept for backwards compatibility
    // do not use these for new cut definitions
    push_back("Muon ID");
    push_back("RelCombinedIso");
    push_back("ECalVeto");
    push_back("HCalVeto");
    push_back("D0");
    push_back("D0 significance");
      
    push_back("vertex");
    push_back("Tracker Layers");
    


    // all on by default
    set("Global muon",               true);
    set("PF muon",                   false);
    set("Tracker muon",              true);
    set("pT",                       15.0);
    set("eta",                       2.1);
    set("RelIso",                    0.05);
    set("Global muon prompt tight",  true);
    set("MinInnerTrackHits",        11);
    set("MinDrMuJet",                0.3);
    set("MaxIpBs2d",                 0.02); //cm

    set("MaxD0",                     0.02); //cm

    set("TrackIso",     3.0);
    set("nChi2",       10.0);
    set("TrackerHits", 10);
    set("MuonHits",     0);
    set("PixelHits",    0);
    set("nMatches",     2);

    // legacy cuts, kept for backwards compatibility
    // do not use these for new cut definitions
    // switching off by default
    set("Muon ID",         false);
    set("RelCombinedIso",  0.15);
    set("ECalVeto",        4.0);
    set("HCalVeto",        6.0);
    set("D0",              0.2);
    set("D0 significance", 0.02);
    set("vertex", 1);
    set("Tracker Layers", false);



    // now set the return values for the ignored parts
    if ( quality_ == OFF ) {
      set("Global muon",               false);
      set("Tracker muon",              false);
      set("pT",                        false);
      set("eta",                       false);
      set("RelIso",                    false);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         false);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }
    
    if ( quality_ == TOP_TIGHT ) {
      set("Global muon",               false);
      set("Tracker muon",              false);
      set("pT",                       15.0);
      set("eta",                       2.1);
      set("RelIso",                    0.05);
      set("Global muon prompt tight",  true);
      set("MinInnerTrackHits",        11);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",       10.0);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        true);
      set("RelCombinedIso", false);
      set("ECalVeto",       4.0);
      set("HCalVeto",       6.0);
      set("D0",             0.2);
      set("D0 significance",false);
      set("vertex", false);
    }


    if ( quality_ == Top_SelV3 ) {
      set("Global muon",               true);
      set("Tracker muon",              true);
      set("pT",                       20.0);
      set("eta",                       2.1);
      set("RelIso",                    0.05);
      set("Global muon prompt tight",  true);
      set("MinInnerTrackHits",        11);
      set("MinDrMuJet",                0.3);
      set("MaxIpBs2d",                 0.02); //cm
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }

    if ( quality_ == TOP_BTAG ) {

      set("Global muon",               true);
      set("Tracker muon",              true);
      set("pT",                       30.0);        //change back to 30
      set("eta",                       2.1);
      set("RelIso",                    0.15);     
      //set("RelIso",                    false);     
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",        10); 
      set("MinDrMuJet",                0.3);
      set("MaxIpBs2d",                 0.02); //cm
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        10);
      set("TrackerHits",  false);
      set("MuonHits",     0);
      set("PixelHits",    1);
      set("nMatches",     1);
      set("Muon ID",        true);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", 1.0);
    }

    if ( quality_ == TOP_2012 ) {

      set("Global muon",               true);
      set("PF muon",                   true);
      set("Tracker muon",              true);
      set("pT",                       20.0);        //change back to 30
      set("eta",                       2.1);
      set("RelIso",                    0.125);     
      //set("RelIso",                    false);     
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",        10); 
      set("MinDrMuJet",                0.3);
      set("MaxIpBs2d",                 0.02); //cm
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        10);
      set("TrackerHits",  false);
      set("MuonHits",     0);
      set("PixelHits",    1);
      set("nMatches",     1);
      set("Muon ID",        true);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", 1.0);
      set("Tracker Layers", 5);
    }

    if ( quality_ == TOP_LOOSE || quality_ == Top_Loose_SelV3 ) {
      set("Global muon",               true);
      set("Tracker muon",              false);
      set("pT",                       10.0);
      set("eta",                       2.5);
      set("RelIso",                    0.2);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         false);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }


    if ( quality_ == Top_SelV3_36X ) {
      set("Global muon",               true);
      set("Tracker muon",              true);
      set("pT",                       20.0);
      set("eta",                       2.1);
      set("RelIso",                    0.05);
      set("Global muon prompt tight",  true);
      set("MinInnerTrackHits",        11);
      //set("MinDrMuJet",                0.3);
      set("MinDrMuJet",                false); //TRUNING OFF FOR PFJETS
      set("MaxIpBs2d",                 false); //cm
      set("MaxD0",                     0.02); //cm
      set("TrackIso",     false);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }


    if ( quality_ == MEENAKSHI ) {
      set("Global muon",               true);
      set("Tracker muon",              true);
      set("pT",                       2.0);
      set("eta",                       2.5);
      set("RelIso",                    false);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",        7);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false); //cm
      set("MaxD0",                     false); //cm
      set("TrackIso",     false);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }



    if ( quality_ == EWK_Wmunu ) {
      set("Global muon",               false);
      set("Tracker muon",              false);
      set("pT",                       15);
      set("eta",                       2.1);
      set("RelIso",                    false);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         false);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        true);
      set("TrackerHits", 10);
      set("MuonHits",     int(0));
      set("PixelHits",    int(0));
      set("nMatches",     2);
      set("Muon ID",        true);
      set("RelCombinedIso", 0.15);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }
    



    if ( quality_ == EWK_Zmumu ) {
      set("Global muon",               false);
      set("Tracker muon",              false);
      set("pT",                       15.0);
      set("eta",                       2.1);
      set("RelIso",                    0.05);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         11);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",       10.0);
      set("TrackerHits",  int(0));
      set("MuonHits",     int(0));
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        false);
      set("RelCombinedIso", false);
      set("ECalVeto",       true);
      set("HCalVeto",       true);
      set("D0",             false);
      set("D0 significance", 3.0);
      set("vertex", false);
    }
    
    if ( quality_ == EWK_Zmumu_loose ) {
      set("Global muon",               false);
      set("Tracker muon",              false);
      set("pT",                       15.0);
      set("eta",                       2.1);
      set("RelIso",                    false);
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         false);
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 false);
      set("MaxD0",                     false);
      set("TrackIso",     3.0);
      set("nChi2",        false);
      set("TrackerHits",  false);
      set("MuonHits",     false);
      set("PixelHits",    false);
      set("nMatches",     false);
      set("Muon ID",        true);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);
    }


    if ( quality_ == DiLepton ) {
      
      set("Global muon",               true);
      set("Tracker muon",              true);
      set("pT",                       20.0);
      set("eta",                       2.4);
      set("RelIso",                    false);     
      set("Global muon prompt tight",  false);
      set("MinInnerTrackHits",         10); 
      set("MinDrMuJet",                false);
      set("MaxIpBs2d",                 0.02); //cm
      set("MaxD0",                     false);
      set("TrackIso",     false);
      set("nChi2",        10);
      set("TrackerHits",  false);
      set("MuonHits",     0);
      set("PixelHits",    1);
      set("nMatches",     1);
      set("Muon ID",        true);
      set("RelCombinedIso", false);
      set("ECalVeto",       false);
      set("HCalVeto",       false);
      set("D0",             false);
      set("D0 significance",false);
      set("vertex", false);

    }
    
    retInternal_ = getBitTemplate();
  }


  bool operator()( const pat::Muon & muon, pat::strbitset & ret )  
  {
    return operator()(muon,ret,0,0,0);
  }

  // this is the actual operator
  bool operator()( const pat::Muon & muon, pat::strbitset & ret,
		   std::vector<edm::Ptr<pat::Jet> > * vp_jets,
		   std::vector<edm::Ptr<reco::Vertex> > * vp_PV,
		   double rho)  
  {
    //if ( version_ == DATA2010 ) return data2010Cuts( muon.p4(), ret );
    if ( version_ == DATA2010 ) return data2010Cuts(muon,ret,vp_jets,vp_PV, rho);
    else {
      return false;
    }
  }
    
  // 
  // cuts inherited from cuts for 2010 early data and modified (a lot)
  // 
 
  bool data2010Cuts( const pat::Muon & muon,
		     pat::strbitset & ret,
		     std::vector<edm::Ptr<pat::Jet> > * vp_jets = 0,
		     std::vector<edm::Ptr<reco::Vertex> > * vp_PV = 0,
		     double rho = 0) 
  {

    ret.set(false);

    double pt      = muon.pt();// - rho*0.4*0.4;
    double eta     = muon.eta();
    //double hcalIso = muon.hcalIso();
    //double ecalIso = muon.ecalIso();
    //double trkIso  = muon.trackIso();
    //double relIso = (ecalIso + hcalIso + trkIso) / pt;
    double relIso = (muon.chargedHadronIso()+muon.neutralHadronIso()+muon.photonIso())/ pt;
    double ecalVeto = muon.isolationR03().emVetoEt;
    double hcalVeto = muon.isolationR03().hadVetoEt;
 
    double SumPt = muon.isolationR03().sumPt;
    double EmEt   = muon.isolationR03().emEt; 
    double HcalEt = muon.isolationR03().hadEt;    
    double relCombinedIso = ( SumPt + EmEt + HcalEt )/pt;

    //double corr_d0;
    double corr_d0     = muon.dB();
    double corr_d0_err = muon.edB();
    double d0_significance = fabs(corr_d0/corr_d0_err);
   

    double chi2 = -1.0;
    double ndof = -1.0;
    double nChi2 = -1.0;

    // FIXME: unused variable nChi2
    (void)nChi2;

    double normChi2 = -1.0; // alternative version of normalized Chi2, directly from global track
    double MinInnerTrackHits = -1.0;
    double muonHits = -1.0;
    double trackerHits = -1.0;
    int pixelHits = -1;
    int nMatches; 
    double minDrMuJet = -1.0;
    //double maxIpBs2d = -1.0;
    double maxIpBs2d = muon.dB(pat::Muon::BS2D);
    int trackerLayers = -1.0;
    
    //reco::Vertex const & pv = vp_PV->at(0);
    
    double Pv_match = 0.5;
    if(vp_PV) Pv_match = fabs( muon.vertex().z() - (*vp_PV->at(0)).z() );

    //nMatches = muon.numberOfMatches(); 
    nMatches = muon.numberOfMatchedStations();

    // now the global track parameters
    reco::TrackRef muTrack = muon.globalTrack();
    if (  muTrack.isNonnull() && muTrack.isAvailable() ) {
      chi2  = muTrack->chi2();
      ndof  = muTrack->ndof();
      normChi2  = muTrack->normalizedChi2();
      muonHits = muTrack->hitPattern().numberOfValidMuonHits();
      trackerLayers = muTrack->hitPattern().trackerLayersWithMeasurement();
      
      if ( ndof  != 0 ){
	nChi2  = chi2 / ndof;
      }
      else nChi2 = -1.0;
    }

    // now inner track parameters
    muTrack = muon.innerTrack();
    if (  muTrack.isNonnull() && muTrack.isAvailable() ) {
      MinInnerTrackHits = muon.innerTrack()->numberOfValidHits();
      pixelHits = muTrack->hitPattern().pixelLayersWithMeasurement();
    }
    
    
    //std::cout << pt << " " << SumPt << " " << corr_d0 << " " << ndof  << " " << pixelHits  << " " <<    MinInnerTrackHits << std::endl;


    while(1){
      if ( muon.isGlobalMuon() || ignoreCut("Global muon")  )  passCut(ret, "Global muon" );
      else break;

      if ( muon.isPFMuon() || ignoreCut("PF muon")  )  passCut(ret, "PF muon" );
      else break;

      if ( muon.isTrackerMuon() || ignoreCut("Tracker muon")  )  passCut(ret, "Tracker muon" );
      else break;

      if ( pt                >  cut("pT", double()) || ignoreCut("pT"))                                passCut(ret, "pT");
      else break;

      if ( fabs(eta)         <  cut("eta", double()) || ignoreCut("eta"))                              passCut(ret, "eta");
      else break;


      // jets
      if (considerCut("MinDrMuJet") && vp_jets){
	// jet loop
	for (std::vector<edm::Ptr<pat::Jet> >::const_iterator 
	       jet = vp_jets->begin();
	     jet != vp_jets->end(); ++jet){
	  
	  double _dr = reco::deltaR(muon, **jet);
  
	  if (minDrMuJet < 0.0 || minDrMuJet > _dr) minDrMuJet = _dr;

	}
      }
      
      if ( minDrMuJet == -1)  minDrMuJet = 1;

      if ( relIso            <  cut("RelIso", double()) || ignoreCut("RelIso"))                        passCut(ret, "RelIso");
      else break;

      if ( (normChi2<10.0 && muonHits>0) || ignoreCut("Global muon prompt tight")  )  passCut(ret, "Global muon prompt tight" );
      else break;

      if ( MinInnerTrackHits >  cut("MinInnerTrackHits", int()) || ignoreCut("MinInnerTrackHits"))     passCut(ret, "MinInnerTrackHits");    
      else break;

      if ( trackerLayers >  cut("Tracker Layers", int()) || ignoreCut("Tracker Layers"))     passCut(ret, "Tracker Layers");    
      else break;


      if ( minDrMuJet        >  cut("MinDrMuJet", double()) || ignoreCut("MinDrMuJet"))                passCut(ret, "MinDrMuJet");
      else break;

      if ( fabs(maxIpBs2d)   <  cut("MaxIpBs2d", double()) || ignoreCut("MaxIpBs2d"))                  passCut(ret, "MaxIpBs2d");
      else break;

      if ( Pv_match   <  cut("vertex", double()) || ignoreCut("vertex"))                  passCut(ret, "vertex");
      else break;

      if ( fabs(corr_d0)     <  cut("MaxD0", double()) || ignoreCut("MaxD0"))                          passCut(ret, "MaxD0");
      else break;

      if ( SumPt             <  cut("TrackIso", double()) || ignoreCut("TrackIso"))                    passCut(ret, "TrackIso");
      else break;

      if ( normChi2             <  cut("nChi2", double()) || ignoreCut("nChi2"))                          passCut(ret, "nChi2");
      else break;

      if ( trackerHits       >  cut("TrackerHits", int()) || ignoreCut("TrackerHits"))                 passCut(ret, "TrackerHits");
      else break;

      if ( muonHits          >  cut("MuonHits", int()) || ignoreCut("MuonHits"))                       passCut(ret, "MuonHits");
      else break;

      if ( pixelHits         >=  cut("PixelHits", int()) || ignoreCut("PixelHits"))                     passCut(ret, "PixelHits");
      else break;

      if ( nMatches         >  cut("nMatches", int()) || ignoreCut("nMatches"))                       passCut(ret, "nMatches");
      else break;

      if ( (muon.isTrackerMuon() == true && muon.isGlobalMuon() == true) || ignoreCut("Muon ID")  )  passCut(ret, "Muon ID" );
      else break;

      if ( relCombinedIso    <  cut("RelCombinedIso", double()) || ignoreCut("RelCombinedIso"))        passCut(ret, "RelCombinedIso");  
      else break;

      if ( ecalVeto          <  cut("ECalVeto",double()) || ignoreCut("ECalVeto"))                     passCut(ret, "ECalVeto");
      else break;

      if ( hcalVeto          <  cut("HCalVeto",double()) || ignoreCut("HCalVeto"))                     passCut(ret, "HCalVeto");
      else break;

      if ( fabs(corr_d0)     <  cut("D0", double()) || ignoreCut("D0"))                                passCut(ret, "D0");
      else break;

      if ( d0_significance   <  cut("D0 significance", double()) || ignoreCut("D0 significance"))      passCut(ret, "D0 significance");
      else break;
    
      break;
    
    }
    
    setIgnored( ret );    

    return (bool)ret;
  }
  
  
  bool data2010Cuts( reco::Candidate::LorentzVector const & correctedP4, 
		     pat::strbitset & ret) 
  {

    ret.set(false);

    // cache some variables
    //double corrPt = correctedP4.pt();

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
    else if (qual == Top_SelV3_36X) return "Top_SelV3_36X";
    else if (qual == TOP_TIGHT) return "TOP_TIGHT";
    else if (qual == TOP_LOOSE) return "TOP_LOOSE";
    else if (qual == EWK_Wmunu) return "EWK_Wmunu";
    else if (qual == EWK_Zmumu) return "EWK_Zmumu";
    else if (qual == EWK_Zmumu_loose) return "EWK_Zmumu_loose";
    else if (qual == MEENAKSHI) return "MEENAKSHI";
    else if (qual == TOP_BTAG) return "TOP_BTAG";
    else if (qual == TOP_2012) return "TOP_2012";
    else if (qual == DiLepton) return "DiLepton";
    else return "N_QUALITY";
  }



 private: // member variables

  std::string legend;
  
  Version_t version_;
  Quality_t quality_;
};

#endif
  
  
