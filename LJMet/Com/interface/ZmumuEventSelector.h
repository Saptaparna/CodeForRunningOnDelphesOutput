#ifndef LJMet_Com_interface_ZmumuEventSelector_h
#define LJMet_Com_interface_ZmumuEventSelector_h

/*   FWLite PAT analyzer-selector for MET commissioning in Z events
*/

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "LJMet/Com/interface/MuonSelectionFunctor.h"
#include "LJMet/Com/interface/ElectronSelector.h"

//#include "LJMet/Com/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

#include "LJMet/Com/interface/MetSelectionFunctor.h"
#include "LJMet/Com/interface/PVSelector.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/PatCandidates/interface/TriggerObject.h"

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

///-------------------------
/// DRIVER FUNCTION
///-------------------------

// -*- C++ -*-

// CMS includes
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "LJMet/Com/interface/LjmetEventSelector.h"


// Root includes
#include "TROOT.h"

using namespace std;


class ZmumuEventSelector : public LjmetEventSelector {
public:


  ZmumuEventSelector( edm::ParameterSet const &  muonSelParams,
		    edm::ParameterSet const & looseMuonSelParams,
		    edm::ParameterSet const & electronSelParams,
		    edm::ParameterSet const & looseElectronSelParams,
		    edm::ParameterSet const & jetSelParams,
		    edm::ParameterSet const & metSelParams,
		    edm::ParameterSet const & pvSelParams,
		    edm::ParameterSet const & params ) : LjmetEventSelector( muonSelParams,	  
									     looseMuonSelParams, 
									     electronSelParams,  
									     looseElectronSelParams,  
									     jetSelParams,	  
									     metSelParams,	  
									     pvSelParams,	  
									     params ){
    // this event selector is for top
    selection = EWK_Zmumu;
    std::string versionStr = params.getParameter<std::string>("version");
    
    if      ( versionStr == "Off" )   version = Off;
    else if ( versionStr == "Loose" ) version = Loose;
    else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of Off, Loose" << std::endl;
    }
    
    initialize( version );
  }
    
    

  ~ZmumuEventSelector();
  
  

  void initialize(Version_t version);



  // main method where the cuts are applied
  bool operator()( edm::EventBase const & event, pat::strbitset & ret);



protected:

};




ZmumuEventSelector::~ZmumuEventSelector(){}




void ZmumuEventSelector::initialize( Version_t version){
    
  push_back("No selection");
  set("No selection");
    
    std::cout << "Applying EWK_Zmumu selection" << std::endl;
      
    push_back("HLT_L1Mu9");
    push_back("HLT_L2Mu9");
    push_back("HLT_Mu9");
    push_back("Trigger names"); // activate to get a list of available trigger names
    push_back("Trigger cuts");
      
    //push_back("PV NDOF");
    //push_back("PV Z");
    //push_back("PV RHO");
    push_back("Primary vertex cuts");
      
    push_back("Dimuon ID");
    push_back("Dimuon pT");
    push_back("Dimuon eta");
    push_back("Dimuon iso");
    push_back("Dimuon");
    push_back("Tight mu multipl", 1);
    push_back("Min loose mu multipl", 1);
    push_back("Max loose mu multipl", 1);
    push_back("Muon cuts");
      
    push_back("Electron veto");
    push_back("Electron cuts");
      
    push_back("One jet or more");
    push_back("Two jets or more");
    push_back("Three jets or more");
    push_back("Min jet multiplicity", 4);
    push_back("Jet cuts");
      
    push_back("MET cuts");
      
    push_back("Opposite charges");
    push_back("Min dimuon mass", 0.0);
    push_back("Muon matched to HLT");
    push_back("Z candidate cuts");
      
      
    // all on by default
    set("HLT_L1Mu9");
    set("HLT_L2Mu9");
    set("HLT_Mu9");
    set("Trigger names", false); // this is a flag, not a cut
    set("Trigger cuts");
    //set("PV NDOF", false);
    //set("PV Z", false);
    //set("PV RHO", false);
    set("Primary vertex cuts");
    set("Dimuon ID");
    set("Dimuon pT");
    set("Dimuon eta");
    set("Dimuon iso");
    set("Dimuon");
    set("Tight mu multipl");
    set("Min loose mu multipl");
    set("Max loose mu multipl");
    set("Muon cuts");
    set("Electron veto");
    set("Electron cuts");
    set("One jet or more");
    set("Two jets or more");
    set("Three jets or more");
    set("Min jet multiplicity");
    set("Jet cuts");
    set("MET cuts");
    set("Opposite charges");
    set("Min dimuon mass");
    set("Muon matched to HLT", "HLT_Mu9");
    set("Z candidate cuts");

    // loose Zmumu selection (MET commissioning)
    if (version == Loose){
      std::cout << "Selection version: Loose" << std::endl;

      set("HLT_L1Mu9", false);
      set("HLT_L2Mu9", false);
      set("HLT_Mu9", true);
      set("Trigger names", false);
      set("Trigger cuts");
	
      set("Primary vertex cuts");
	
      set("Dimuon ID");
      set("Dimuon pT");
      set("Dimuon eta");
      set("Dimuon iso");
      set("Dimuon");
      set("Tight mu multipl", false);
      set("Min loose mu multipl", false);
      set("Max loose mu multipl", false);
      set("Muon cuts");
	
      set("Electron veto", false);
      set("Electron cuts", false);
	
      set("One jet or more", false);
      set("Two jets or more", false);
      set("Three jets or more", false);
      set("Min jet multiplicity", false);
      set("Jet cuts", true);
	
      set("MET cuts", true);
	
      set("Opposite charges", true);
      set("Min dimuon mass", 60.0);
      set("Muon matched to HLT", "HLT_Mu9");
      set("Z candidate cuts");
    }
    else if (version == Off){
      std::cout << "Selection version: Off" << std::endl;

      set("HLT_L1Mu9", false);
      set("HLT_L2Mu9", false);
      set("HLT_Mu9", false);
      set("Trigger names", false);
      set("Trigger cuts", false);
      set("Primary vertex cuts", false);
      set("Dimuon ID", false);
      set("Dimuon pT", false);
      set("Dimuon eta", false);
      set("Dimuon iso", false);
      set("Dimuon", false);
      set("Tight mu multipl", false);
      set("Min loose mu multipl", false);
      set("Max loose mu multipl", false);
      set("Muon cuts", false);
      set("Electron veto", false);
      set("Electron cuts", false);
      set("One jet or more", false);
      set("Two jets or more", false);
      set("Three jets or more", false);
      set("Min jet multiplicity", false);
      set("Jet cuts", false);
      set("MET cuts", false);
      set("Opposite charges", false);
      set("Min dimuon mass", false);
      set("Muon matched to HLT", false);
      set("Z candidate cuts", false);
    }

    muonMatchedToHlt_firstpass = true;

} // initialize()





bool ZmumuEventSelector::operator()( edm::EventBase const & event, pat::strbitset & ret){

  pat::strbitset retMuon      = muonSel_->getBitTemplate();
  pat::strbitset retDimuon      = muonSel_->getBitTemplate();
  pat::strbitset retLooseMuon = looseMuonSel_->getBitTemplate();
  pat::strbitset retElectron  = electronSel_->getBitTemplate();
  pat::strbitset retJet       = jetSel_->getBitTemplate();
  pat::strbitset retMet       = metSel_->getBitTemplate();
  pat::strbitset retPv        = pvSel_->getBitTemplate();
     
  _dimuon_mass = -1.0; // unphysical default


  while(1){ // standard infinite while loop trick to avoid nested ifs

    passCut(ret, "No selection");
      
    //
    //_____ Trigger cuts __________________________________
    //
    _hlt.clear();
    //std::map<std::string, int> _iHlt;
    if ( considerCut("Trigger cuts") ) {
      event.getByLabel( triggerSrc_, hltresults_ );

      for ( vector<pat::TriggerPath>::const_iterator trig = hltresults_->begin();
	    trig!=hltresults_->end();
	    ++trig){
	_hlt[trig->name()]=trig->wasAccept();
	//_iHlt[trig->name()] = trig -> indexModule(trig->name());

	if ( considerCut("Trigger names") ) {
	  std::cout << trig->name() << std::endl;
	}
      }

      // if we are only dumping trigger names, no need to go further
      if ( considerCut("Trigger names") ) return 0;

      //HLT_Mu9
      //HLT_L1_HFtech
      //HLT_L2Mu9
      //if ( _hlt["HLT_Mu9"] == true  ) {
      //  passTrig = true;
      //}

      if ( ignoreCut("HLT_Mu9") || _hlt["HLT_Mu9"] == true ) passCut(ret, "HLT_Mu9");
      else break;

      if ( ignoreCut("HLT_L2Mu9") || _hlt["HLT_L2Mu9"] == true ) passCut(ret, "HLT_L2Mu9");
      else break;

      if ( ignoreCut("HLT_L1Mu9") || _hlt["HLT_L1Mu9"] == true ) passCut(ret, "HLT_L1Mu9");
      else break;

      passCut(ret, "Trigger cuts"); // trigger cuts total
	
    } // end of trigger cuts
      

      //
      //_____ Primary vertex cuts __________________________________
      //
    if ( considerCut("Primary vertex cuts") ) {
      retPv.set(false);
      bool pass = (*pvSel_)( event, retPv );

      if (pass){}
      else break;

      passCut(ret, "Primary vertex cuts"); // PV cuts total

    } // end of PV cuts
      

      //
      //_____ Muon cuts __________________________________
      //      
    int nMuons = 0;
    int nLooseMuons = 0;
    int nMuonsId = 0;
    int nMuonsPt = 0;
    int nMuonsEta = 0;
    int nMuonsIso = 0;
    if ( considerCut("Muon cuts") ) {
      event.getByLabel( muonSrc_, h_muons_ );

      // loop over muons
      int muonIndex = 0;
      for ( vector<pat::Muon>::const_iterator mu = h_muons_->begin();
	    mu != h_muons_->end(); mu++){
	  
	retMuon.set(false);
	bool pass = (*muonSel_)( *mu, retMuon );

	if (retMuon.test("Muon ID")){
	  ++nMuonsId;
	  //if (nMuonsId == 1) muon0_ = edm::Ptr<pat::Muon>( h_muons_, muonIndex);
	  //if (nMuonsId == 2) muon1_ = edm::Ptr<pat::Muon>( h_muons_, muonIndex);
	}
	if (retMuon.test("PT")) ++nMuonsPt;
	if (retMuon.test("ETA")) ++nMuonsEta;
	if (retMuon.test("TrackIso")) ++nMuonsIso;
	  
	if ( pass && nMuons!=cut("Tight mu multipl",int()) ){ // looking for exactly N tight muons
	  ++nMuons;
	  // save the first good muon
	  if (nMuons == 1) muon0_ = edm::Ptr<pat::Muon>( h_muons_, muonIndex);
	}
	else{ // check for loose muon
	  retLooseMuon.set(false);
	  bool pass_loose = (*looseMuonSel_)( *mu, retLooseMuon );
	    
	  if ( pass_loose ){
	    ++nLooseMuons;
	    // save the first loose muon as the second muon in the event
	    if (nLooseMuons == 1) muon1_ = edm::Ptr<pat::Muon>( h_muons_, muonIndex);
	  }
	}

	++muonIndex;
      } // end of muon loop

      if( nMuonsId > 1 || ignoreCut("Dimuon ID") ) passCut(ret, "Dimuon ID");
      else break;
	
      if( nMuonsPt > 1 || ignoreCut("Dimuon pT") ) passCut(ret, "Dimuon pT");
      else break;
	
      if( nMuonsEta > 1 || ignoreCut("Dimuon eta") ) passCut(ret, "Dimuon eta");
      else break;
	
      if( nMuonsIso > 1 || ignoreCut("Dimuon iso") ) passCut(ret, "Dimuon iso");
      else break;
	
      if( nMuons+nLooseMuons == 2 || ignoreCut("Dimuon") ) passCut(ret, "Dimuon");
      else break;
	
      if( nMuons == cut("Tight mu multipl", int()) ) passCut(ret, "Tight mu multipl");
      else break;
	
      if( nLooseMuons >= cut("Min loose mu multipl", int()) ) passCut(ret, "Min loose mu multipl");
      else break;

      if( nLooseMuons <= cut("Max loose mu multipl", int()) ) passCut(ret, "Max loose mu multipl");
      else break;

      passCut(ret, "Muon cuts"); // muon cuts total

    } // end of muon cuts
      
      
      //
      //_____ Electron cuts __________________________________
      //      
    if ( considerCut("Electron cuts") ) {
      event.getByLabel( electronSrc_, h_electrons_ );
	
      // loop over electrons
      int nElectrons = 0;
      for ( vector<pat::Electron>::const_iterator el = h_electrons_->begin();
	    el != h_electrons_->end(); el++){
	retElectron.set(false);
	bool pass = (*electronSel_)( *el, retElectron );
	  
	if ( pass ){
	  ++nElectrons;
	  // save the first good electron
	  if (nElectrons == 1) electron0_ = edm::Ptr<pat::Electron>( h_electrons_, 0);
	}
      } // end of loop over electron
	
      if( nElectrons == 0 || ignoreCut("Electron veto") ) passCut(ret, "Electron veto");
      else break;

      passCut(ret, "Electron cuts"); // electron cuts total

    } // end of electron cuts


      //
      //_____ Jet cuts __________________________________
      //      
    if ( considerCut("Jet cuts") ) {
      event.getByLabel( jetSrc_, h_jets_ );
	
      // loop over jets
      int nJets = 0;
      int nGoodJets = 0;
      good_jets_ . clear();
      for (vector<pat::Jet>::const_iterator jet = h_jets_->begin();
	   jet != h_jets_->end(); jet++){
	  
	retJet.set(false);
	bool pass = (*jetSel_)( *jet, retJet );
	  
	  
	if ( pass &&
	     jet->pt() > 30.0 &&
	     fabs(jet->eta()) < 2.4){
	  ++nJets;
	    
	  // save the first two good jets
	  if (nJets == 1) jet0_ = edm::Ptr<pat::Jet>( h_jets_, 0);
	  if (nJets == 2) jet1_ = edm::Ptr<pat::Jet>( h_jets_, 1);
	    
	  // save all the good jets
	  good_jets_.push_back(edm::Ptr<pat::Jet>(h_jets_, nGoodJets));

	}
	++nGoodJets;

      }

      if ( nJets >= 1 || ignoreCut("One jet or more") ) passCut(ret, "One jet or more");
      else break; 
	
      if ( nJets >= 2 || ignoreCut("Two jets or more") ) passCut(ret, "Two jets or more");
      else break; 
	
      if ( nJets >= 3 || ignoreCut("Three jets or more") ) passCut(ret, "Three jets or more");
      else break; 
	
      if ( nJets >= cut("Min jet multiplicity",int()) || ignoreCut("Min jet multiplicity") ) passCut(ret, "Min jet multiplicity");
      else break; 
	
      passCut(ret, "Jet cuts");

    } // end of jet cuts


      //
      //_____ MET cuts __________________________________
      //      
    if ( considerCut("MET cuts") ) {
      event.getByLabel( metSrc_, h_met_ );
      event.getByLabel( metTCSrc_, h_tcMet_ );
      event.getByLabel( metPFSrc_, h_pfMet_ );
      event.getByLabel( metMuonSrc_, h_caloMet_ );
      event.getByLabel( metTypeIISrc_, h_caloMet2_ );
      
      // default MET
      if ( signed (h_met_->size()) >= 0 ) {
	pat::MET const & met = h_met_->at(0);
	retMet.set(false);
	bool pass = (*metSel_)( met, retMet );
	if (pass){
	  met_ = edm::Ptr<pat::MET>( h_met_, 0);
	}
	else break;
      }
      else break;

      // tcMet
      if ( signed (h_tcMet_->size()) >= 0 ) {
	tcMet_ = edm::Ptr<pat::MET>( h_tcMet_, 0);
      }
      else break;

      // pfMet
	if ( signed (h_pfMet_->size()) >= 0 ) {
	pfMet_ = edm::Ptr<pat::MET>( h_pfMet_, 0);
      }
      else break;

      // caloMet
	if ( signed (h_caloMet_->size()) >= 0 ) {
	caloMet_ = edm::Ptr<reco::CaloMET>( h_caloMet_, 0);
      }
      else break;

      // caloMet corr type 2
	if ( signed (h_caloMet2_->size()) >= 0 ) {
	caloMet2_ = edm::Ptr<reco::CaloMET>( h_caloMet2_, 0);
      }
      else break;

      passCut(ret, "MET cuts");

    } // end of MET cuts


      //
      //_____ Z candidate cuts __________________________________
      //      
    if ( considerCut("Z candidate cuts") ) {

      // opposite muon charges
      if ( nMuons + nLooseMuons > 1 && 
	   ( muon0_ -> charge() + muon1_ -> charge() ) == 0 ) passCut(ret, "Opposite charges");
      else break;

      // dimuon mass
      double _m2 = (muon0_->p4() + muon1_->p4()).M2();
      if ( _m2 > 0 &&
	   sqrt(_m2) > cut("Min dimuon mass", double()) ){

	passCut(ret, "Min dimuon mass");
	_dimuon_mass = sqrt(_m2);

      }
      else break;

      //
      //_____ trigger muon matching
      //
      edm::InputTag _triggerEventSrc("patTriggerEvent");
      event.getByLabel( _triggerEventSrc, h_triggerEvent_ );

      // get the trigger path (like HLT_Mu9)
      const pat::TriggerPath * _path = h_triggerEvent_ -> path( "HLT_Mu9" );
      if (!_path) break;

      // get all trigger filters
      const std::vector<pat::TriggerFilter> * filters = h_triggerEvent_ -> filters();

      // get filter indices used in the path
      std::vector<unsigned> iFilters = _path -> filterIndices();

      // loop over filters in the path
      if (muonMatchedToHlt_firstpass == true){
	std::cout << "Filters used in the trigger path: " << std::endl;
      }

      std::vector<unsigned> keys;

      for (std::vector<unsigned>::const_iterator iFilter = iFilters.begin();
	   iFilter != iFilters.end();
	   ++iFilter){
	const pat::TriggerFilter & filter = (*filters)[*iFilter];

	if (muonMatchedToHlt_firstpass == true){
	  std::cout << "Filter: " << filter.label() << std::endl;
	}
	std::vector<unsigned> moreKeys = filter.objectKeys();
	keys . insert( keys.end(), moreKeys.begin(), moreKeys.end() );
      }
      muonMatchedToHlt_firstpass = false;

      //std::cout << keys.size() << std::endl;

      // all trigger objects in the event
      const std::vector<pat::TriggerObject> * objects = h_triggerEvent_->objects();

      // loop over all keys and match to corresponding objects
      double _dR = 0.2;
      double _dPtRel = 1.0;
      unsigned _nMatched = 0;
      //std::cout << std::endl << "DEBUG: list of trigger objects and the muons" << std::endl;
      for (std::vector<unsigned>::const_iterator key = keys.begin();
	   key != keys.end();
	   ++key){
	//std::cout  << "object pt = " << (*objects)[*key].pt() << std::endl;
	//std::cout  << "muon0 pt = " << muon0_->pt() << std::endl;
	//std::cout  << "muon1 pt = " << muon1_->pt() << std::endl;
	if ( areMatched(*muon0_, (*objects)[*key], _dR, _dPtRel) ) ++_nMatched;
	else if ( areMatched(*muon1_, (*objects)[*key], _dR, _dPtRel) ) ++_nMatched;
      }

      //
      //_____ end of trigger muon matching
	
      if( _nMatched > 0 || ignoreCut("Muon matched to HLT") ) passCut(ret, "Muon matched to HLT");
      else break;

      passCut(ret, "Z candidate cuts"); // Z cuts total

    } // end of Z candidate cuts

    break;
  } // end of while loop

    
    
  return (bool)ret;

  setIgnored(ret);
    
  return false;

}// end of operator()


#endif
