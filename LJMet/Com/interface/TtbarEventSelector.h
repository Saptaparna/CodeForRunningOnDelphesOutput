// -*- C++ -*-
//
// FWLite PAT analyzer-selector for light stop analysis
//
// Gena Kukartsev, March 2012
//
#ifndef LJMet_Com_interface_TtbarEventSelector_h
#define LJMet_Com_interface_TtbarEventSelector_h



#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "PhysicsTools/SelectorUtils/interface/PFElectronSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFMuonSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"



using namespace std;
using trigger::TriggerObject;




class TtbarEventSelector : public BaseEventSelector {

 public:

  
  TtbarEventSelector( std::map<std::string, edm::ParameterSet const> par);  
  
  
  ~TtbarEventSelector();
  
  
  void initialize( void );
  
  
  // main method where the cuts are applied
  bool operator()( edm::EventBase const & event, pat::strbitset & ret);
  

  virtual void ProcessEvent( edm::EventBase const & event,
			     LjmetEventContent & ec );


  boost::shared_ptr<PFJetIDSelectionFunctor> const & jetSel()        const { return jetSel_;}
  boost::shared_ptr<PFMuonSelector>          const & muonSel()       const { return muonSel_;}
  boost::shared_ptr<PFMuonSelector>          const & looseMuonSel()  const { return looseMuonSel_;}
  boost::shared_ptr<PFElectronSelector>      const & pfElectronSel() const { return pfElectronSel_;}
  boost::shared_ptr<PVSelector>              const & pvSel()         const { return pvSel_;}




 protected:

  std::string legend;
  bool bFirstEntry;

  // containers for config parameter values
  std::map<std::string,bool>           mbPar;
  std::map<std::string,int>            miPar;
  std::map<std::string,double>         mdPar;
  std::map<std::string,std::string>    msPar;
  std::map<std::string, edm::InputTag> mtPar;

  boost::shared_ptr<PFMuonSelector>  muonSel_;
  boost::shared_ptr<PFMuonSelector>  looseMuonSel_;
  boost::shared_ptr<PFElectronSelector>      pfElectronSel_;
  boost::shared_ptr<PFJetIDSelectionFunctor> jetSel_;
  boost::shared_ptr<PVSelector>            pvSel_;

  edm::Handle<edm::TriggerResults >           mhEdmTriggerResults;
  edm::Handle<vector<pat::Jet> >              mhJets;
  edm::Handle<vector<pat::Muon> >             mhMuons;
  edm::Handle<vector<pat::Electron> >         mhElectrons;
  edm::Handle<vector<pat::MET> >              mhMet;
  edm::Handle<double>                         h_rho;
  edm::Handle<std::vector<reco::Vertex> >     h_primVtx;

  std::vector<edm::Ptr<reco::Vertex> >  good_pvs_;

  edm::Ptr<pat::Muon>     muon0_;
  edm::Ptr<pat::Muon>     muon1_;
  edm::Ptr<pat::Electron> electron0_;
  edm::Ptr<pat::Electron> electron1_;
};




TtbarEventSelector::TtbarEventSelector( std::map<std::string, edm::ParameterSet const> par) :
  BaseEventSelector( par ){
  
  
  // legend for identifying output messages
  legend = "[";
  legend.append("TtbarEventSelector");
  legend.append("]: ");

  
  std::string _key;


  _key = "pfMuonSelector";
  if ( par.find(_key)!=par.end() ){
    muonSel_ = boost::shared_ptr<PFMuonSelector>( new PFMuonSelector(par[_key]) );
  }
  else {
    std::cout << legend << "muon selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  _key = "looseMuonSelector";
  if ( par.find(_key)!=par.end() ){
    looseMuonSel_ = boost::shared_ptr<PFMuonSelector>( new PFMuonSelector(par[_key]) );
  }
  else {
    std::cout << legend << "loose muon selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  _key = "pfElectronSelector";
  if ( par.find(_key)!=par.end() ){
    pfElectronSel_ = boost::shared_ptr<PFElectronSelector>( new PFElectronSelector(par[_key]) );
  }
  else {
    std::cout << legend << "electron selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  _key = "pfJetIDSelector";
  if ( par.find(_key)!=par.end() ){
    jetSel_ = boost::shared_ptr<PFJetIDSelectionFunctor>( new PFJetIDSelectionFunctor(par[_key]) );
  }
  else {
    std::cout << legend << "jet ID selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  
  _key = "pvSelector";
  if ( par.find(_key)!=par.end() ){
    pvSel_ = boost::shared_ptr<PVSelector>( new PVSelector(par[_key]) );
  }
  else {
    std::cout << legend << "PV selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  
  _key = "event_selector";
  if ( par.find(_key)!=par.end() ){
    mbPar["trigger_cut"]              = par[_key].getParameter<bool>         ("trigger_cut");
    mbPar["dump_trigger"]             = par[_key].getParameter<bool>         ("dump_trigger");
    msPar["trigger_path"]             = par[_key].getParameter<std::string>  ("trigger_path");
    mbPar["pv_cut"]                   = par[_key].getParameter<bool>         ("pv_cut");
    mbPar["hbhe_cut"]                 = par[_key].getParameter<bool>         ("hbhe_cut");
    mbPar["jet_cuts"]                 = par[_key].getParameter<bool>         ("jet_cuts");
    mdPar["jet_minpt"]                = par[_key].getParameter<double>       ("jet_minpt");
    mdPar["jet_maxeta"]               = par[_key].getParameter<double>       ("jet_maxeta");
    miPar["min_jet"]                  = par[_key].getParameter<int>          ("min_jet");
    miPar["max_jet"]                  = par[_key].getParameter<int>          ("max_jet");
    mbPar["muon_cuts"]                = par[_key].getParameter<bool>         ("muon_cuts");
    miPar["min_tight_muon"]           = par[_key].getParameter<int>          ("min_tight_muon");
    mdPar["tight_muon_minpt"]         = par[_key].getParameter<double>       ("tight_muon_minpt");
    mdPar["tight_muon_maxeta"]        = par[_key].getParameter<double>       ("tight_muon_maxeta");
    mdPar["tight_muon_mindeltaR_jet"] = par[_key].getParameter<double>       ("tight_muon_mindeltaR_jet");
    miPar["max_tight_muon"]           = par[_key].getParameter<int>          ("max_tight_muon");
    mdPar["loose_muon_minpt"]         = par[_key].getParameter<double>       ("loose_muon_minpt");
    mdPar["loose_muon_maxeta"]        = par[_key].getParameter<double>       ("loose_muon_maxeta");
    mbPar["loose_muon_veto"]          = par[_key].getParameter<bool>         ("loose_muon_veto");
    mdPar["electron_minEt"]           = par[_key].getParameter<double>       ("electron_minEt");
    mdPar["electron_maxeta"]          = par[_key].getParameter<double>       ("electron_maxeta");
    mbPar["electron_veto"]            = par[_key].getParameter<bool>         ("electron_veto");
    mbPar["met_cuts"]                 = par[_key].getParameter<bool>         ("met_cuts");
    mbPar["btag_cuts"]                = par[_key].getParameter<bool>         ("btag_cuts");
    msPar["btagger"]                  = par[_key].getParameter<std::string>  ("btagger");
    mdPar["btag_min_discr"]           = par[_key].getParameter<double>       ("btag_min_discr");
    mbPar["btag_1"]                   = par[_key].getParameter<bool>         ("btag_1");
    mbPar["btag_2"]                   = par[_key].getParameter<bool>         ("btag_2");
    mbPar["btag_3"]                   = par[_key].getParameter<bool>         ("btag_3");

    mtPar["trigger_collection"]       = par[_key].getParameter<edm::InputTag>("trigger_collection");
    mtPar["pv_collection"]            = par[_key].getParameter<edm::InputTag>("pv_collection");
    mtPar["jet_collection"]           = par[_key].getParameter<edm::InputTag>("jet_collection");
    mtPar["muon_collection"]          = par[_key].getParameter<edm::InputTag>("muon_collection");
    mtPar["electron_collection"]      = par[_key].getParameter<edm::InputTag>("electron_collection");
    mtPar["met_collection"]           = par[_key].getParameter<edm::InputTag>("met_collection");
  }   
  else {
    std::cout << legend << "event selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  
  
  std::cout << legend << "initializing light stop selection" << std::endl;

  
  initialize();
  
}




TtbarEventSelector::~TtbarEventSelector(){}




void TtbarEventSelector::initialize( void ){

  bFirstEntry = true;
  
  push_back("No selection");
  set("No selection");
  
  push_back("Trigger");
  push_back("Primary vertex");
  push_back("HBHE noise and scraping filter");
  push_back("Min tight muon");
  push_back("Max tight muon");
  push_back("Loose muon veto");
  push_back("Electron veto");
  push_back("One jet or more");
  push_back("Two jets or more");
  push_back("Three jets or more");
  push_back("Min jet multiplicity");
  push_back("Max jet multiplicity");
  push_back("Min MET");
  push_back("1 btag or more");
  push_back("2 btag or more");
  push_back("3 btag or more");
  push_back("All cuts");          // sanity check


  
  // TOP PAG sync selection v3

  set("Trigger", mbPar["trigger_cut"]);
  set("Primary vertex", mbPar["pv_cut"]);
  set("HBHE noise and scraping filter", mbPar["hbhe_cut"]);

  if (mbPar["muon_cuts"]){
    set("Min tight muon", miPar["min_tight_muon"]);
    set("Max tight muon", miPar["max_tight_muon"]);
    set("Loose muon veto", mbPar["loose_muon_veto"]);
  }
  else{
    set("Min tight muon", false);
    set("Max tight muon", false);
    set("Loose muon veto", false);
  }

  set("Electron veto", mbPar["electron_veto"]);
      
  if (mbPar["jet_cuts"]){
    set("One jet or more", true);
    set("Two jets or more", true);
    set("Three jets or more", true);
    set("Min jet multiplicity", miPar["min_jet"]);
    set("Max jet multiplicity", miPar["max_jet"]);
  }
  else{
    set("One jet or more", false);
    set("Two jets or more", false);
    set("Three jets or more", false);
    set("Min jet multiplicity", false);
    set("Max jet multiplicity", false);
  }

  if (mbPar["met_cuts"]) set("Min MET", mdPar["min_met"]);

  if (mbPar["btag_cuts"]){
    set("1 btag or more", mbPar["btag_1"]);
    set("2 btag or more", mbPar["btag_2"]);
    set("3 btag or more", mbPar["btag_3"]);
  }
  else{
    set("1 btag or more", false);
    set("2 btag or more", false);
    set("3 btag or more", false);
  }

  set("All cuts", true);


} // initialize() 




bool TtbarEventSelector::operator()( edm::EventBase const & event, pat::strbitset & ret){
  
  pat::strbitset retMuon      = muonSel_->getBitTemplate();
  pat::strbitset retLooseMuon = looseMuonSel_->getBitTemplate();
  pat::strbitset retElectron  = pfElectronSel_->getBitTemplate();
  pat::strbitset retJet       = jetSel_->getBitTemplate();
  
  
  while(1){ // standard infinite while loop trick to avoid nested ifs
    
    passCut(ret, "No selection");
    
    //
    //_____ Trigger cuts __________________________________
    //
    if ( considerCut("Trigger") ) {


      event.getByLabel( mtPar["trigger_collection"], mhEdmTriggerResults );
      //const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
      const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);


      bool passTrig = false;
      unsigned int _tIndex = trigNames.triggerIndex(msPar["trigger_path"]);
      unsigned int _tSize = mhEdmTriggerResults->size();


      // dump trigger names
      if (bFirstEntry && mbPar["dump_trigger"]){
	for (unsigned int i=0; i<_tSize; i++){
	  std::string trigName = trigNames.triggerName(i);
	  
	  cout << i << "   " << trigName;
	  bool fired = mhEdmTriggerResults->accept(trigNames.triggerIndex(trigName));
	  cout <<", FIRED = "<<fired<<endl;
	} 
	std::cout << std::endl << legend << "requested trigger path ["
		  << msPar["trigger_path"] << "]" << std::endl;
	std::cout << legend << "path index: " << _tIndex
		  << " out of [0, "
		  << _tSize-1 << "]" << std::endl << std::endl;
      }


      if ( _tIndex<_tSize){
	passTrig = mhEdmTriggerResults->accept(_tIndex);
      }


      if ( ignoreCut("Trigger") || passTrig ) passCut(ret, "Trigger");
      else break;

    } // end of trigger cuts
    
    

    //
    //_____ Primary vertex cuts __________________________________
    //
    if ( considerCut("Primary vertex") ) {

      if ( (*pvSel_)(event) ){
	passCut(ret, "Primary vertex"); // PV cuts total
      }

    } // end of PV cuts



    
    //
    //_____ HBHE noise and scraping filter________________________
    //
    if ( considerCut("HBHE noise and scraping filter") ) {

      passCut(ret, "HBHE noise and scraping filter"); // PV cuts total

    } // end of PV cuts



    
    //======================================================
    //
    // jet loop
    //
    // We have to do it here, before the muon cuts, because
    // the muon cuts require dR(mu, jet)
    //


    int _n_good_jets = 0;
    int _n_jets = 0;
    //int njetsPF = 0;
    
    if ( miPar["min_tight_muon"] > 0 ) {

      mvSelJets.clear();
      mvAllJets.clear();

      event.getByLabel( mtPar["jet_collection"], mhJets );
      for (vector<pat::Jet>::const_iterator _ijet = mhJets->begin();
	   _ijet != mhJets->end(); ++_ijet){
	
	retJet.set(false);
	bool _pass = false;

	// jet cuts
	while(1){

	  // quality cuts
	  if ( (*jetSel_)( *_ijet, retJet ) ){ }
	  else break; // fail

	  if ( _ijet->pt()>mdPar["jet_minpt"] ){ }
	  else break; // fail

	  if ( fabs(_ijet->eta())<mdPar["jet_maxeta"] ){ }
	  else break;

	  _pass = true;
	  break;
	}
	
	if ( _pass ){
	  // save all the good jets
	  ++_n_good_jets;
	  mvSelJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 

	  //const std::vector<std::pair<std::string, float> > & vpd = _ijet->getPairDiscri();
	  //for ( std::vector<std::pair<std::string, float> >::const_iterator i = vpd.begin();
	  //	i != vpd.end(); ++i){
	  //  std::cout << legend << i->first << ", " << i->second << std::endl;
	  //}
	  
	}
	
	//Fill all the jets regardless if they pass or not
	mvAllJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 
	
	++_n_jets;
	
      } // end of loop over jets
      
    } // end of if (min_tight_muon)


    // get rho correction
    /*
    event.getByLabel( rhoCorrectionSrc_, h_rho); 
    double rho = (*h_rho);
    double rho = 0.0;
    */

    //
    //_____ Muon cuts ________________________________
    //      
    if ( mbPar["muon_cuts"] ) {

      // now, finally, get muons
      event.getByLabel( mtPar["muon_collection"], mhMuons );
      
      // loop over muons
      int _n_muons = 0;
      int _n_tight_muons = 0;
      int _n_loose_muons = 0;

      mvSelMuons.clear();
      mvLooseMuons.clear();
	
      for ( vector<pat::Muon>::const_iterator _imu = mhMuons->begin();
	    _imu != mhMuons->end(); _imu++){
	
	retMuon.set(false);
	bool pass = false;

	//muon cuts
	while(1){

	  // quality cuts
	  // safety for nonexistent track - selector was not safe
	  if ( _imu->globalTrack().isNonnull() && 
	       _imu->globalTrack().isAvailable() ){ }
	  else break; // fail

	  if ( (*muonSel_)( *_imu, retMuon ) ){ }
	  else break; // fail

	  if (_imu->pt()>mdPar["tight_muon_minpt"]){ }
	  else break;
	  
	  if ( fabs(_imu->eta())<mdPar["tight_muon_maxeta"] ){ }
	  else break;

	  // deltaR(mu,jet) cut
	  double min_deltar = numeric_limits<double>::max();
	  for (std::vector<edm::Ptr<pat::Jet> >::const_iterator _ijet = mvSelJets.begin();
	       _ijet!=mvSelJets.end(); ++_ijet){
	    min_deltar = std::min( min_deltar, reco::deltaR(*_imu, *(*_ijet)) );
	  }
	  if ( min_deltar>0.3 ){ }
	  else break;
	  
	  // vertex-pv Z distance
	  reco::Vertex const & pv = pvSel_->vertices()->at(0);
	  if( fabs(_imu->vertex().z() - pv.z()) < 1.0){ }
	  else  break;

	  pass = true; // success
	  break;
	}

	if ( pass ){
	  ++_n_tight_muons;

	  // save the first good muon
	  if (_n_tight_muons == 1){
	    muon0_ = edm::Ptr<pat::Muon>( mhMuons, _n_muons);
	      
	  }

	  // save the second good muon
	  else if (_n_tight_muons == 2){
	    muon1_ = edm::Ptr<pat::Muon>( mhMuons, _n_muons);
	  }

	  // save every good muon
	  mvSelMuons.push_back( edm::Ptr<pat::Muon>( mhMuons, _n_muons) );

	}	  
	else{ // check for loose muon

	  retLooseMuon.set(false);
	  bool pass_loose = false;

	  while(1){

	    // safety for nonexistent track - selector was not safe
	    if ( _imu->globalTrack().isNonnull() && 
		 _imu->globalTrack().isAvailable() ){ }
	    else break;

	    // quality cuts
	    if ( (*looseMuonSel_)( *_imu, retLooseMuon ) ){ }
	    else break;

	    // muon pt cut
	    if (_imu->pt()>mdPar["loose_muon_minpt"]){ }
	    else break;
	    
	    // muon eta cut
	    if ( fabs(_imu->eta())<mdPar["loose_muon_maxeta"] ){ }
	    else break;

	    pass_loose = true;
	    break;
	  }

	  if ( pass_loose ){
	    mvLooseMuons.push_back( edm::Ptr<pat::Muon>( mhMuons, _n_muons) );
	    ++_n_loose_muons;
	  }
	  
	}
	  
	
	++_n_muons;
      } // end of the muon loop

      //std::cout << "DEBUG1: n tight muons=" << _n_tight_muons << std::endl;
      //std::cout << "DEBUG2: n tight muons cut=" << cut("Min tight muon", int()) << std::endl;

      if( _n_tight_muons >= cut("Min tight muon", int()) || ignoreCut("Min tight muon") ) passCut(ret, "Min tight muon");
      else break;
	
      if( _n_tight_muons <= cut("Max tight muon", int()) || ignoreCut("Max tight muon") ) passCut(ret, "Max tight muon");
      else break;
	
      if( (_n_loose_muons == 0)  || ignoreCut("Loose muon veto") ) passCut(ret, "Loose muon veto");
      else break;

    } // end of muon cuts
      
    //
    //_____ Electron cuts __________________________________
    //      
    if ( considerCut("Electron veto") ) {

      event.getByLabel( mtPar["electron_collection"], mhElectrons );

      // loop over electrons
      int n_elec = 0;
      int nSelElectrons = 0;

      mvSelElectrons.clear();

      for ( vector<pat::Electron>::const_iterator _iel = mhElectrons->begin();
	    _iel != mhElectrons->end(); _iel++){

	retElectron.set(false);
	bool pass = false;
	
	while(1){

	  // quality cuts
	  if ( (*pfElectronSel_)(*_iel, retElectron) ){ }
	  else break;

	  // electron Et cut
	  if (_iel->et()>mdPar["electron_minEt"]){ }
	  else break;
	  
	  // electron eta cut
	  if ( fabs(_iel->eta())<mdPar["electron_maxeta"] ){ }
	  else break;

	  pass = true;
	  break;
	}

	if ( pass ){
	  ++nSelElectrons;

	  // save the first good muon
	  if (nSelElectrons == 1){
	    electron0_ = edm::Ptr<pat::Electron>( mhElectrons, n_elec);	      
	  }
	  // save the second good muon
	  else if (nSelElectrons == 2){
	    electron1_ = edm::Ptr<pat::Electron>( mhElectrons, n_elec);	      
	  }
	  // save every good electorn
	  //mvSelElectrons.push_back( edm::Ptr<pat::Electron>( mhElectrons, n_elec) );

	}	  
	else{
	  // save "bad" electrons which are not vetoed
	  mvSelElectrons.push_back( edm::Ptr<pat::Electron>( mhElectrons, n_elec) );
	}
	  
	++n_elec;
      } // end of the electron loop
      
      if( nSelElectrons == 0 || ignoreCut("Electron veto")) passCut(ret, "Electron veto");
      else break;

    } // end of electron cuts




      //
    if ( mbPar["jet_cuts"] ) {

      if ( ignoreCut("One jet or more") || _n_good_jets >= 1 ) passCut(ret, "One jet or more");
      else break; 
	
      if ( ignoreCut("Two jets or more") || _n_good_jets >= 2 ) passCut(ret, "Two jets or more");
      else break; 
	
      if ( ignoreCut("Three jets or more") || _n_good_jets >= 3 ) passCut(ret, "Three jets or more");
      else break; 
	
      if ( ignoreCut("Min jet multiplicity") || _n_good_jets >= cut("Min jet multiplicity",int()) ) passCut(ret, "Min jet multiplicity");
      else break; 
	
      if ( ignoreCut("Max jet multiplicity") || _n_good_jets <= cut("Max jet multiplicity",int()) ) passCut(ret, "Max jet multiplicity");
      else break; 

    } // end of jet cuts
    
    
    
    
    //
    //_____ MET cuts __________________________________
    //      
    event.getByLabel( mtPar["met_collection"], mhMet );
    mpMet = edm::Ptr<pat::MET>( mhMet, 0);
    if ( mbPar["met_cuts"] ) {


      // pfMet
      if ( mpMet.isNonnull() && mpMet.isAvailable() ) {

	pat::MET const & met = mhMet->at(0);
	
	while(1){
	  
	  if ( ignoreCut("Min MET") || met.et()>cut("Min MET", double()) ) passCut(ret, "Min MET");
	  else break;

	  break;
	}
      }
    } // end of MET cuts


    
    //
    //_____ Btagging cuts _____________________
    //
    if ( mbPar["btag_cuts"] ) {
      
      mvSelBtagJets.clear();
      int nBtagJets = 0;
      double bTagCut = mdPar["btag_min_discr"];
      for (std::vector<edm::Ptr<pat::Jet> >::const_iterator _ijet = mvSelJets.begin(); _ijet != mvSelJets.end(); ++_ijet){
	
	if ( (*_ijet)->bDiscriminator( msPar["btagger"] ) > bTagCut ){
	  ++nBtagJets;
	  mvSelBtagJets.push_back(*_ijet); 
	}

	// DEBUG: list all available taggers with values
	//const std::vector<std::pair<std::string, float> > & vpd = (*_ijet)->getPairDiscri();
	//for ( std::vector<std::pair<std::string, float> >::const_iterator i = vpd.begin();
	//      i != vpd.end(); ++i){
	//  std::cout << legend << i->first << ", " << i->second << std::endl;
	//}

      }
			           
      if ( nBtagJets >= 1 || ignoreCut("1 btag or more") )  passCut(ret, "1 btag or more");
      else break;
      if ( nBtagJets >= 2 || ignoreCut("2 btag or more") )  passCut(ret, "2 btag or more");
      else break;
      if ( nBtagJets >= 3 || ignoreCut("3 btag or more") )  passCut(ret, "3 btag or more");
      else break;

    }
    

    passCut(ret, "All cuts");
    break;

  } // end of while loop


  bFirstEntry = false;
    

  return (bool)ret;
  
  setIgnored(ret);

  return false;
}// end of operator()




void TtbarEventSelector::ProcessEvent( edm::EventBase const & event,
				       LjmetEventContent & ec ){
  //
  // Compute analysis-specific quantities in the event, return via reference as maps
  //

  //m_double_branch["pi"] = 3.14;

  //m_double_branch["ele_et"] = -1.0;
  //if (mvSelElectrons.size()>0) m_double_branch["ele_et"] = mvSelElectrons[0]->et();

  

  return;
}


#endif
