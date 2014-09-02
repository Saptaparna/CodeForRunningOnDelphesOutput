/*
  Calculator for the SLiTT analysis
  
  Author: Gena Kukartsev, 2012
  Modified by Saptaparna Bhattacharya for CMS analysis with 8 TeV data and for use with Delphes. Saptaparna has modified the code to include various electron and muon ID requirements, has added full generator level history and also modified the code to include jet substructure variables. 

*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "LJMet/Com/interface/TopElectronSelector.h"

#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"

#include "TClonesArray.h"
#include "LJMet/DelphesData/interface/DelphesClasses.h"

using std::cout;
using std::endl;

class LjmetFactory;

class DileptonCalc : public BaseCalc{
  
public:
  
  DileptonCalc();
  virtual ~DileptonCalc(){}
  
  virtual int BeginJob();
  virtual int AnalyzeEvent(ExRootTreeReader * reader, Int_t entry,
			   BaseEventSelector * selector);
  virtual int EndJob(){return 0;}

private:
  
  bool                      isMc;
  std::string               dataType;
  edm::InputTag             rhoSrc_it;
  edm::InputTag             pvCollection_it;
  edm::InputTag             genParticles_it;
  std::vector<unsigned int> keepPDGID;
  std::vector<unsigned int> keepMomPDGID;
  bool keepFullMChistory;
  
  double rhoIso;

  boost::shared_ptr<TopElectronSelector>     electronSelL_, electronSelM_, electronSelT_;
  std::vector<reco::Vertex> goodPVs;
  int findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi);
  double mdeltaR(double eta1, double phi1, double eta2, double phi2);
  void fillMotherInfo(TClonesArray * genparticles, Int_t mInd,
		      int i, vector <int> & momid, vector <int> & momstatus,
		      vector<double> & mompt, vector<double> & mometa,
		      vector<double> & momphi, vector<double> & momenergy);
};


static int reg = LjmetFactory::GetInstance()->Register(new DileptonCalc(), "DileptonCalc");


DileptonCalc::DileptonCalc(){
}

int DileptonCalc::BeginJob(){
  
  if (mPset.exists("dataType"))     dataType = mPset.getParameter<std::string>("dataType");
  else                              dataType = "None";
  
  if (mPset.exists("rhoSrc"))       rhoSrc_it = mPset.getParameter<edm::InputTag>("rhoSrc");
  else                              rhoSrc_it = edm::InputTag("kt6PFJetsForIsolation", "rho", "PAT");
  
  if (mPset.exists("pvCollection")) pvCollection_it = mPset.getParameter<edm::InputTag>("pvCollection");
  else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (mPset.exists("isMc"))         isMc = mPset.getParameter<bool>("isMc");
  else                              isMc = false;
  
  if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
  else                              genParticles_it = edm::InputTag("prunedGenParticles");
  
  if (mPset.exists("keepPDGID"))    keepPDGID = mPset.getParameter<std::vector<unsigned int> >("keepPDGID");
  else                              keepPDGID.clear();
  
  if (mPset.exists("keepMomPDGID")) keepMomPDGID = mPset.getParameter<std::vector<unsigned int> >("keepMomPDGID");
  else                              keepMomPDGID.clear();

  if (mPset.exists("keepFullMChistory")) keepFullMChistory = mPset.getParameter<bool>("keepFullMChistory");
  else                                   keepFullMChistory = false;
  cout << "keepFullMChistory "     <<    keepFullMChistory << endl;
  
  if ( mPset.exists("cutbasedIDSelectorLoose")){
    electronSelL_ = boost::shared_ptr<TopElectronSelector>( 
	new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorLoose")) );
  }
  else {
    std::cout << "DileptonCalc: Loose electron selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  if ( mPset.exists("cutbasedIDSelectorMedium")){
    electronSelM_ = boost::shared_ptr<TopElectronSelector>( 
	new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorMedium")) );
  }
  else {
    std::cout << "DileptonCalc: Medium electron selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }
  if ( mPset.exists("cutbasedIDSelectorTight")){
    electronSelT_ = boost::shared_ptr<TopElectronSelector>( 
	new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorTight")) );
  }
  else {
    std::cout << "DileptonCalc: Tight electron selector not configured, exiting"
	      << std::endl;
    std::exit(-1);
  }


  return 0;
}

int DileptonCalc::AnalyzeEvent(ExRootTreeReader * reader, Int_t entry,
			       BaseEventSelector * selector){
  //
  // compute event variables here
  //



  // Get pointers to branches used in this analysis
  TClonesArray *b_Event = reader->UseBranch("Event");
  TClonesArray *b_Jet = reader->UseBranch("Jet");
  TClonesArray *b_Electron = reader->UseBranch("Electron");
  TClonesArray *b_Muon = reader->UseBranch("Muon");
  TClonesArray *b_MET = reader->UseBranch("MissingET");
  TClonesArray *b_Truth = reader->UseBranch("Particle");
  TClonesArray *b_CAJet = reader->UseBranch("CAJet");

  // Load selected branches with data from specified event
  reader->ReadEntry(entry);
  //reader->LoadTree(entry);
 
  //
  // _____ Get objects from the selector _____________________
  //
  std::vector<int>     const & vSelMuons     = selector->GetSelectedMuons();
  std::vector<int>     const & vSelElectrons = selector->GetSelectedElectrons();
  std::vector<int>     const & vSelJets      = selector->GetSelectedJets();
  std::vector<int>     const & vSelCAJets      = selector->GetSelectedCAJets();
  MissingET *          const pMet            = selector->GetMet();
  std::vector<unsigned int>             const & vSelTriggers  = selector->GetSelectedTriggers();
  
  //
  // _____ Primary dataset (from python cfg) _____________________
  //
  //
  int dataEE = 0;
  int dataEM = 0;
  int dataMM = 0;
  
  if      (dataType == "EE" or dataType == "ElEl") dataEE = 1; 
  else if (dataType == "EM" or dataType == "ElMu") dataEM = 1;
  else if (dataType == "MM" or dataType == "MuMu") dataMM = 1;
  else if (dataType == "All" or dataType == "ALL") {
    dataEE = 1; dataEM = 1; dataMM = 1;
  }

  SetValue("dataEE", dataEE);
  SetValue("dataEM", dataEM);
  SetValue("dataMM", dataMM);
  
  //
  // ____ Trigger ____________________________
  //
  int passEE = 0;
  int passEM = 0;
  int passMM = 0;
  
  if (vSelTriggers.size() == 3){    
    passEE = (int)vSelTriggers.at(0);
    passEM = (int)vSelTriggers.at(1);
    passMM = (int)vSelTriggers.at(2);
  }
  
  SetValue("trigEE", passEE);
  SetValue("trigEM", passEM);
  SetValue("trigMM", passMM);
  
  //-------EventWeight for HT binned samples------

  double eventWeight = 0.0;
  LHEFEvent* event=(LHEFEvent*)b_Event->At(0);
  eventWeight=event->Weight;
  //cout << "eventWeight = " << eventWeight << endl;

  SetValue("eventWeight", eventWeight);

  //
  //_____ Event kinematics __________________
  //
  
  //Primary vertices
  int _npv = 1; // no info in Delphes
  SetValue("nPV", _npv);
  
  //
  //_____ Electrons _________________________
  //
  
  //Four vector
  std::vector <double> elPt;
  std::vector <double> elEta;
  std::vector <double> elPhi;
  std::vector <double> elEnergy;
  
  //Quality criteria
  std::vector <double> elRelIso;
  std::vector <double> elDxy;
  std::vector <int>    elNotConversion;
  std::vector <int>    elChargeConsistent;
  std::vector <int>    elIsEBEE; 
  std::vector <int>    elQuality;
  std::vector <int>    elCharge;
  
  //ID requirement
  std::vector <double> elDeta;
  std::vector <double> elDphi;
  std::vector <double> elSihih;
  std::vector <double> elHoE;
  std::vector <double> elD0;
  std::vector <double> elDZ;
  std::vector <double> elOoemoop;
  std::vector <int>    elMHits;
  std::vector <int>    elVtxFitConv;

  //Extra info about isolation
  std::vector <double> elChIso;
  std::vector <double> elNhIso;
  std::vector <double> elPhIso;
  std::vector <double> elAEff;
  std::vector <double> elRhoIso;

  //mother-information
  //Generator level information -- MC matching
  vector<double> elGen_Reco_dr;
  vector<int> elPdgId;
  vector<int> elStatus;
  vector<int> elMatched;
  vector<int> elNumberOfMothers;
  vector<double> elMother_pt;
  vector<double> elMother_eta;
  vector<double> elMother_phi;
  vector<double> elMother_energy;
  vector<int> elMother_id;
  vector<int> elMother_status;
  //Matched gen electron information:
  vector<double> elMatchedPt;
  vector<double> elMatchedEta;
  vector<double> elMatchedPhi;
  vector<double> elMatchedEnergy;
  vector<double> elIsolation; 
  //edm::Handle<double> rhoHandle;
  //event.getByLabel(rhoSrc_it, rhoHandle);
  //rhoIso = std::max(*(rhoHandle.product()), 0.0);
  rhoIso = 0.0; // dummy

  //pat::strbitset retElectron  = electronSelL_->getBitTemplate();
  bool retElectronT,retElectronM,retElectronL;



  //
  //_____Electrons______
  //
  
  for (std::vector<int>::const_iterator iel = vSelElectrons.begin();
       iel != vSelElectrons.end(); ++iel){   

    Electron * _el = (Electron *) b_Electron->At(*iel);
    
    //Four vector
    elPt     . push_back(_el->PT); //Must check: why ecalDrivenMomentum?
    elEta    . push_back(_el->Eta);
    elPhi    . push_back(_el->Phi);
    elEnergy . push_back(_el->P4().Energy());  
    elIsolation.push_back(1.0);   

  
    //Isolation
    double AEff  = 0.0;
    double chIso = 0.0;
    double nhIso = 0.0;
    double phIso = 0.0;
    double relIso = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) ) / _el->PT;
      
    elChIso  . push_back(chIso);
    elNhIso  . push_back(nhIso);
    elPhIso  . push_back(phIso);
    elAEff   . push_back(AEff);
    elRhoIso . push_back(rhoIso);

    elRelIso . push_back(relIso);

    //Conversion rejection
    int nLostHits = 0;
    //double dist   = (*iel)->convDist();
    //double dcot   = (*iel)->convDcot();
    //int notConv   = nLostHits == 0 and (fabs(dist) > 0.02 or fabs(dcot) > 0.02);
    int notConv   = 1;
    elCharge.push_back(_el->Charge); 
    elNotConversion . push_back(notConv);
      
    retElectronL = true; //(*electronSelL_)(**iel, event, retElectron);
    retElectronM = true; //(*electronSelM_)(**iel, event, retElectron);
    retElectronT = true; //(*electronSelT_)(**iel, event, retElectron);
      
    elQuality.push_back((retElectronT<<2) + (retElectronM<<1) + retElectronL);

    //IP: for some reason this is with respect to the first vertex in the collection
    elDxy.push_back( 0.0 );
    elDZ.push_back( 0.0 );
    elChargeConsistent.push_back(1);
    elIsEBEE.push_back(1);
    elDeta.push_back(1);
    elDphi.push_back(1);
    elSihih.push_back(1);
    elHoE.push_back(_el->EhadOverEem);
    elD0.push_back(0.0);
    elOoemoop.push_back(0.01); 
    elMHits.push_back(12);
    elVtxFitConv.push_back(1);
    //std::cout << isMc << " " << keepFullMChistory << std::endl;
    if(isMc && keepFullMChistory){
      cout << "start\n";
      // this does not seem to work but not needed anyway
      GenParticle * _gen = (GenParticle *)_el->Particle.GetObject();
      double closestDR = 10000.;
      if (_gen) {
	closestDR = mdeltaR( _el->Eta, _el->Phi, _gen->Eta, _gen->Phi);
      }
      cout << "closestDR "<<closestDR <<endl;
      if(closestDR < 0.3){
	elGen_Reco_dr.push_back(closestDR);
	elPdgId.push_back(_gen->PID);
	elStatus.push_back(_gen->Status);
	elMatched.push_back(1);
	elMatchedPt.push_back( _gen->PT);
	elMatchedEta.push_back(_gen->Eta);
	elMatchedPhi.push_back(_gen->Phi);
	elMatchedEnergy.push_back(_gen->E);
	int oldSize = elMother_id.size();
	fillMotherInfo(b_Truth, _gen->M1, 0, elMother_id, elMother_status, elMother_pt, elMother_eta, elMother_phi, elMother_energy);
	elNumberOfMothers.push_back(elMother_id.size()-oldSize);
      }
      if(closestDR >= 0.3){
	elNumberOfMothers.push_back(-1);
	elGen_Reco_dr.push_back(-1.0);
	elPdgId.push_back(-1);
	elStatus.push_back(-1);
	elMatched.push_back(0);
	elMatchedPt.push_back(-1000.0);
	elMatchedEta.push_back(-1000.0);
	elMatchedPhi.push_back(-1000.0);
	elMatchedEnergy.push_back(-1000.0);
	
      }
      
      
    }//closing the isMC checking criteria       
    /*
    */
  }
  
  //Four vector
  SetValue("elPt"     , elPt);
  SetValue("elEta"    , elEta);
  SetValue("elPhi"    , elPhi);
  SetValue("elEnergy" , elEnergy);
  
  SetValue("elCharge", elCharge);
  
  SetValue("elIsolation", elIsolation);

  //Quality requirements
  SetValue("elRelIso" , elRelIso); //Isolation
  SetValue("elDxy"    , elDxy);    //Dxy
  SetValue("elNotConversion" , elNotConversion);  //Conversion rejection
  SetValue("elChargeConsistent", elChargeConsistent);
  SetValue("elIsEBEE", elIsEBEE);
  SetValue("elQuality", elQuality);

  //ID cuts 
  SetValue("elDeta", elDeta);
  SetValue("elDphi", elDphi);
  SetValue("elSihih", elSihih);
  SetValue("elHoE", elHoE);
  SetValue("elD0", elD0);
  SetValue("elDZ", elDZ);
  SetValue("elOoemoop", elOoemoop);
  SetValue("elMHits", elMHits);
  SetValue("elVtxFitConv", elVtxFitConv);

  //Extra info about isolation
  SetValue("elChIso" , elChIso);
  SetValue("elNhIso" , elNhIso);
  SetValue("elPhIso" , elPhIso);
  SetValue("elAEff"  , elAEff);
  SetValue("elRhoIso", elRhoIso);

  //MC matching -- mother information
  SetValue("elNumberOfMothers", elNumberOfMothers);
  SetValue("elGen_Reco_dr", elGen_Reco_dr);
  SetValue("elPdgId", elPdgId);
  SetValue("elStatus", elStatus);
  SetValue("elMatched",elMatched);
  SetValue("elMother_pt", elMother_pt);
  SetValue("elMother_eta", elMother_eta);
  SetValue("elMother_phi", elMother_phi);
  SetValue("elMother_energy", elMother_energy);
  SetValue("elMother_status", elMother_status);
  SetValue("elMother_id", elMother_id);
  //Matched gen muon information:
  SetValue("elMatchedPt", elMatchedPt);
  SetValue("elMatchedEta", elMatchedEta);
  SetValue("elMatchedPhi", elMatchedPhi);
  SetValue("elMatchedEnergy", elMatchedEnergy);
  
  



  //
  //_____ Muons _____________________________
  //
  
  std::vector <int> muCharge;
  std::vector <int> muGlobal;
  
  //Four vector
  std::vector <double> muPt;
  std::vector <double> muEta;
  std::vector <double> muPhi;
  std::vector <double> muEnergy;
  std::vector <double> muIsolation;
  
  //Quality criteria
  std::vector <double> muChi2;
  std::vector <double> muDxy;
  std::vector <double> muDz;
  std::vector <double> muRelIso;
  
  std::vector <int> muNValMuHits;
  std::vector <int> muNMatchedStations;
  std::vector <int> muNValPixelHits;
  std::vector <int> muNTrackerLayers;


  //Extra info about isolation
  std::vector <double> muChIso;
  std::vector <double> muNhIso;
  std::vector <double> muGIso;
  std::vector <double> muPuIso;

  //Generator level information -- MC matching
  vector<double> muGen_Reco_dr;
  vector<int> muPdgId;
  vector<int> muStatus;
  vector<int> muMatched;
  vector<int> muNumberOfMothers;
  vector<double> muMother_pt;
  vector<double> muMother_eta;
  vector<double> muMother_phi;
  vector<double> muMother_energy;
  vector<int> muMother_id;
  vector<int> muMother_status;
  //Matched gen muon information:
  vector<double> muMatchedPt;
  vector<double> muMatchedEta;
  vector<double> muMatchedPhi;
  vector<double> muMatchedEnergy;


      
  for (std::vector<int>::const_iterator imu = vSelMuons.begin();
       imu != vSelMuons.end(); ++imu){   
    
    Muon * _mu = (Muon *) b_Muon->At(*imu);

    //charge
    muCharge.push_back(_mu->Charge);
    
    //Four vector
    muPt     . push_back(_mu->PT);
    muEta    . push_back(_mu->Eta);
    muPhi    . push_back(_mu->Phi);
    muEnergy . push_back(_mu->P4().Energy());  
    muIsolation.push_back(1.0);   
 
    muGlobal.push_back((true<<2)+true);
    //Chi2
    muChi2 . push_back(0.0);
      
    //Isolation
    double chIso  = 0.0;
    double nhIso  = 0.0;
    double gIso   = 0.0;
    double puIso  = 0.0;
    double relIso = 0.0;
    muRelIso . push_back(relIso);
    
    muChIso . push_back(chIso);
    muNhIso . push_back(nhIso);
    muGIso  . push_back(gIso);
    muPuIso . push_back(puIso);
    
    //IP: for some reason this is with respect to the first vertex in the collection
    muDxy . push_back(0.0);
    muDz  . push_back(0.0);
    //Numbers of hits
    muNValMuHits       . push_back(12);
    muNMatchedStations . push_back(2);
    muNValPixelHits    . push_back(12);
    muNTrackerLayers   . push_back(5);
    
    if(isMc && keepFullMChistory){
      /*      
      edm::Handle<reco::GenParticleCollection> genParticles;
      event.getByLabel(genParticles_it, genParticles);
      int matchId = findMatch(*genParticles, 13, (*imu)->eta(), (*imu)->phi());
      double closestDR = 10000.;
      if (matchId>=0) {
	const reco::GenParticle & p = (*genParticles).at(matchId);
	closestDR = mdeltaR( (*imu)->eta(), (*imu)->phi(), p.eta(), p.phi());
	if(closestDR < 0.3){
	  muGen_Reco_dr.push_back(closestDR);
	  muPdgId.push_back(p.pdgId());
	  muStatus.push_back(p.status());
	  muMatched.push_back(1);
	  muMatchedPt.push_back( p.pt());
	  muMatchedEta.push_back(p.eta());
	  muMatchedPhi.push_back(p.phi());
	  muMatchedEnergy.push_back(p.energy());
	  int oldSize = muMother_id.size();
	  fillMotherInfo(p.mother(), 0, muMother_id, muMother_status, muMother_pt, muMother_eta, muMother_phi, muMother_energy);
	  muNumberOfMothers.push_back(muMother_id.size()-oldSize);
	}
      } 
      if(closestDR >= 0.3){
	muNumberOfMothers.push_back(-1);
	muGen_Reco_dr.push_back(-1.0);
	muPdgId.push_back(-1);
	muStatus.push_back(-1);
	muMatched.push_back(0);
	muMatchedPt.push_back(-1000.0);
	muMatchedEta.push_back(-1000.0);
	muMatchedPhi.push_back(-1000.0);
	muMatchedEnergy.push_back(-1000.0);
	
      }
      */
    }
  }
  
  
  SetValue("muCharge", muCharge);
  SetValue("muGlobal", muGlobal);
  //Four vector
  SetValue("muPt"     , muPt);
  SetValue("muEta"    , muEta);
  SetValue("muPhi"    , muPhi);
  SetValue("muEnergy" , muEnergy);
  SetValue("muIsolation", muIsolation);

  //Quality criteria
  SetValue("muChi2"   , muChi2);
  SetValue("muDxy"    , muDxy);
  SetValue("muDz"     , muDz);
  SetValue("muRelIso" , muRelIso);

  SetValue("muNValMuHits"       , muNValMuHits);
  SetValue("muNMatchedStations" , muNMatchedStations);
  SetValue("muNValPixelHits"    , muNValPixelHits);
  SetValue("muNTrackerLayers"   , muNTrackerLayers);
  
  //Extra info about isolation
  SetValue("muChIso", muChIso);
  SetValue("muNhIso", muNhIso);
  SetValue("muGIso" , muGIso);
  SetValue("muPuIso", muPuIso);

   //MC matching -- mother information
  SetValue("muGen_Reco_dr", muGen_Reco_dr);
  SetValue("muPdgId", muPdgId);
  SetValue("muStatus", muStatus);
  SetValue("muMatched",muMatched);
  SetValue("muMother_pt", muMother_pt);
  SetValue("muMother_eta", muMother_eta);
  SetValue("muMother_phi", muMother_phi);
  SetValue("muMother_energy", muMother_energy);
  SetValue("muMother_status", muMother_status);
  SetValue("muMother_id", muMother_id);
  SetValue("muNumberOfMothers", muNumberOfMothers);
  //Matched gen muon information:
  SetValue("muMatchedPt", muMatchedPt);
  SetValue("muMatchedEta", muMatchedEta);
  SetValue("muMatchedPhi", muMatchedPhi);
  SetValue("muMatchedEnergy", muMatchedEnergy);  




  //
  //_____ Jets ______________________________
  //


  //Four vector
  std::vector <double> CATopJetPt;
  std::vector <double> CATopJetEta;
  std::vector <double> CATopJetPhi;
  std::vector <double> CATopJetEnergy;

  std::vector <double> CATopJetCSV;
  //   std::vector <double> CATopJetRCN;

  //Identity
  std::vector <int> CATopJetIndex;
  std::vector <int> CATopJetnDaughters;

  //Top-like properties
  std::vector <double> CATopJetTopMass;
  std::vector <double> CATopJetMinPairMass;

  //Daughter four vector and index
  std::vector <double> CATopDaughterPt;
  std::vector <double> CATopDaughterEta;
  std::vector <double> CATopDaughterPhi;
  std::vector <double> CATopDaughterEnergy;

  std::vector <int> CATopDaughterMotherIndex;

  /*
  for (std::vector<pat::Jet>::const_iterator ijet = topJets->begin(); ijet != topJets->end(); ijet++){

    int index = (int)(ijet-topJets->begin());

    CATopJetPt     . push_back(ijet->pt());
    CATopJetEta    . push_back(ijet->eta());
    CATopJetPhi    . push_back(ijet->phi());
    CATopJetEnergy . push_back(ijet->energy());
    
    CATopJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
//     CATopJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    
			       
    CATopJetIndex      . push_back(index);
    CATopJetnDaughters . push_back((int)ijet->numberOfDaughters());

    reco::CATopJetTagInfo* jetInfo = (reco::CATopJetTagInfo*) ijet->tagInfo("CATop");

    CATopJetTopMass     . push_back(jetInfo->properties().topMass);
    CATopJetMinPairMass . push_back(jetInfo->properties().minMass);

    for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
      CATopDaughterPt     . push_back(ijet->daughter(ui)->pt());
      CATopDaughterEta    . push_back(ijet->daughter(ui)->eta());
      CATopDaughterPhi    . push_back(ijet->daughter(ui)->phi());
      CATopDaughterEnergy . push_back(ijet->daughter(ui)->energy());        
      
      CATopDaughterMotherIndex . push_back(index);      
    }
  }
  */

  //Four vector
  SetValue("CATopJetPt"     , CATopJetPt);
  SetValue("CATopJetEta"    , CATopJetEta);
  SetValue("CATopJetPhi"    , CATopJetPhi);
  SetValue("CATopJetEnergy" , CATopJetEnergy);

  SetValue("CATopJetCSV"    , CATopJetCSV);
  //   SetValue("CATopJetRCN"    , CATopJetRCN);

  //Identity
  SetValue("CATopJetIndex"      , CATopJetIndex);
  SetValue("CATopJetnDaughters" , CATopJetnDaughters);

  //Properties
  SetValue("CATopJetTopMass"     , CATopJetTopMass);
  SetValue("CATopJetMinPairMass" , CATopJetMinPairMass);

  //Daughter four vector and index
  SetValue("CATopDaughterPt"     , CATopDaughterPt);
  SetValue("CATopDaughterEta"    , CATopDaughterEta);
  SetValue("CATopDaughterPhi"    , CATopDaughterPhi);
  SetValue("CATopDaughterEnergy" , CATopDaughterEnergy);

  SetValue("CATopDaughterMotherIndex"      , CATopDaughterMotherIndex);


  //Four vector
  std::vector <double> CAWJetPt;
  std::vector <double> CAWJetEta;
  std::vector <double> CAWJetPhi;
  std::vector <double> CAWJetEnergy;

  std::vector <double> CAWJetCSV;
//   std::vector <double> CAWJetRCN;

  //Identity
  std::vector <int> CAWJetIndex;
  std::vector <int> CAWJetnDaughters;

  //Mass
  std::vector <double> CAWJetMass;
  
  //Daughter four vector and index
  std::vector <double> CAWDaughterPt;
  std::vector <double> CAWDaughterEta;
  std::vector <double> CAWDaughterPhi;
  std::vector <double> CAWDaughterEnergy;

  std::vector <int> CAWDaughterMotherIndex;

  /*
  for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++){

    int index = (int)(ijet-CAWJets->begin());

    //Four vector
    CAWJetPt     . push_back(ijet->pt());
    CAWJetEta    . push_back(ijet->eta());
    CAWJetPhi    . push_back(ijet->phi());
    CAWJetEnergy . push_back(ijet->energy());        

    CAWJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
//     CAWJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    

    //Identity
    CAWJetIndex      . push_back(index);
    CAWJetnDaughters . push_back((int)ijet->numberOfDaughters());

    //Mass
    CAWJetMass . push_back(ijet->mass());
    
    for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
      CAWDaughterPt     . push_back(ijet->daughter(ui)->pt());
      CAWDaughterEta    . push_back(ijet->daughter(ui)->eta());
      CAWDaughterPhi    . push_back(ijet->daughter(ui)->phi());
      CAWDaughterEnergy . push_back(ijet->daughter(ui)->energy());        
      
      CAWDaughterMotherIndex . push_back(index);      
    }
  }
  */

  //Four vector
  SetValue("CAWJetPt"     , CAWJetPt);
  SetValue("CAWJetEta"    , CAWJetEta);
  SetValue("CAWJetPhi"    , CAWJetPhi);
  SetValue("CAWJetEnergy" , CAWJetEnergy);

  SetValue("CAWJetCSV"    , CAWJetCSV);
//   SetValue("CAWJetRCN"    , CAWJetRCN);

  //Identity
  SetValue("CAWJetIndex"      , CAWJetIndex);
  SetValue("CAWJetnDaughters" , CAWJetnDaughters);
  
  //Mass
  SetValue("CAWJetMass"     , CAWJetMass);

  //Daughter four vector and index
  SetValue("CAWDaughterPt"     , CAWDaughterPt);
  SetValue("CAWDaughterEta"    , CAWDaughterEta);
  SetValue("CAWDaughterPhi"    , CAWDaughterPhi);
  SetValue("CAWDaughterEnergy" , CAWDaughterEnergy);

  SetValue("CAWDaughterMotherIndex" , CAWDaughterMotherIndex);


  //Get all CA8 jets (not just for W and Top)
  //edm::InputTag CA8JetColl = edm::InputTag("goodPatJetsCA8PF");
  //edm::Handle<std::vector<pat::Jet> > CA8Jets;
  //event.getByLabel(CA8JetColl, CA8Jets);

  //Four vector
  std::vector <double> CA8JetPt;
  std::vector <double> CA8JetEta;
  std::vector <double> CA8JetPhi;
  std::vector <double> CA8JetEnergy;

  std::vector <double> CA8JetCSV;
  //   std::vector <double> CA8JetRCN;

  /*
  for (std::vector<pat::Jet>::const_iterator ijet = CA8Jets->begin(); ijet != CA8Jets->end(); ijet++){
    
  //Four vector
    CA8JetPt     . push_back(ijet->pt());
    CA8JetEta    . push_back(ijet->eta());
    CA8JetPhi    . push_back(ijet->phi());
    CA8JetEnergy . push_back(ijet->energy());

    CA8JetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
    //     CA8JetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    
  }
  */

  //Four vector
  SetValue("CA8JetPt"     , CA8JetPt);
  SetValue("CA8JetEta"    , CA8JetEta);
  SetValue("CA8JetPhi"    , CA8JetPhi);
  SetValue("CA8JetEnergy" , CA8JetEnergy);

  SetValue("CA8JetCSV"    , CA8JetCSV);
  //   SetValue("CA8JetRCN"    , CA8JetRCN);
 
  //Get AK5 Jets
  //Four vector
  std::vector <double> AK5JetPt;
  std::vector <double> AK5JetEta;
  std::vector <double> AK5JetPhi;
  std::vector <double> AK5JetEnergy;

  std::vector <int>    AK5JetTBag;
  std::vector <double> AK5JetRCN;
  std::vector <int>  AK5JetnChHad;
  std::vector <int>  AK5JetnNeuHad; 

  for (std::vector<int>::const_iterator ijet = vSelJets.begin();
       ijet != vSelJets.end(); ++ijet){   

    Jet * _jet = (Jet *) b_Jet->At(*ijet);
    
    //Four vector
    TLorentzVector lv = _jet->P4();

    AK5JetPt     . push_back(lv.Pt());
    AK5JetEta    . push_back(lv.Eta());
    AK5JetPhi    . push_back(lv.Phi());
    AK5JetEnergy . push_back(lv.Energy());
    
    AK5JetTBag   . push_back(_jet->BTag);
    AK5JetRCN    . push_back(10.0);    
    AK5JetnChHad    . push_back(10);  //saving charged hadron multiplicity.
    AK5JetnNeuHad   . push_back(1);  //saving neutral hadron multiplicity.
}

  //Four vector
  SetValue("AK5JetPt"     , AK5JetPt);
  SetValue("AK5JetEta"    , AK5JetEta);
  SetValue("AK5JetPhi"    , AK5JetPhi);
  SetValue("AK5JetEnergy" , AK5JetEnergy);

  SetValue("AK5JetTBag"   , AK5JetTBag);
  SetValue("AK5JetRCN"    , AK5JetRCN);
  SetValue("AK5JetnChHad", AK5JetnChHad);
  SetValue("AK5JetnNeuHad", AK5JetnNeuHad);

  std::vector <double> CAJetPt;
  std::vector <double> CAJetEta;
  std::vector <double> CAJetPhi;
  std::vector <double> CAJetEnergy;

  std::vector <int> CAJetBTag;
  std::vector <int> CAJetWTag;
  std::vector <int> CAJetTopTag;
  
  std::vector <double> CAJetMass;
  std::vector <double> CAJetTrimmedMass;
  std::vector <double> CAJetNsubJets;
  std::vector <double> CAJetTau1;
  std::vector <double> CAJetTau2;
  std::vector <double> CAJetTau3; 
  std::vector <double> CAJetMassDrop;  
 
  for (std::vector<int>::const_iterator icajet = vSelCAJets.begin();
       icajet != vSelCAJets.end(); ++icajet) {

    Jet * _cajet = (Jet *) b_CAJet->At(*icajet);
    
    TLorentzVector lv = _cajet->P4();

    CAJetPt.push_back(lv.Pt());
    CAJetEta.push_back(lv.Eta()); 
    CAJetPhi.push_back(lv.Phi());
    CAJetEnergy.push_back(lv.Energy());
    
    CAJetBTag.push_back(_cajet->BTag);
    CAJetWTag.push_back(_cajet->WTag);
    CAJetTopTag.push_back(_cajet->TopTag);  
    CAJetMass.push_back(_cajet->Mass);
    CAJetTrimmedMass.push_back(_cajet->TrimmedMass);
    CAJetNsubJets.push_back(_cajet->NSubJets);
    CAJetTau1.push_back(_cajet->Tau1);
    CAJetTau2.push_back(_cajet->Tau2);
    CAJetTau3.push_back(_cajet->Tau3);
    CAJetMassDrop.push_back(_cajet->MassDrop);
    

   }

   //Four vector
  SetValue("CAJetPt"     , CAJetPt);
  SetValue("CAJetEta"    , CAJetEta);
  SetValue("CAJetPhi"    , CAJetPhi);
  SetValue("CAJetEnergy" , CAJetEnergy);

  SetValue("CAJetBTag"   , CAJetBTag);
  SetValue("CAJetWTag"    , CAJetWTag);
  SetValue("CAJetTopTag", CAJetTopTag);

  SetValue("CAJetMass"   , CAJetMass);
  SetValue("CAJetTrimmedMass"    , CAJetTrimmedMass);
  SetValue("CAJetNsubJets", CAJetNsubJets);
  SetValue("CAJetTau1", CAJetTau1);
  SetValue("CAJetTau2", CAJetTau2);
  SetValue("CAJetTau3", CAJetTau3);
  SetValue("CAJetMassDrop", CAJetMassDrop);



  // MET
  double _met = -9999.0;
  double _met_phi = -9999.0;
  // Corrected MET
  double _corr_met = -9999.0;
  double _corr_met_phi = -9999.0;
  //double _calo_met = -9999.0;
  //double _calo_met_phi = -9999.0;


  if(b_MET->GetEntries() > 0){

    MissingET * met = (MissingET *) b_MET->At(0);
    _met = met->MET;
    _met_phi = met->Phi;

    _corr_met = _met;
    _corr_met_phi = _met_phi;

  }

  SetValue("met", _met);
  SetValue("met_phi", _met_phi);
  SetValue("corr_met", _corr_met);
  SetValue("corr_met_phi", _corr_met_phi);
  //SetValue("calo_met", _calo_met);
  //SetValue("calo_met_phi", _calo_met_phi);

  //"Cleaned" leptons
  //Four vectors only
  std::vector <double> cleanedElPt;
  std::vector <double> cleanedElEta;
  std::vector <double> cleanedElPhi;
  std::vector <double> cleanedElEnergy;

  std::vector <double> cleanedMuPt;
  std::vector <double> cleanedMuEta;
  std::vector <double> cleanedMuPhi;
  std::vector <double> cleanedMuEnergy;


  for (std::vector<int>::const_iterator iel = vSelElectrons.begin();
       iel != vSelElectrons.end(); ++iel){   

    Electron * _el = (Electron *) b_Electron->At(*iel);

    cleanedElPt     . push_back(_el->PT);
    cleanedElEta    . push_back(_el->Eta);
    cleanedElPhi    . push_back(_el->Phi);
    cleanedElEnergy . push_back(_el->P4().Energy());        
  }

  for (std::vector<int>::const_iterator imu = vSelMuons.begin();
       imu != vSelMuons.end(); ++imu){   
    
    Muon * _mu = (Muon *) b_Muon->At(*imu);

    cleanedMuPt     . push_back(_mu->PT);
    cleanedMuEta    . push_back(_mu->Eta);
    cleanedMuPhi    . push_back(_mu->Phi);
    cleanedMuEnergy . push_back(_mu->P4().Energy());        
  }

  SetValue("cleanedElPt"     , cleanedElPt);
  SetValue("cleanedElEta"    , cleanedElEta);
  SetValue("cleanedElPhi"    , cleanedElPhi);
  SetValue("cleanedElEnergy" , cleanedElEnergy);

  SetValue("cleanedMuPt"     , cleanedMuPt);
  SetValue("cleanedMuEta"    , cleanedMuEta);
  SetValue("cleanedMuPhi"    , cleanedMuPhi);
  SetValue("cleanedMuEnergy" , cleanedMuEnergy);


  //
  //_____ Gen Info ______________________________
  //

  //Four vector
  std::vector <double> genPt;
  std::vector <double> genEta;
  std::vector <double> genPhi;
  std::vector <double> genEnergy;

  //Identity
  std::vector <int> genID;
  std::vector <int> genIndex;
  std::vector <int> genStatus;
  std::vector <int> genMotherID;
  std::vector <int> genMotherIndex;

  double higgsWeight = 1.0;
/*
  double eventWeight = 0.0;
  LHEFEvent* event=(LHEFEvent*)b_Event->At(0);
  eventWeight=event->Weight; 
  cout << "eventWeight = " << eventWeight << endl;
*/

  /*

  if (isMc){
    
  //From Mike Luk: April 10th 2013 
      double higgsBBSf     = 1.47415;
      double higgsTauTauSf = 1.32511;
      double higgsMuMuSf   = 1.30178;
      double higgsCCSf     = 1.35842;

      double higgsGGSf         = 0.25024;
      double higgsGammaGammaSf = 5.55457;
      double higgsZGammaSf     = 1.61765;
      double higgsWWSf = 1.22012;
      double higgsZZSf = 1.38307;


    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByLabel(genParticles_it, genParticles);
    
    for(size_t i = 0; i < genParticles->size(); i++){
      const reco::GenParticle & p = (*genParticles).at(i);

      int id = p.pdgId();
      // find higgs (25)                                                                                                            
      if(abs(id) == 25 && p.status() == 3){
        int absDauIds = 0;
        size_t nDaughters = p.numberOfDaughters();
        // get all daughters                                                                                                        
        for(size_t j = 0; j < nDaughters; ++ j) {
          int dauId = (p.daughter(j))->pdgId();
          absDauIds += abs(dauId);
        }// daughters                                                                                                               

        // for each higgs, find the decay products and weight accordingly                                                           
        if(absDauIds==10) higgsWeight *= higgsBBSf;
        if(absDauIds==30) higgsWeight *= higgsTauTauSf;
        if(absDauIds==26) higgsWeight *= higgsMuMuSf;
        if(absDauIds==8)  higgsWeight *= higgsCCSf;
        if(absDauIds==42) higgsWeight *= higgsGGSf;
        if(absDauIds==44) higgsWeight *= higgsGammaGammaSf;
        if(absDauIds==45) higgsWeight *= higgsZGammaSf;
        if(absDauIds==48) higgsWeight *= higgsWWSf;
        if(absDauIds==46) higgsWeight *= higgsZZSf;
      } // if higgs                                                                                                                 

      //Find status 3 particles
      if (p.status() == 3){
	reco::Candidate* mother = (reco::Candidate*) p.mother();
	if (not mother)            continue;
	
	bool bKeep = false;
	for (unsigned int uk = 0; uk < keepMomPDGID.size(); uk++){
	  if (abs(mother->pdgId()) == (int) keepMomPDGID.at(uk)){
	    bKeep = true;
	    break;
	  }
	}
	
	if (not bKeep){
	  for (unsigned int uk = 0; uk < keepPDGID.size(); uk++){
	    if (abs(p.pdgId()) == (int) keepPDGID.at(uk)){
	      bKeep = true;
	      break;
	    }
	  }
	}
	
	if (not bKeep) continue;
	
	//Find index of mother
	int mInd = 0;
	for(size_t j = 0; j < genParticles->size(); j++){
	  const reco::GenParticle & q = (*genParticles).at(j);
	  if (q.status() != 3) continue;
	  if (mother->pdgId() == q.pdgId() and fabs(mother->eta() - q.eta()) < 0.01 and fabs(mother->pt() - q.pt()) < 0.01){
	    mInd = (int) j;
	    break;
	  }
	}
	
	//Four vector
	genPt     . push_back(p.pt());
	genEta    . push_back(p.eta());
	genPhi    . push_back(p.phi());
	genEnergy . push_back(p.energy());
	
	//Identity
	genID            . push_back(p.pdgId());
	genIndex         . push_back((int) i);
	genStatus        . push_back(p.status());
	genMotherID      . push_back(mother->pdgId());
	genMotherIndex   . push_back(mInd);
      }
    }//End loop over gen particles
  }  //End MC-only if

  */

  //Four vector
  SetValue("genPt"     , genPt);
  SetValue("genEta"    , genEta);
  SetValue("genPhi"    , genPhi);
  SetValue("genEnergy" , genEnergy);

  //Identity
  SetValue("genID"            , genID);
  SetValue("genIndex"         , genIndex);
  SetValue("genStatus"        , genStatus);
  SetValue("genMotherID"      , genMotherID);
  SetValue("genMotherIndex"   , genMotherIndex);
  SetValue("higgsWeight",higgsWeight);


  return 0;
}

int DileptonCalc::findMatch(const reco::GenParticleCollection & genParticles,
			    int idToMatch, double eta, double phi){

  float dRtmp = 1000;
  float closestDR = 10000.;
  int closestGenPart = -1;

  for(size_t j = 0; j < genParticles.size(); ++ j) {
    const reco::GenParticle & p = (genParticles).at(j);
    dRtmp = mdeltaR( eta, phi, p.eta(), p.phi());
    if ( dRtmp < closestDR && abs(p.pdgId()) == idToMatch){// && dRtmp < 0.3) {
      closestDR = dRtmp;
      closestGenPart = j;
    }//end of requirement for matching
  }//end of gen particle loop 
  return closestGenPart;
}


double DileptonCalc::mdeltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}

void DileptonCalc::fillMotherInfo(TClonesArray * genparticles, Int_t mInd, int i,
				  vector <int> & momid, vector <int> & momstatus,
				  vector<double> & mompt, vector<double> & mometa,
				  vector<double> & momphi, vector<double> & momenergy)
{

  GenParticle * mother = (GenParticle *)genparticles->At(mInd);

  if(mother) {
    momid.push_back(mother->PID);
    momstatus.push_back(mother->Status);
    mompt.push_back(mother->PT);
    mometa.push_back(mother->Eta);
    momphi.push_back(mother->Phi);
    momenergy.push_back(mother->E);
    if(i<10)fillMotherInfo(genparticles, mother->M1, i+1, momid, momstatus, mompt, mometa, momphi, momenergy);
  }


} 
