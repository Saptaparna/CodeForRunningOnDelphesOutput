#ifndef LJMet_Com_interface_LjmetEventSelector_h
#define LJMet_Com_interface_LjmetEventSelector_h

/*
   Interface class for FWLite PAT analyzer-selectors
   Specific selectors must implement the () operator

   Author: Gena Kukartsev, 2010,2012
*/





#include <cmath>      //necessary for absolute function fabs()
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
//#include "LJMet/Com/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
#include "LJMet/Com/interface/ElectronSelector.h"
#include "LJMet/Com/interface/MetSelectionFunctor.h"
#include "LJMet/Com/interface/MuonSelectionFunctor.h"
#include "LJMet/Com/interface/PVSelector.h"
#include "LJMet/Com/interface/PvObjectSelector.h"

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


// Root includes
#include "TROOT.h"

using namespace std;


class LjmetEventSelector : public EventSelector {
 public:
  

  enum Selection_t { Top, TopBTag, EWK_Zmumu, DiLepton, TprimeInc };
  enum Version_t { Off, test, SelV1, SelV2, SelV3, Tight, Loose, DiLepton_MuMu, DiLepton_ElEl, DiLepton_MuEl, Wprime_Mu, Wprime_Mu_forPrompt, Wprime_El, Tprime_Inc};
  


  LjmetEventSelector( std::map<std::string, edm::ParameterSet const > params){
  }




  LjmetEventSelector( edm::ParameterSet const & muonSelParams,
		      edm::ParameterSet const & looseMuonSelParams,
		      edm::ParameterSet const & electronSelParams,
		      edm::ParameterSet const & looseElectronSelParams,
		      edm::ParameterSet const & jetSelParams,
		      edm::ParameterSet const & metSelParams,
		      edm::ParameterSet const & pvSelParams,
		      edm::ParameterSet const & pvObjSelParams,
		      edm::ParameterSet const & params ) :

    muonSel_     (new MuonSelectionFunctor (muonSelParams)),
    looseMuonSel_(new MuonSelectionFunctor (looseMuonSelParams)),
    electronSel_ (new ElectronSelector     (electronSelParams)),
    looseElectronSel_ (new ElectronSelector     (looseElectronSelParams)),
    jetSel_      (new JetIDSelectionFunctor(jetSelParams)),
    metSel_      (new MetSelectionFunctor  (metSelParams)),
    pvSel_       (new PVSelector           (pvSelParams)),
    pvObjSel_    (new PvObjectSelector     (pvObjSelParams)),
    pvSrc_       (params.getParameter<edm::InputTag>("pvSrc")),
    muonSrc_     (params.getParameter<edm::InputTag>("muonSrc")),
    electronSrc_ (params.getParameter<edm::InputTag>("electronSrc")),
    jetSrc_      (params.getParameter<edm::InputTag>("jetSrc")),
    jetSrcPF_    (params.getParameter<edm::InputTag>("jetSrcPF")),
      //jetSrcJPT_   (params.getParameter<edm::InputTag>("jetSrcJPT")),
    metSrc_      (params.getParameter<edm::InputTag>("metSrc")),
    metTCSrc_    (params.getParameter<edm::InputTag>("metTCSrc")),
    metPFSrc_    (params.getParameter<edm::InputTag>("metPFSrc")),
    metPfTypeISrc_ (params.getParameter<edm::InputTag>("metPfTypeISrc")),
    metMuonSrc_  (params.getParameter<edm::InputTag>("metMuonSrc")),
    metTypeIISrc_(params.getParameter<edm::InputTag>("metTypeIISrc")),
    rhoCorrectionSrc_  (params.getParameter<edm::InputTag>("rhoCorrection")),
    triggerSrc_  (params.getParameter<edm::InputTag>("triggerSrc")),
    BeamspotSrc_ (params.getParameter<edm::InputTag>("BeamspotSrc")) {
    
    std::string selectionStr = params.getParameter<std::string>("selection");
    std::string versionStr = params.getParameter<std::string>("version");
    
  }
    
    
  // legacy constructor for backwards compatibility,
  // only PV event selector is used here   
  LjmetEventSelector( edm::ParameterSet const & muonSelParams,
		      edm::ParameterSet const & looseMuonSelParams,
		      edm::ParameterSet const & electronSelParams,
		      edm::ParameterSet const & looseElectronSelParams,
		      edm::ParameterSet const & jetSelParams,
		      edm::ParameterSet const & metSelParams,
		      edm::ParameterSet const & pvSelParams,
		      edm::ParameterSet const & params ) :

    muonSel_     (new MuonSelectionFunctor (muonSelParams)),
    looseMuonSel_(new MuonSelectionFunctor (looseMuonSelParams)),
    electronSel_ (new ElectronSelector     (electronSelParams)),
    looseElectronSel_ (new ElectronSelector     (looseElectronSelParams)),
    jetSel_      (new JetIDSelectionFunctor(jetSelParams)),
    metSel_      (new MetSelectionFunctor  (metSelParams)),
    pvSel_       (new PVSelector           (pvSelParams)),
    muonSrc_     (params.getParameter<edm::InputTag>("muonSrc")),
    electronSrc_ (params.getParameter<edm::InputTag>("electronSrc")),
    jetSrc_      (params.getParameter<edm::InputTag>("jetSrc")),
    jetSrcPF_    (params.getParameter<edm::InputTag>("jetSrcPF")),
      //jetSrcJPT_   (params.getParameter<edm::InputTag>("jetSrcJPT")),
    metSrc_      (params.getParameter<edm::InputTag>("metSrc")),
    metTCSrc_    (params.getParameter<edm::InputTag>("metTCSrc")),
    metPFSrc_    (params.getParameter<edm::InputTag>("metPFSrc")),
    metPfTypeISrc_ (params.getParameter<edm::InputTag>("metPfTypeISrc")),
    metMuonSrc_  (params.getParameter<edm::InputTag>("metMuonSrc")),
    metTypeIISrc_(params.getParameter<edm::InputTag>("metTypeIISrc")),
    rhoCorrectionSrc_  (params.getParameter<edm::InputTag>("rhoCorrection")),
    triggerSrc_  (params.getParameter<edm::InputTag>("triggerSrc")),
    BeamspotSrc_ (params.getParameter<edm::InputTag>("BeamspotSrc")) {
    
    std::string selectionStr = params.getParameter<std::string>("selection");
    std::string versionStr = params.getParameter<std::string>("version");
    
  }
    
    
  virtual ~LjmetEventSelector();
  
  
  
  virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret) = 0;



  bool areMatched ( const reco::Candidate & c1,
		    const reco::Candidate & c2,
		    double DR,
		    double DPtRel ) {
    unsigned int nPass=0;
    if (  (deltaR(c1, c2) < DR)   && (fabs(c2.pt() - c1.pt())/ c2.pt()<DPtRel)) {
      nPass++;
    }
    return (nPass>0);
  }



  boost::shared_ptr<MetSelectionFunctor>   const & metSel()       const { return metSel_;}
  boost::shared_ptr<JetIDSelectionFunctor> const & jetSel()       const { return jetSel_;}
  boost::shared_ptr<MuonSelectionFunctor>  const & muonSel()      const { return muonSel_;}
  boost::shared_ptr<MuonSelectionFunctor>  const & looseMuonSel() const { return looseMuonSel_;}
  boost::shared_ptr<ElectronSelector>      const & electronSel()  const { return electronSel_;}
  boost::shared_ptr<ElectronSelector>      const & looseElectronSel()  const { return looseElectronSel_;}
  boost::shared_ptr<PVSelector>            const & pvSel()        const { return pvSel_;}
  boost::shared_ptr<PvObjectSelector>      const & pvObjSel()     const { return pvObjSel_;}
  
  pat::MET      const & met()      const { return *met_; }
  pat::MET      const & tcMet()    const { return *tcMet_; }
  pat::MET      const & pfMet()    const { return *pfMet_; }
  reco::MET     const & pfTypeIMet() const { return *pfTypeIMet_; }
  reco::CaloMET const & caloMet()  const { return *caloMet_; }
  reco::CaloMET const & caloMet2() const { return *caloMet2_; }

  vector<pat::Jet>            const & allJets()  const { return *h_jets_; }
  vector<edm::Ptr<reco::Vertex> >  const   selected_PVs() const { return good_pvs_; }
  vector<edm::Ptr<pat::Jet> >  const   goodJets() const { return good_jets_; }
  vector<edm::Ptr<pat::Jet> >  const   goodJetsPF() const { return good_jetsPF_; }
  vector<edm::Ptr<pat::Jet> >  const   allJetsPF() const { return all_jetsPF_; }
  vector<edm::Ptr<pat::Jet> >  const   goodJetsJPT() const { return good_jetsJPT_; }
  vector<edm::Ptr<pat::Muon> > const  goodMuons() const { return good_muons_; }
  vector<edm::Ptr<pat::Electron> > const  goodElectrons() const { return good_electrons_; }
  pat::Jet                    const & jet0()     const { return *jet0_; }
  pat::Jet                    const & jet1()     const { return *jet1_; }
  pat::Jet                    const & jet2()     const { return *jet2_; }
  pat::Jet                    const & jet3()     const { return *jet3_; }

  vector<pat::Muon> const & allMuons() const { return *h_muons_; }
  pat::Muon         const & muon0()    const { return *muon0_; }
  pat::Muon         const & muon1()    const { return *muon1_; }
  double                    dimuonMass()     { return _dimuon_mass; }

  vector<pat::Electron> const & allElectrons () const { return *h_electrons_; }
  pat::Electron         const & electron0()     const { return *electron0_; }
  pat::Electron         const & electron1()    const { return *electron1_; }


protected:

  boost::shared_ptr<MuonSelectionFunctor>  muonSel_;
  boost::shared_ptr<MuonSelectionFunctor>  looseMuonSel_;
  boost::shared_ptr<ElectronSelector>      electronSel_;
  boost::shared_ptr<ElectronSelector>      looseElectronSel_;
  boost::shared_ptr<JetIDSelectionFunctor> jetSel_;
  boost::shared_ptr<MetSelectionFunctor>   metSel_;
  boost::shared_ptr<PVSelector>            pvSel_;
  boost::shared_ptr<PvObjectSelector>      pvObjSel_;

  edm::InputTag pvSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag jetSrcPF_;
  //edm::InputTag jetSrcJPT_;
  edm::InputTag metSrc_;
  edm::InputTag metTCSrc_;
  edm::InputTag metPFSrc_;
  edm::InputTag metPfTypeISrc_;
  edm::InputTag metMuonSrc_;
  edm::InputTag metTypeIISrc_;
  edm::InputTag rhoCorrectionSrc_;
  edm::InputTag triggerSrc_;
  edm::InputTag BeamspotSrc_;


  edm::Handle<std::vector<reco::Vertex> >     h_primVtx;
  edm::Handle<vector<pat::MET> >              h_met_;
  edm::Handle<vector<pat::MET> >              h_tcMet_;
  edm::Handle<vector<pat::MET> >              h_pfMet_;
  edm::Handle<vector<reco::MET> >             h_pfTypeIMet_;
  edm::Handle<vector<reco::CaloMET> >         h_caloMet_;
  edm::Handle<vector<reco::CaloMET> >         h_caloMet2_;
  edm::Handle<vector<pat::Jet> >              h_jets_;
  edm::Handle<vector<pat::Jet> >              h_jetsPF_;
  edm::Handle<vector<pat::Jet> >              h_jetsJPT_;
  edm::Handle<vector<pat::Muon> >             h_muons_;
  edm::Handle<vector<pat::Electron> >         h_electrons_;
  edm::Handle<double>                         h_rho;
  edm::Handle<std::vector<pat::TriggerPath> > hltresults_;
  edm::Handle<edm::TriggerResults>            hlt_;
  edm::Handle<pat::TriggerEvent>              triggerEvent;
  edm::Handle<pat::TriggerEvent >             h_triggerEvent_;
  edm::Handle<edm::TriggerResults >           mhEdmTriggerResults;
  edm::Handle<reco::BeamSpot >                Beamspot_;


  edm::Ptr<pat::MET>      met_;
  edm::Ptr<pat::MET>      tcMet_;
  edm::Ptr<pat::MET>      pfMet_;
  edm::Ptr<reco::MET>     pfTypeIMet_;
  edm::Ptr<reco::CaloMET> caloMet_;
  edm::Ptr<reco::CaloMET> caloMet2_;
  edm::Ptr<pat::Jet>      jet0_;
  edm::Ptr<pat::Jet>      jet1_;
  edm::Ptr<pat::Jet>      jet2_;  
  edm::Ptr<pat::Jet>      jet3_;
  edm::Ptr<pat::Muon>     muon0_;
  edm::Ptr<pat::Muon>     muon1_;
  edm::Ptr<pat::Electron> electron0_;
  edm::Ptr<pat::Electron> electron1_;

  double                           _dimuon_mass;

  std::vector<edm::Ptr<reco::Vertex> > good_pvs_;
  std::map<std::string,bool>       _hlt;
  std::vector<edm::Ptr<pat::Jet> > good_jets_;
  std::vector<edm::Ptr<pat::Jet> > good_jetsPF_;
  std::vector<edm::Ptr<pat::Jet> > all_jetsPF_;
  std::vector<edm::Ptr<pat::Jet> > good_jetsJPT_;
  std::vector<edm::Ptr<pat::Muon> > good_muons_;
  std::vector<edm::Ptr<pat::Muon> > loose_muons_;
  std::vector<edm::Ptr<pat::Electron> > good_electrons_;
  std::vector<edm::Ptr<pat::Electron> > loose_electrons_;

  bool muonMatchedToHlt_firstpass;

  Selection_t selection;
  Version_t version;

};

LjmetEventSelector::~LjmetEventSelector(){}

#endif
