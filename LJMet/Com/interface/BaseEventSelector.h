#ifndef LJMet_Com_interface_BaseEventSelector_h
#define LJMet_Com_interface_BaseEventSelector_h

/*
   Interface class for FWLite PAT analyzer-selectors
   Specific selectors must implement the () operator

   Author: Gena Kukartsev, 2010,2012
*/

#include <cmath>
#include <iostream>

#include "TROOT.h"
#include "TVector3.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
#include "LJMet/Com/interface/BTagSFUtil.h"
#include "LJMet/Com/interface/BtagHardcodedConditions.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "FWCore/Framework/interface/Event.h"

#include "LJMet/DelphesData/interface/ExRootTreeReader.h"
#include "LJMet/DelphesData/interface/DelphesClasses.h"

class BaseEventSelector : public EventSelector {
  //
  // Base class for all event selector plugins
  //


  friend class LjmetFactory;


 public:
  

  BaseEventSelector();
    
    
  virtual ~BaseEventSelector();
  

  virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
  //virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret) = 0;
  virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret){return false;} // dummy
  virtual bool operator()( ExRootTreeReader * reader, Int_t entry, pat::strbitset & ret) = 0;
  virtual void EndJob();


  //virtual void AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec );
  virtual void AnalyzeEvent( ExRootTreeReader * reader, Int_t entry, LjmetEventContent & ec );


  std::string  GetName();
  double       GetPerp(TVector3 & v1, TVector3 & v2);
  bool         AreMatched ( const reco::Candidate & c1,
			    const reco::Candidate & c2,
			    double DR,
			    double DPtRel );
  
  
  std::vector<int>            const & GetAllJets()           const;
  std::vector<int>            const & GetSelectedJets()      const;
  std::vector<int>            const & GetSelectedCAJets()      const;
  std::vector<int>            const & GetLooseJets()         const;
  std::vector<int>            const & GetSelectedBtagJets()  const;
  std::vector<std::pair<TLorentzVector,bool>> const & GetCorrJetsWithBTags()  const;
  std::vector<edm::Ptr<pat::Muon> >           const & GetAllMuons()          const;
  std::vector<int>           const & GetSelectedMuons()     const;
  std::vector<edm::Ptr<pat::Muon> >           const & GetLooseMuons()        const;
  std::vector<edm::Ptr<pat::Electron> >       const & GetAllElectrons()      const;
  std::vector<int>       const & GetSelectedElectrons() const;
  std::vector<edm::Ptr<pat::Electron> >       const & GetLooseElectrons()    const;
  MissingET *                 const GetMet()               const;
  edm::Ptr<reco::PFMET>                       const & GetType1CorrMet()      const;
  TLorentzVector                              const & GetCorrectedMet()      const;
  std::vector<unsigned int>                   const & GetSelectedTriggers()  const;
  std::vector<edm::Ptr<reco::Vertex> >        const & GetSelectedPVs()  const;
  double const & GetTestValue() const;


  void SetMc(bool isMc);
  bool IsMc();


  // LJMET event content setters
  void Init( void );
  void SetEventContent(LjmetEventContent * pEc);
  void SetHistogram(std::string name, int nbins, double low, double high);
  void SetHistValue(std::string name, double value);
  void SetTestValue(double & test);

  void SetCorrectedMet(TLorentzVector & met);
  void SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector,bool>> & jets);

  bool isJetTagged(int iJet, ExRootTreeReader * reader, Int_t entry, bool applySF = true);
  TLorentzVector correctJet(int iJet, ExRootTreeReader * reader, Int_t entry);
  TLorentzVector correctMet( ExRootTreeReader * reader, Int_t entry);

 protected:

  std::vector<int>      mvAllJets;
  std::vector<int>      mvSelJets;
  std::vector<int>      mvAllCAJets;
  std::vector<int>      mvSelCAJets;
  std::vector<int>      mvLooseJets;
  std::vector<std::pair<TLorentzVector,bool> >      mvCorrJetsWithBTags;
  std::vector<int>      mvSelBtagJets;
  std::vector<int>      mvSelBtagCAJets;
  std::vector<edm::Ptr<pat::Muon> >     mvAllMuons;
  std::vector<int>                  mvSelMuons;
  std::vector<edm::Ptr<pat::Muon> >     mvLooseMuons;
  std::vector<edm::Ptr<pat::Electron> > mvAllElectrons;
  std::vector<int> mvSelElectrons;
  std::vector<edm::Ptr<pat::Electron> > mvLooseElectrons;
  MissingET *           mpMet;
  edm::Ptr<reco::PFMET>                 mpType1CorrMet;
  TLorentzVector                        correctedMET_p4;
  std::vector<unsigned int>             mvSelTriggers;
  std::vector<edm::Ptr<reco::Vertex> >  mvSelPVs;
  double                                mTestValue;

  // containers for config parameter values
  std::map<std::string,bool>           mbPar;
  std::map<std::string,int>            miPar;
  std::map<std::string,double>         mdPar;
  std::map<std::string,std::string>    msPar;
  std::map<std::string, edm::InputTag> mtPar;
  std::map<std::string,std::vector<std::string> > mvsPar;

  std::string mName;
  std::string mLegend;

  bool mbIsMc;
  


 private:

  void init();
  void setName(std::string name);
  void BeginEvent( ExRootTreeReader * reader, Int_t entry, LjmetEventContent & ec );
  void EndEvent( ExRootTreeReader * reader, Int_t entry, LjmetEventContent & ec );

  BTagSFUtil mBtagSfUtil;
  BtagHardcodedConditions mBtagCond;
  double bTagCut;
  JetCorrectionUncertainty *jecUnc;
  FactorizedJetCorrector *JetCorrector;

  LjmetEventContent * mpEc;

  int mNCorrJets;
  int mNBtagSfCorrJets;
};



#endif
