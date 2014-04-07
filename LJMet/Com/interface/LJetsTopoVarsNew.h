/* -*- mode: c++ -*- */

#ifndef LJETSTOPOVARSNEW
#define LJETSTOPOVARSNEW

#include <string>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "LJMet/Com/interface/TMBLorentzVector.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "LJMet/Com/interface/METzCalculator.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

class LJetsTopoVarsNew
{
public:
  LJetsTopoVarsNew()
    :_ht(          std::vector<double>(22, 0.) ),
     _htOK(        false),
     _evtTopo(     std::vector<double>( 3, 0.) ),
     _evtTopoOK(   false),
     _kt(          std::vector<double>( 3, 0.) ),
     _ktOK(        false),
     _mt(          std::vector<double>( 2, 0.) ),
     _mtOK(        false){};

    //LJetsTopoVarsNew(std::vector<TLorentzVector> & jets,
    LJetsTopoVarsNew(const std::vector<std::pair<TLorentzVector,bool> > jets, 
                     TLorentzVector & lepton,
                     TLorentzVector & met,
                     bool isMuon,
                     bool bestTop,
                     double wmass = 80.398)
        :m_isMuon(isMuon),
         WMassPdg(wmass),
         _ht(          std::vector<double>(22, 0.) ),
         _htOK(        false),
         _evtTopo(     std::vector<double>( 3, 0.) ),
         _evtTopoOK(   false),
         _kt(          std::vector<double>( 3, 0.) ),
         _ktOK(        false),
         _mt(          std::vector<double>( 2, 0.) ),
         _mtOK(        false){
 
       setEvent(jets, lepton, met, isMuon, bestTop);

  };

  virtual ~LJetsTopoVarsNew(){ }; 
	
  // Initiate LJetsTopoVarsNew using one lepton, one MET and
  // 4 leading jets momenta. Note that if fewer than 4 jets are supplied,
  // some variables are not well-defined. Every effort is made to process
  // such situations correctly. Still, the user should use caution.
  //int setEvent(std::vector<TLorentzVector> & jets,
    int setEvent(const std::vector<std::pair<TLorentzVector,bool> > jets,
               TLorentzVector & lepton,
               TLorentzVector & met,
               bool isMuon, 
               bool bestTop);

  int setEventMetFixed(TLorentzVector&,TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&,double min_dr_jet_lepton=-0.01);
    

  double aplanarity() const;
  double centrality() const;
  double sphericity() const;
  double ht() const;
  double htpluslepton() const;
  double methtpluslepton() const;
  double h() const;
  double ktMinPrime() const;
  double dphiLepMet() const;
  double minDijetMass() const;
  double maxJetEta() const;
  double Et3() const;
  double minDijetDeltaR() const;
  double LeptonJet_DeltaR();           //btw lepton and leading 2 jet, minimum
  double Jet1Jet2_DeltaR();     
  double Jet1Jet2_DeltaPhi();   
  double Jet1Jet2_M();     
  double Jet1Jet2_Pt();   
  double Jet1Jet2W_M();     
  double Jet1Jet2W_Pt();   
		
  double Hz();        //scalar sum of longitudinal energies of first four jets, the muon, and the neutrino
  double HT2();       //scalar sum of transverse energies of the 2nd, 3rd, and 4th jets
  double HT2prime();  //HT2/Hz
  double W_MT();      //W transverse mass
  double W_M();       //W mass
  double W_Pt();       //W 
  double DphiJMET();  //Delta phi btw leading jet and MET

  //muon only (throws exception if not)
  double Muon_DeltaR();           //btw muon and jet, minimum
  //double Muon_etHaloScaled();     //TMBMuon::etHalo()/pT
  //double Muon_etTrkConeScaled();  //TMBMuon::etTrkCone()/pT

  //electron only (throws exception if not)
  //double Electron_iso();    //TMBEMCluster::iso()
  //double Electron_lhood();  //TMBEMCluster::Lhood8(); is this right?
		
  // Ht
  double getHt()        {if(!_htOK) calcHt(); return _ht[ 0];} 
  double getHtp()       {if(!_htOK) calcHt(); return _ht[ 1];} 
  double getHtpp()      {if(!_htOK) calcHt(); return _ht[ 2];}
  double getHt2()       {if(!_htOK) calcHt(); return _ht[ 3];}
  double getHt2p()      {if(!_htOK) calcHt(); return _ht[ 4];}
  double getHt2pp()     {if(!_htOK) calcHt(); return _ht[ 5];} 
  double getHt3()       {if(!_htOK) calcHt(); return _ht[ 6];}
  double getHt3p()      {if(!_htOK) calcHt(); return _ht[ 7];} 
  double getHt3pp()     {if(!_htOK) calcHt(); return _ht[ 8];} 
  double getCen()       {if(!_htOK) calcHt(); return _ht[ 9];} 
  double getNJW()       {if(!_htOK) calcHt(); return _ht[10];}
  double getJetEtaMax() {if(!_htOK) calcHt(); return _ht[11];}
  double getMdijetMin() {if(!_htOK) calcHt(); return _ht[12];}
  double getMtjets()    {if(!_htOK) calcHt(); return _ht[13];}
  double getSqrtsT()    {if(!_htOK) calcHt(); return _ht[14];}
  double getMtAurelio() {if(!_htOK) calcHt(); return _ht[15];}
  double getPzOverHT()  {if(!_htOK) calcHt(); return _ht[16];}
  double getMevent()    {if(!_htOK) calcHt(); return _ht[17];}
  double getM123inv()   {if(!_htOK) calcHt(); return _ht[18];}
  double getEta2Sum()   {if(!_htOK) calcHt(); return _ht[19];}
  double getMwRec()     {if(!_htOK) calcHt(); return _ht[20];}
  double getH()         {if(!_htOK) calcHt(); return _ht[21];}

  // event topo
  double getSph()   {if(!_evtTopoOK) calcEvtTopo(); return _evtTopo[0];}
  double getApl()   {if(!_evtTopoOK) calcEvtTopo(); return _evtTopo[1];}
  double getAplMu() {if(!_evtTopoOK) calcEvtTopo(); return _evtTopo[2];}

  // Kt
  double getKtminp()        {if(!_ktOK) calcKt(); return _kt[0];}
  double getKtminpReduced() {if(!_ktOK) calcKt(); return _kt[1];}
  double getDrMinJetJet()   {if(!_ktOK) calcKt(); return _kt[2];}

  // mT
  double getDphiMuMet() {if(!_mtOK) calcMt(); return _mt[0];}
  double getMt()        {if(!_mtOK) calcMt(); return _mt[1];}

  // Momentum tensor eigenvalues
  TVectorD getEigen() {if(!_evtTopoOK) calcEvtTopo(); return eigenval;}

  // Neutrino
  TMBLorentzVector GetNeutrino() {return _neutrino;}

  // From W' analysis
  double dphiLepJ1() ;
  double dphiLepJ2() ;
  double dphiLepJ3() ;
  double dphiLepJ4() ;
  double dphiLepLeadBTagJet() ;
  double dphiLepSecLeadBTagJet() ;
  double dphiLepLightJet() ;
  double LeadBTagJet_DiscVal();
  double SecLeadBTagJet_DiscVal();
  double Lead12BTagJet_DiscVal();
  double LightJet_DiscVal();
  double BestJetJet2W_M(); //new    
  double Jet1TagJet2TagW_M(); //new    
  double Cos_LightjetJetLepton_BestTop();
  double Cos_LightjetJetLepton_BTagTop();
  double Cos_BestJetLepton_BestTop();		
   
  double SecBestTop();//new
  double SecBestBTagTop();//new
  double BestTop();
  double BestTop_Pt();
  double BestJet_Pt();
  double BestJet_Eta();
  double BestJet_Phi();
  double BTagTopMass();
  double BTagTop_Pt();
  double SecBTagTopMass();
  double SecBTagTop_Pt();  
  double HT_AllJets_MinusBestJet();
  double H_AllJets_MinusBestJet();
  double J1_NotBestJet_Eta();
  double J2_NotBestJet_Eta();
  double AllJets_MinusBestJet_Pt();
  double AllJets_M();
  double AllJetsW_M();//sqrt_shat
  double J1_NotBestJet_Pt();
  double J2_NotBestJet_Pt();
  double J1_NotBestJet_Phi();

  int LeptonMETxy();

  void SetBestTop_JetIndex(unsigned int index) {
    _BestTop_JetIndex = index;
  }
  unsigned int  GetBestTop_JetIndex() {
    return _BestTop_JetIndex;
  }
  //
  void SetBestTop(TMBLorentzVector BestTop) {
    _BestTop = BestTop;
  }
  TMBLorentzVector  GetBestTop() {
    return _BestTop;
  }
	//
  void SetBTagTop(TMBLorentzVector BTagTop) {
    _TopLeadingBTaggedJet = BTagTop;
  }
  TMBLorentzVector  GetBTagTop() {
    return _TopLeadingBTaggedJet;
  }
  void SetSecBTagTop(TMBLorentzVector SecBTagTop) {
    _TopSecLeadingBTaggedJet = SecBTagTop;
  }
  TMBLorentzVector  GetSecBTagTop() {
    return _TopSecLeadingBTaggedJet;
  }
  //
  void SetGoodJetsMinusBestJet(std::vector<TMBLorentzVector> GoodJetsMinusBestJet) {
    _GoodJetsMinusBestJet = GoodJetsMinusBestJet;
  }
  std::vector<TMBLorentzVector> GetGoodJetsMinusBestJet() {
    return _GoodJetsMinusBestJet;
  }
  
  void SetLeptonMETxy(std::vector<TMBLorentzVector> LeptonMETxy) {
    _LeptonMETxy = LeptonMETxy;
  }
  
  std::vector<TMBLorentzVector> GetLeptonMETxy() {
    return _LeptonMETxy;
  }

  //
  void SetGoodJetsMinusLeadingBTaggedJet(std::vector<TMBLorentzVector> GoodJetsMinusLeadingBTaggedJet) {
    _GoodJetsMinusLeadingBTaggedJet = GoodJetsMinusLeadingBTaggedJet;
  }
  std::vector<TMBLorentzVector> GetGoodJetsMinusLeadingBTaggedJet() {
    return _GoodJetsMinusLeadingBTaggedJet;
  }


private:
  std::vector<TMBLorentzVector> m_jets;
  TMBLorentzVector m_met;
  TMBLorentzVector m_lepton;
  TMBLorentzVector _neutrino;

  int nJets;

  unsigned int _BestTop_JetIndex; // index of the jet that gives best top mass
  std::vector<TMBLorentzVector> _GoodJetsMinusBestJet;
  std::vector<TMBLorentzVector> _LeptonMETxy;
  TMBLorentzVector _BestTop;
  TMBLorentzVector _TopLeadingBTaggedJet;
  TMBLorentzVector _TopSecLeadingBTaggedJet;
  std::vector<TMBLorentzVector> _GoodJetsMinusLeadingBTaggedJet;
  int number_of_jets ;
  int number_of_tagged_jets ;
  int number_of_untagged_jets ;
	
  unsigned int tagged_jet_highpt_index;
  unsigned int second_tagged_jet_highpt_index;
  unsigned int untagged_jet_highpt_index;
  unsigned int second_untagged_jet_index;

  double m_discCutValue;
  std::vector<double> m_v_discCutValue;
  double tagged_jet_highpt_DiscVal;
  double second_tagged_jet_highpt_DiscVal;
  double untagged_jet_highpt_DiscVal;
  double second_untagged_jet_DiscVal;
   	
  TVectorD eigenval;

  bool m_isMuon;
  double WMassPdg;

		
  void calcHt();
  std::vector<double> _ht;
  bool _htOK;

  void calcEvtTopo();
  std::vector<double> _evtTopo;
  bool _evtTopoOK;

  void calcKt();
  std::vector<double> _kt;
  bool _ktOK;

  void calcMt();
  std::vector<double> _mt;
  bool _mtOK;

  METzCalculator fzCalculator;


  //ClassDef(LJetsTopoVarsNew,1) // L+jets topological and kinematic variables
};

#endif //LJETSTOPOVARSNEW
