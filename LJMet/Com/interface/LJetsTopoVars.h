/* -*- mode: c++ -*- */

#ifndef LJETSTOPOVARS
#define LJETSTOPOVARS

#include <string>
#include <vector>
#include "LJMet/Com/interface/TMBLorentzVector.h"
#include "TVectorD.h"
#include "TLorentzVector.h"


class LJetsTopoVars
{
public:
  LJetsTopoVars()
    :_ht(          std::vector<double>(22, 0.) ),
     _htOK(        false),
     _evtTopo(     std::vector<double>( 3, 0.) ),
     _evtTopoOK(   false),
     _kt(          std::vector<double>( 3, 0.) ),
     _ktOK(        false),
     _mt(          std::vector<double>( 2, 0.) ),
     _mtOK(        false){};

  LJetsTopoVars(std::vector<TLorentzVector> & jets,
		TLorentzVector & lepton,
		TLorentzVector & met,
		bool isMuon,
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

    setEvent(jets, lepton, met, isMuon);

  };

  virtual ~LJetsTopoVars(){ };
	
  // Initiate LJetsTopoVars using one lepton, one MET and
  // 4 leading jets momenta. Note that if fewer than 4 jets are supplied,
  // some variables are not well-defined. Every effort is made to process
  // such situations correctly. Still, the user should use caution.
  int setEvent(std::vector<TLorentzVector> & jets,
	       TLorentzVector & lepton,
	       TLorentzVector & met,
	       bool isMuon);

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

private:
  std::vector<TMBLorentzVector> m_jets;
  TMBLorentzVector m_met;
  TMBLorentzVector m_lepton;

  int nJets;

  TVectorD eigenval;

  bool m_isMuon;
  double WMassPdg;

  TMBLorentzVector _neutrino;
		
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

  //ClassDef(LJetsTopoVars,1) // L+jets topological and kinematic variables
};

#endif //LJETSTOPOVARS
