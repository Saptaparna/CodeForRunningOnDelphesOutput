#include "LJMet/Com/interface/LJetsTopoVars.h"
#include "LJMet/Com/interface/TopTopologicalVariables.h"
#include "LJMet/Com/interface/AnglesUtil.h"

#include "TMatrixDSymEigen.h"

#include <iostream>
#include <stdexcept>


using namespace std;


// Initiate LJetsTopoVars using one lepton, one MET and
// 4 leading jets momenta. Note that if fewer than 4 jets are supplied,
// some variables are not well-defined. Every effort is made to process
// such situations correctly. Still, the user should use caution.
//
// Returns number of jets used: 4 or less. If more than 4 jets were
// supplied, only the first 4 are used, the rest are ignored.
int LJetsTopoVars::setEvent(std::vector<TLorentzVector> & jets,
			    TLorentzVector & lepton,
			    TLorentzVector & met,
			    bool isMuon){
  nJets = 0; // will return this as result

  m_jets.clear();
  
  eigenval.ResizeTo(3);
  eigenval.Zero();
  
  m_met = TMBLorentzVector(met);
   
  m_lepton = TMBLorentzVector(lepton);

  // loop over jets
  for (std::vector<TLorentzVector>::const_iterator jet=jets.begin(); (jet!=jets.end()) && (m_jets.size()!=4); jet++){
    TMBLorentzVector _j(*jet);
    m_jets.push_back(_j);
    ++nJets;
  }

  double nu_px = m_met.Px();
  double nu_py = m_met.Py();

  //set all OK flags to FALSE;
  _htOK = false;
  _evtTopoOK = false;
  _ktOK = false;
  _mtOK = false;

  //
  // calculate neutrino lorentz vector (from Tobi's TopSvtAnalysis)
  //
  double nu_pz = 0.;
  double nu_e  = sqrt(pow(nu_px,2)+pow(nu_py,2));
  
  double Mw    = WMassPdg;
  double l_px  = m_lepton.Px();
  double l_py  = m_lepton.Py();
  double l_pz  = m_lepton.Pz();
  double l_pt  = m_lepton.Pt();
  double l_e   = m_lepton.E();
  double Mt    = sqrt(pow(l_pt+nu_e ,2)-
		      pow(l_px+nu_px,2)-
		      pow(l_py+nu_py,2));
  
  double A;

  // FIXME: do we want this Mt fixed to Mw stuff?
  if (Mt<Mw) A = pow(Mw,2)/2.;
  else       {           // assume Mt=Mw, rescale MET accordingly
    A = pow(Mt,2)/2.;
    double k = nu_e*l_pt - nu_px*l_px - nu_py*l_py;
    k = (k == 0. ? 0.00001 : k);
    double scf = 0.5*pow(Mw,2)/k ;
    nu_px *= scf;
    nu_py *= scf;
    nu_e = sqrt(pow(nu_px,2)+pow(nu_py,2));
  }    
  
  double B = nu_px*l_px + nu_py*l_py;
  double C = TMath::Max(1. + pow(nu_e,2) * (pow(l_pz,2)-pow(l_e,2)) / pow(A+B,2)  , 0.);
  C = sqrt(C);
  double S1= (-(A+B)*l_pz + (A+B)*l_e*C) / (pow(l_pz,2)-pow(l_e,2));
  double S2= (-(A+B)*l_pz - (A+B)*l_e*C) / (pow(l_pz,2)-pow(l_e,2));
  
  // choose solution with smallest |l_pz| a la Run I
  nu_pz = fabs (S1) < fabs (S2) ? S1 : S2 ;
  
  //NOTE: neutrino PX, PY are not necessarily metPX, metPY any more!!!
  _neutrino.SetPxPyPzE(nu_px,nu_py,nu_pz,nu_e);
  
  return nJets;
}


int LJetsTopoVars::setEventMetFixed(TLorentzVector& Jet1, TLorentzVector& Jet2, TLorentzVector& Jet3, TLorentzVector& Jet4, TLorentzVector& NewMet,TLorentzVector& Muon1, double min_dr_jet_lepton)
{
 
  using namespace std;

  cout<<"!!!!!!!!!!!! Im inside TopoVar"<<endl;

  int removed_jets = 0;

  m_jets.clear();
  eigenval.ResizeTo(3);
  eigenval.Zero();

   m_met = TMBLorentzVector(NewMet[0],NewMet[1],NewMet[2],NewMet[3],TMBLorentzVector::kXYZE);
  // cout<<"!!! NewMet eneryg = "<<m_met(3)<<endl;
  m_lepton= TMBLorentzVector(Muon1[0],Muon1[1],Muon1[2],Muon1[3],TMBLorentzVector::kXYZE);
  TMBLorentzVector jets[4] = {Jet1,Jet2,Jet3,Jet4};
 
  // cout << "jets_1_px = "<<jets[0][0];
  // cout << "jets_3_energy = "<<jets[2][3];

   for (int i = 0; i<4; i++){
    //cout << "LJetsTopoVars::setEvent(): jet pt() = " << jet -> pt() << endl;
    TMBLorentzVector _j(jets[i][0],jets[i][1],jets[i][2],jets[i][3],TMBLorentzVector::kXYZE);
    if (m_lepton.DeltaR(_j) > min_dr_jet_lepton){
      m_jets.push_back(_j);
      //jet++;
      //  cout << "!!!!!!! i = "<<i<<endl;
    }
    else{
      removed_jets++;
    }
  }
  
   // cout<<"mjets1_pt"<<m_jets[0].Pt()<<endl;
   // cout<<"mjets3_energy"<<m_jets[2].E()<<endl;
   //cout <<"mjet1_px = "<<m_jets[0].Px()<<endl;


   cout<< "m_lepton.px, py, pz, energy = "<<m_lepton[0]<<", "<<m_lepton[1]<<", "<<m_lepton[2]<<", "<<m_lepton[3]<<endl;
   cout<< "m_met.px, py, pz, energy = "<<m_met[0]<<", "<<m_met[1]<<", "<<m_met[2]<<", "<<m_met[3]<<endl;
   cout<< "m_jet1.px, py, pz, energy = "<<m_jets[0][0]<<", "<<m_jets[0][1]<<", "<<m_jets[0][2]<<", "<<m_jets[0][3]<<endl; 
   cout<< "m_jet2.px, py, pz, energy = "<<m_jets[1][0]<<", "<<m_jets[1][1]<<", "<<m_jets[1][2]<<", "<<m_jets[1][3]<<endl;
   cout<< "m_jet3.px, py, pz, energy = "<<m_jets[2][0]<<", "<<m_jets[2][1]<<", "<<m_jets[2][2]<<", "<<m_jets[2][3]<<endl;
   cout<< "m_jet4.px, py, pz, energy = "<<m_jets[3][0]<<", "<<m_jets[3][1]<<", "<<m_jets[3][2]<<", "<<m_jets[3][3]<<endl;

   // cout<<"IM stick inside jets"<<endl;

   //  cout<<"m_jets size = "<<m_jets.size()<<endl;
	double nu_px = m_met.Px();
	double nu_py = m_met.Py();

	//set all OK flags to FALSE;
	_htOK = false;
	_evtTopoOK = false;
	_ktOK = false;
	_mtOK = false;

	//
	// calculate neutrino lorentz vector (from Tobi's TopSvtAnalysis)
	//
	double nu_pz = 0.;
	double nu_e  = sqrt(pow(nu_px,2)+pow(nu_py,2));

	double Mw    = WMassPdg;
	double l_px  = m_lepton.Px();
	double l_py  = m_lepton.Py();
	double l_pz  = m_lepton.Pz();
	double l_pt  = m_lepton.Pt();
	double l_e   = m_lepton.E();
	double Mt    = sqrt(pow(l_pt+nu_e ,2)-
	               pow(l_px+nu_px,2)-
	               pow(l_py+nu_py,2));

	double A;
	// FIXME: do we need this Mt to Mw fix?
	if (Mt<Mw) A = pow(Mw,2)/2.;
	else       {           // assume Mt=Mw, rescale MET accordingly
		A = pow(Mt,2)/2.;
		double k = nu_e*l_pt - nu_px*l_px - nu_py*l_py;
		k = (k == 0. ? 0.00001 : k);
		double scf = 0.5*pow(Mw,2)/k ;
		nu_px *= scf;
		nu_py *= scf;
		nu_e = sqrt(pow(nu_px,2)+pow(nu_py,2));
	}    
  	double B = nu_px*l_px + nu_py*l_py;
	double C = TMath::Max(1. + pow(nu_e,2) * (pow(l_pz,2)-pow(l_e,2)) / pow(A+B,2)  , 0.);
	C = sqrt(C);
	double S1= (-(A+B)*l_pz + (A+B)*l_e*C) / (pow(l_pz,2)-pow(l_e,2));
	double S2= (-(A+B)*l_pz - (A+B)*l_e*C) / (pow(l_pz,2)-pow(l_e,2));

	// choose solution with smallest |l_pz| a la Run I
	nu_pz = fabs (S1) < fabs (S2) ? S1 : S2 ;

	//NGO: NOTE: neutrino PX, PY are not necessarily metPX, metPY any more!!!
	_neutrino.SetPxPyPzE(nu_px,nu_py,nu_pz,nu_e);

	cout<<"!!!!!!!!!!!!!!!"<<endl;
	 cout<< "im beofre variable defintion"<<endl;
	return removed_jets;
 
}








double LJetsTopoVars::aplanarity() const
{
   vector<TMBLorentzVector> objects(m_jets);
   objects.push_back(m_lepton);
   TopTopologicalVariables jetsPlusLepton(objects);
   return jetsPlusLepton.Aplanarity();
}

double LJetsTopoVars::centrality() const
{
   TopTopologicalVariables jets(m_jets);
   return jets.Centrality();
}

double LJetsTopoVars::sphericity() const
{
   vector<TMBLorentzVector> objects(m_jets);
   objects.push_back(m_lepton);
   TopTopologicalVariables jetsPlusLepton(objects);
   return jetsPlusLepton.Sphericity();
}

double LJetsTopoVars::ht() const
{
   TopTopologicalVariables jets(m_jets);
   return jets.Ht();
}

double LJetsTopoVars::htpluslepton() const
{
   vector<TMBLorentzVector> objects(m_jets);
   objects.push_back(m_lepton);
   TopTopologicalVariables jetsPlusLepton(objects);
   return jetsPlusLepton.Ht();
}

double LJetsTopoVars::methtpluslepton() const
{
   vector<TMBLorentzVector> objects(m_jets);
   objects.push_back(m_lepton);
   objects.push_back(m_met);
   TopTopologicalVariables metjetsPlusLepton(objects);
   return metjetsPlusLepton.Ht();
}

double LJetsTopoVars::h() const
{
   TopTopologicalVariables jets(m_jets);
   return jets.H();
}


double LJetsTopoVars::ktMinPrime() const
{
   TopTopologicalVariables jets(m_jets);
   float ktmin = jets.KtMin();
   float etw = m_met.Pt() + m_lepton.Pt();
   return ktmin/etw;
}

double LJetsTopoVars::dphiLepMet() const
{
   return kinem::delta_phi(m_met.Phi(), m_lepton.Phi());
}

double LJetsTopoVars::minDijetMass() const
{
   TopTopologicalVariables jets(m_jets);
   return jets.MinimumPairMass();
}

double LJetsTopoVars::maxJetEta() const 
{
	double jetEta = 0;
	for (unsigned int i=0; i<m_jets.size(); i++) {
	  if(TMath::Abs(m_jets.at(i).Eta()) > TMath::Abs(jetEta) ) jetEta = TMath::Abs(m_jets.at(i).Eta());
	}
	return jetEta;
}


double LJetsTopoVars::Et3() const 
{
	double Et3 = 0;
	for (unsigned int i=2; i<m_jets.size(); i++) {
	  Et3+=m_jets.at(i).Pt();
	}
	return Et3;
}

double LJetsTopoVars::minDijetDeltaR() const
{

  int nJet = m_jets.size();

  double dRmin = 9999.;
  double eTmin = 9999.;
  (void)eTmin; // silence warning
  for(int i=0;i<nJet-1;i++){
    for(int j=i+1;j<nJet;j++){
      double dR = m_jets[i].DeltaR(m_jets[j]);
      if(dR<dRmin){
	dRmin = dR;
	eTmin = std::min(m_jets[i].Pt(),m_jets[j].Pt());
      }
    }
  }
  if(dRmin>100.) {dRmin=0.;}
  
  return dRmin;
}


double LJetsTopoVars::Hz() {
	vector<TMBLorentzVector> objects;
	objects.assign(m_jets.begin(), m_jets.end());
	objects.push_back(m_lepton);
	objects.push_back(_neutrino);
	double pz = 0;
	for (vector<TMBLorentzVector>::iterator obj = objects.begin(); obj!=objects.end(); ++obj) pz += abs((*obj).Pz());
	return pz;
}

double LJetsTopoVars::HT2() {
	vector<TMBLorentzVector> objects;
	objects.assign(++m_jets.begin(), m_jets.end());
	TopTopologicalVariables topo(objects);
	return topo.Ht();
}

double LJetsTopoVars::HT2prime() {
	return HT2()/Hz();
}

double LJetsTopoVars::W_MT() {
	vector<TMBLorentzVector> objects;
	//objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.TransverseMass();
}

double LJetsTopoVars::W_Pt() {
	vector<TMBLorentzVector> objects;
	//objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.Pt();
}

double LJetsTopoVars::W_M() {
	vector<TMBLorentzVector> objects;
	objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	//	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.M();
}

double LJetsTopoVars::Jet1Jet2_M() {
  if(m_jets.size()>=2) {
    vector<TMBLorentzVector> objects;
    objects.push_back(m_jets.at(0));
    objects.push_back(m_jets.at(1));
    TopTopologicalVariables topo(objects);
    return topo.M();
  } else return -1;
}

double LJetsTopoVars::Jet1Jet2_Pt() {
  if(m_jets.size()>=2) {	
    vector<TMBLorentzVector> objects;
    objects.push_back(m_jets.at(0));
    objects.push_back(m_jets.at(1));
    TopTopologicalVariables topo(objects);
    return topo.Pt();
  } else return -1;
}

double LJetsTopoVars::Jet1Jet2_DeltaR() {
  if(m_jets.size()>=2) {  
    return m_jets.at(0).DeltaR(m_jets.at(1));
  } else return -1;
}

double LJetsTopoVars::Jet1Jet2_DeltaPhi() {
  if(m_jets.size()>=2) {  
    return TMath::Abs(m_jets.at(0).DeltaPhi(m_jets.at(1)));
  } else return -1;
}



double LJetsTopoVars::Jet1Jet2W_M() {
  if(m_jets.size()>=2) {
    vector<TMBLorentzVector> objects;
    objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
    //	objects.push_back(m_met);
    objects.push_back(m_lepton);
    objects.push_back(m_jets.at(0));
    objects.push_back(m_jets.at(1));
    TopTopologicalVariables topo(objects);
    return topo.M();
  } else return -1;
}

double LJetsTopoVars::Jet1Jet2W_Pt() {
  if(m_jets.size()>=2) {	
    vector<TMBLorentzVector> objects;
    objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
    //	objects.push_back(m_met);
    objects.push_back(m_lepton);
    objects.push_back(m_jets.at(0));
    objects.push_back(m_jets.at(1));
    TopTopologicalVariables topo(objects);
    return topo.Pt();
  } else return -1;
}

double LJetsTopoVars::DphiJMET() {
	return kinem::delta_phi(m_met.Phi(), m_jets.at(0).Phi());
}

double LJetsTopoVars::LeptonJet_DeltaR() {
  if(m_jets.size()>=2) {  
    return m_lepton.DeltaR(m_jets.at(0))< m_lepton.DeltaR(m_jets.at(1)) ? m_lepton.DeltaR(m_jets.at(0)) : m_lepton.DeltaR(m_jets.at(1));
  } else {
    return m_lepton.DeltaR(m_jets.at(0));
  }
}

double LJetsTopoVars::Muon_DeltaR() {
	//is this already stored in the muon somewhere?
	double DeltaR = 1e99;
	for (unsigned int i=0; i<m_jets.size(); i++) DeltaR = min(DeltaR, m_lepton.DeltaR(m_jets.at(i)));
	return DeltaR;
}


//
//_____________________________________________________________________
void LJetsTopoVars::calcHt(){
  
  //ht[ 0] = Ht 
  //ht[ 1] = Htp 
  //ht[ 2] = Htpp
  //ht[ 3] = Ht2
  //ht[ 4] = Ht2p
  //ht[ 5] = Ht2pp 
  //ht[ 6] = Ht3
  //ht[ 7] = Ht3p 
  //ht[ 8] = Ht3pp 
  //ht[ 9] = centrality
  //ht[10] = NJW;
  //ht[11] = eta_max
  //ht[12] = MdijetMin
  //ht[13] = Mtjets (definition from Jean-Roche)
  //ht[14] = sqrtsT (=Tobi's M_{T} in Note) from Tobi
  //ht[15] = MtAurelio
  //ht[16] = pZoverHT
  //ht[17] = Mevent
  //ht[18] = M123inv
  //ht[19] = Eta2Sum (Eta^2 sum)
  //ht[20] = mWrec
  //ht[21] = H = sum(jetE)

  //reset
  for(unsigned int i=0;i<_ht.size();i++) _ht[i]=0.;


  double h        = 0.;
  double hz       = 0.;
  double hx       = 0.;
  double hy       = 0.;
  double hzSigned = 0.;
  _ht[12]         =-1.;
  double mtjets   = 0.;
  TMBLorentzVector Mevent;
  int nJet = m_jets.size();
  for(int i=0;i<nJet;i++){
    hz += TMath::Abs(m_jets[i].Pz());
    hx += m_jets[i].Px();
    hy += m_jets[i].Py();
    h  += m_jets[i].E();
    _ht[0] += m_jets[i].Pt();
    Mevent += m_jets[i];
    hzSigned += m_jets[i].Pz();
    if(i>0) _ht[3] += m_jets[i].Pt();
    if(i>1) _ht[6] += m_jets[i].Pt();
    if(TMath::Abs(m_jets[i].Eta())>_ht[11] && i<4){
      _ht[11] = TMath::Abs(m_jets[i].Eta());
    }
    for(int j=i+1; j<nJet; j++){
      double mDijet = (m_jets[i]+m_jets[j]).Mag();
      if(_ht[12]<0. || mDijet<_ht[12]){ _ht[12]=mDijet; }
    }
    mtjets += 
      TMath::Power(m_jets[i].E(),2)  -
      TMath::Power(m_jets[i].Px(),2) -
      TMath::Power(m_jets[i].Py(),2);

    _ht[19] += m_jets[i].Eta()*m_jets[i].Eta();
  }
  _ht[21] = h; 

  // the "M_T"s
  if(mtjets > 0.){ _ht[13]=TMath::Sqrt(mtjets); } 
  _ht[14] = _ht[0]*_ht[0] - hx*hx - hy*hy;
  if(_ht[14]>0.) _ht[14] = TMath::Sqrt(_ht[14]);
  _ht[15] = h*h - hzSigned*hzSigned;
  if(_ht[15]>0.) _ht[15] = TMath::Sqrt(_ht[15]);

  if(_ht[0]>0.) _ht[16] = hzSigned/_ht[0];

  double hzNoLep = hz;
  hz += TMath::Abs(m_lepton.Pz());
  hz += TMath::Abs(_neutrino.Pz());
  if(hz!=0.){
    _ht[1]=_ht[0]/hz;
    _ht[4]=_ht[3]/hz;
    _ht[7]=_ht[6]/hz;
  }
  if(hzNoLep!=0.){
    _ht[2]=_ht[0]/hzNoLep;
    _ht[5]=_ht[3]/hzNoLep;
    _ht[8]=_ht[6]/hzNoLep;
  }
  if(h>0.){
    _ht[9] = _ht[0]/h;
  }

  //
  // NJW
  //
  double NJW=0;
  for(Int_t ijet=0; ijet<nJet-1; ijet++){
    double emin=55.;
    double emax=55.;
    if(m_jets[ijet  ].Pt() < 55.){emax=m_jets[ijet  ].Pt();}
    if(m_jets[ijet+1].Pt() < 55.){emin=m_jets[ijet+1].Pt();}
    NJW += 0.5*(emax*emax-emin*emin)*(ijet+1);
  }
  double elo=15.;
  if(m_jets[nJet-1].Pt() > elo){elo=m_jets[nJet-1].Pt();}
  NJW += 0.5*(elo*elo-(15.*15.))*(nJet);
  NJW /= ((55*55)-100.)/2.0;
  _ht[10] = NJW;
  

  // total event invariant mass
  Mevent += m_lepton;
  Mevent += _neutrino;
  _ht[17] = Mevent.Mag();

  // sum of dijet invariant masses for three highest jets
  // and mWrec
  if(nJet>2){
    double min=1e10;
    for(int i=0;i<2;i++){
      for(int j=i+1; j<3; j++){
	double m = (m_jets[i]+m_jets[j]).Mag();
	_ht[18] += m;
	double diff = TMath::Abs(WMassPdg-m);
	if(diff<min){
	  min = diff;
	  _ht[20] = m;
	}
      }
    }
  }

  _htOK = true;
}

//
//_____________________________________________________________________
void LJetsTopoVars::calcEvtTopo(){

  //evtTopo[0] = sphericity
  //evtTopo[1] = aplanarity
  //evtTopo[2] = aplanarity including muon

  int nJet = m_jets.size();

  // calculate tensor
  //
  double psum = 0.;
  for(int k=0;k<nJet;k++){
    psum += m_jets[k].Vect().Mag32();
  }
  
  TMatrixDSym M(3);
  for(int i=0;i<3;i++){
    for(int j=i;j<3;j++){
      M(i,j)=0.;
      for(int k=0;k<nJet;k++){
	M(i,j) += m_jets[k](i) * m_jets[k](j);
      }
      M(i,j)/=psum;
      if(i!=j){M(j,i) = M(i,j);} 
    }
  }
  
  //
  // get eigenvalues
  //
  TMatrixDSymEigen eigenMatrix(M);
  const TVectorD *eigen = &eigenMatrix.GetEigenValues();

  eigenval.ResizeTo(eigen->GetNrows());
  eigenval = *eigen;
  
  //NGO fix eigenvalues to zero if too small
  //otherwise ev might be marginally below zero!
  for(int i=0;i<3;i++){
    if(fabs(eigenval[i])<1e-10) eigenval[i]=0.;
  }

  _evtTopo[0] = (3./2.) * (eigenval[1]+eigenval[2]);
  _evtTopo[1] = (3./2.) *              eigenval[2];


  //
  // some consistency checks
  //
  if( eigenval[0]<eigenval[1] || 
      eigenval[0]<eigenval[2] ||
      eigenval[1]<eigenval[2]){
    cout << "ERROR: Eigenvals not ordered!" << endl;
    std::exit(1);
  }
  if(_evtTopo[0]<0. || _evtTopo[1]<0.){
    cout << "ERROR: SPHERICITY: " << _evtTopo[0] << endl;
    cout << "ERROR: APLANARITY: " << _evtTopo[1] << endl;
    M.Print();
  }


  //
  //
  // include muon in calculation
  //
  // ------------------------------------------------------
  std::vector<TMBLorentzVector> jetMu(m_jets);
  jetMu.push_back(m_lepton);
  nJet = jetMu.size();

  // calculate tensor
  psum = 0.;
  for(int k=0;k<nJet;k++){
    psum += jetMu[k].Vect().Mag32();
  }
  
  for(int i=0;i<3;i++){
    for(int j=i;j<3;j++){
      M(i,j)=0.;
      for(int k=0;k<nJet;k++){
	M(i,j) += jetMu[k](i) * jetMu[k](j);
      }
      M(i,j)/=psum;
      if(i!=j){M(j,i) = M(i,j);} 
    }
  }
  
  // get eigenvalues
  TMatrixDSymEigen       eigenMatrix_01(M);
  TVectorD eigenval_01 = eigenMatrix_01.GetEigenValues();
    for(int i=0;i<3;i++){
    if(fabs(eigenval_01[i])<1e-10) eigenval_01[i]=0.;
  }
  _evtTopo[2] = (3./2.) * eigenval_01[2];


  _evtTopoOK = true;
}


//
//_____________________________________________________________________
void LJetsTopoVars::calcKt()
{

  //kt[0] = Ktminp
  //kt[1] = Ktminpreduced
  //kt[2] = dRmin(jet,jet);

  int nJet = m_jets.size();

  double dRmin = 9999.;
  double eTmin = 9999.;
  for(int i=0;i<nJet-1;i++){
    for(int j=i+1;j<nJet;j++){
      double dR = m_jets[i].DeltaR(m_jets[j]);
      if(dR<dRmin){
	dRmin = dR;
	eTmin = std::min(m_jets[i].Pt(),m_jets[j].Pt());
      }
    }
  }
  if(dRmin>100.) {dRmin=0.;}
  
  _kt[0] = dRmin*eTmin/(m_lepton.Pt()+_neutrino.Pt());
  _kt[1] = dRmin*eTmin;
  _kt[2] = dRmin;

  _ktOK = true;
}


//
//_____________________________________________________________________
void LJetsTopoVars::calcMt()
{
  
  //mt[0] = dPhi(muon,MET)
  //mt[1] = mT

  //
  // NGO: which met should be used in calculating the mt??
  // the raw or the one form the neutrino pz calculation, which 
  // might be scaled???
  //
  double met = TMath::Sqrt(TMath::Power(_neutrino.Px(),2.)+
			   TMath::Power(_neutrino.Py(),2));
  
    
  _mt[0] = TMath::Abs(m_lepton.DeltaPhi(_neutrino));
  _mt[1] = TMath::Sqrt(2*m_lepton.Pt()*met*(1.-TMath::Cos(_mt[0])));
  
  _mtOK = true;
}
