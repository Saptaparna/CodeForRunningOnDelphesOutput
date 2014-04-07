#include "LJMet/Com/interface/LJetsTopoVarsNew.h"
#include "LJMet/Com/interface/TopTopologicalVariables.h"
#include "LJMet/Com/interface/AnglesUtil.h"
#include "LJMet/Com/interface/TopAngleUtils.h"

#include "TMatrixDSymEigen.h"

#include <iostream>
#include <stdexcept>


using namespace std;
using namespace top_cafe;


// Initiate LJetsTopoVarsNew using one lepton, one MET and
// 4 leading jets momenta. Note that if fewer than 4 jets are supplied,
// some variables are not well-defined. Every effort is made to process
// such situations correctly. Still, the user should use caution.
//
// Returns number of jets used: 4 or less. If more than 4 jets were
// supplied, only the first 4 are used, the rest are ignored.
//int LJetsTopoVarsNew::setEvent(std::vector<TLorentzVector> & jets,
//			    TLorentzVector & lepton,
//			    TLorentzVector & met,
//			    bool isMuon){

int LJetsTopoVarsNew::setEvent(const vector<std::pair<TLorentzVector,bool> > jets,
                               TLorentzVector & lepton,
                               TLorentzVector & met,
                               bool isMuon,
                               bool bestTop){

    m_jets.clear();
  
    eigenval.ResizeTo(3);
    eigenval.Zero();
  
    m_met = TMBLorentzVector(met);
    m_lepton = TMBLorentzVector(lepton);

    // loop over jets
    nJets = 0; // will return this as result

    int jetindex = 0;
    number_of_jets = 0;
    number_of_tagged_jets = 0;
    number_of_untagged_jets = 0;
    tagged_jet_highpt_index = 0;
    second_tagged_jet_highpt_index = 0;
    untagged_jet_highpt_index = 0;
    second_untagged_jet_index = 0;
  
    tagged_jet_highpt_DiscVal= -100;
    second_tagged_jet_highpt_DiscVal= -100;
    untagged_jet_highpt_DiscVal= -100;
    second_untagged_jet_DiscVal= -100;
    int tagged_jet_highpt = -1;
    int second_tagged_jet_highpt = -1;
    int untagged_jet_highpt = -1;
    int second_untagged_jet = -1;

    int cnt = -1;

    //for (std::vector<TLorentzVector>::const_iterator jet=jets.begin(); (jet!=jets.end()) && (m_jets.size()!=4); jet++){
    for (vector<std::pair<TLorentzVector,bool>>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet){

        ++cnt;

        TMBLorentzVector _j((*jet).first.Pt(),(*jet).first.Eta(),(*jet).first.Phi(),(*jet).first.Energy(),TMBLorentzVector::kPtEtaPhiE); 
        m_jets.push_back(_j);

        bool tagged=false;

        if  ((*jet).second) {
            number_of_tagged_jets++;
            tagged=true;
        }
        else number_of_untagged_jets++;

            
        if(tagged) {
            if (tagged_jet_highpt < 0) { 
                tagged_jet_highpt_index = jetindex; 
                //tagged_jet_highpt_DiscVal = (*jet)->bDiscriminator("combinedSecondaryVertexBJetTags");
                tagged_jet_highpt = 1; 
            }
            else if (second_tagged_jet_highpt < 0){ 
                second_tagged_jet_highpt_index = jetindex; 
                //second_tagged_jet_highpt_DiscVal = (*jet)->bDiscriminator("combinedSecondaryVertexBJetTags");
                second_tagged_jet_highpt =1; 
            }
        } 
        else { 
            if (untagged_jet_highpt < 0) { 
                untagged_jet_highpt_index = jetindex;
                //untagged_jet_highpt_DiscVal = (*jet)->bDiscriminator("combinedSecondaryVertexBJetTags"); 
                untagged_jet_highpt =1; 
            }
            else if ( second_untagged_jet < 0 ){ 
                second_untagged_jet_index = jetindex; 
                //second_untagged_jet_DiscVal = (*jet)->bDiscriminator("combinedSecondaryVertexBJetTags"); 
                second_untagged_jet =1;
            } 
        }

        double nu_px = m_met.Px();
        double nu_py = m_met.Py();

        //set all OK flags to FALSE;
        _htOK = false;
        _evtTopoOK = false;
        _ktOK = false;
        _mtOK = false;

        if (bestTop) {
            /****************************************************************/
            /// alternative method estimate Pz of neutrino//////////////
            /****************************************************************/
            TLorentzVector p4Nu, p4OtherNu;
            fzCalculator.SetMET(m_met);
            fzCalculator.SetLepton(m_lepton);
            if (m_isMuon) {
                fzCalculator.SetLeptonType("muon");
            } else {
                fzCalculator.SetLeptonType("electron");
            }
    
            double pzNu = fzCalculator.Calculate();
            p4Nu = TLorentzVector();
            p4OtherNu = TLorentzVector();
    
            p4Nu.SetPxPyPzE(m_met.Px(), m_met.Py(), pzNu, sqrt(m_met.Px()*m_met.Px()+m_met.Py()*m_met.Py()+pzNu*pzNu));
            double pzOtherNu = fzCalculator.getOther();
            p4OtherNu.SetPxPyPzE( m_met.Px(), m_met.Py(),pzOtherNu,sqrt(m_met.Px()*m_met.Px()+m_met.Py()*m_met.Py()+pzOtherNu*pzOtherNu));
            //std::cout<<"metx "<<m_met.Px()<<" mety "<<m_met.Py()<<" pzNu "<<pzNu<<" pzOtherNu "<<pzOtherNu<<std::endl;
            if ( fzCalculator.IsComplex() ) {
                double ptNu1 = fzCalculator.getPtneutrino(1);
                double ptNu2 = fzCalculator.getPtneutrino(2);
                TLorentzVector p4Nu1tmp;
                TLorentzVector p4Nu2tmp;
	
                p4Nu1tmp.SetPxPyPzE( ptNu1*m_met.Px()/m_met.Pt(), ptNu1*m_met.Py()/m_met.Pt(), pzNu, sqrt(ptNu1*ptNu1+pzNu*pzNu));
                p4Nu2tmp.SetPxPyPzE( ptNu2*m_met.Px()/m_met.Pt(), ptNu2*m_met.Py()/m_met.Pt(), pzNu, sqrt(ptNu2*ptNu2+pzNu*pzNu));
	
                TLorentzVector Wtmp;
                Wtmp = m_lepton + p4Nu1tmp;
                double Wm1 = 0;
                double Wm2 = 0;
                Wm1 = Wtmp.M();
                Wtmp = m_lepton + p4Nu2tmp;
                Wm2 = Wtmp.M();
                if ( fabs( Wm1 - 80.4) < fabs( Wm2 - 80.4) ) p4Nu = p4Nu1tmp;
                else p4Nu = p4Nu2tmp;
	
                p4OtherNu = p4Nu; // since we chose the real part, the two solutions are the same.
            }
        
            TLorentzVector p4LepW = m_lepton + p4Nu;
            TLorentzVector p4OtherLepW = m_lepton + p4OtherNu;
        
            TLorentzVector Top1;
            TLorentzVector Top2;
            double TopMass1=0.0;
            double TopMass2=0.0;
            double BestTopMass1 = -9999.0;
            double BestTopMass2 = -9999.0;
            for (unsigned int i=0; i< m_jets.size(); i++ ) {
                Top1 = p4LepW + m_jets[i];
                Top2 = p4OtherLepW + m_jets[i];
    
                //MakeTop(W,Jets[i], Top);
                TopMass1 = Top1.M();
                TopMass2 = Top2.M();
                if ( fabs(172.5-TopMass1) <  fabs(172.5-BestTopMass1) ) BestTopMass1 = TopMass1;
                if ( fabs(172.5-TopMass2) <  fabs(172.5-BestTopMass2) ) BestTopMass2 = TopMass2;

            } // loop over jets
    
            if (fabs(172.5-BestTopMass1) < fabs(172.5-BestTopMass2)) _neutrino.SetPxPyPzE(p4Nu.Px(),p4Nu.Py(),p4Nu.Pz(),p4Nu.E());    
            else _neutrino.SetPxPyPzE(p4OtherNu.Px(),p4OtherNu.Py(),p4OtherNu.Pz(),p4OtherNu.E());
        }

        else {
            double nu_pz = 0.;
            double nu_e  = sqrt(pow(nu_px,2)+pow(nu_py,2));
    
            double Mw    = 80.4;  // NGO fix this!(read from one place)
            double l_px  = m_lepton.Px();
            double l_py  = m_lepton.Py();
            double l_pz  = m_lepton.Pz();
            double l_pt  = m_lepton.Pt();
            double l_e   = m_lepton.E();
            double Mt    = sqrt(pow(l_pt+nu_e ,2)-
                                pow(l_px+nu_px,2)-
                                pow(l_py+nu_py,2));
            //if(debug) cout << "making W Mt = " << Mt << endl;
            double A;
            if (Mt<Mw) A = pow(Mw,2)/2.;
            else       {           // assume Mt=Mw, rescale MET accordingly (NGO???)
                A = pow(Mt,2)/2.;
                double k = nu_e*l_pt - nu_px*l_px - nu_py*l_py;
                k = (k == 0. ? 0.00001 : k);
                //if(debug) cout << "K = " << k << endl;
                double scf = 0.5*pow(Mw,2)/k ;
                //if(debug) cout << "neutrio px py before  = " << nu_px  << " " << nu_py <<  endl;
                nu_px *= scf;
                nu_py *= scf;
                //if(debug) cout << "neutrio px py after  = " << nu_px  << " " << nu_py <<  endl;
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
        }

        ++nJets;


    }

    return nJets;
}


int LJetsTopoVarsNew::setEventMetFixed(TLorentzVector& Jet1, TLorentzVector& Jet2, TLorentzVector& Jet3, TLorentzVector& Jet4, TLorentzVector& NewMet,TLorentzVector& Muon1, double min_dr_jet_lepton)
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
        //cout << "LJetsTopoVarsNew::setEvent(): jet pt() = " << jet -> pt() << endl;
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
  
    //cout<<"mjets1_pt"<<m_jets[0].Pt()<<endl;
    //cout<<"mjets3_energy"<<m_jets[2].E()<<endl;
    //cout <<"mjet1_px = "<<m_jets[0].Px()<<endl;


    //cout<< "m_lepton.px, py, pz, energy = "<<m_lepton[0]<<", "<<m_lepton[1]<<", "<<m_lepton[2]<<", "<<m_lepton[3]<<endl;
    //cout<< "m_met.px, py, pz, energy = "<<m_met[0]<<", "<<m_met[1]<<", "<<m_met[2]<<", "<<m_met[3]<<endl;
    //cout<< "m_jet1.px, py, pz, energy = "<<m_jets[0][0]<<", "<<m_jets[0][1]<<", "<<m_jets[0][2]<<", "<<m_jets[0][3]<<endl; 
    //cout<< "m_jet2.px, py, pz, energy = "<<m_jets[1][0]<<", "<<m_jets[1][1]<<", "<<m_jets[1][2]<<", "<<m_jets[1][3]<<endl;
    //cout<< "m_jet3.px, py, pz, energy = "<<m_jets[2][0]<<", "<<m_jets[2][1]<<", "<<m_jets[2][2]<<", "<<m_jets[2][3]<<endl;
    //cout<< "m_jet4.px, py, pz, energy = "<<m_jets[3][0]<<", "<<m_jets[3][1]<<", "<<m_jets[3][2]<<", "<<m_jets[3][3]<<endl;

    //cout<<"IM stick inside jets"<<endl;

    //cout<<"m_jets size = "<<m_jets.size()<<endl;
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


double LJetsTopoVarsNew::aplanarity() const
{
    vector<TMBLorentzVector> objects(m_jets);
    objects.push_back(m_lepton);
    TopTopologicalVariables jetsPlusLepton(objects);
    return jetsPlusLepton.Aplanarity();
}

double LJetsTopoVarsNew::centrality() const
{
    TopTopologicalVariables jets(m_jets);
    return jets.Centrality();
}

double LJetsTopoVarsNew::sphericity() const
{
    vector<TMBLorentzVector> objects(m_jets);
    objects.push_back(m_lepton);
    TopTopologicalVariables jetsPlusLepton(objects);
    return jetsPlusLepton.Sphericity();
}

double LJetsTopoVarsNew::ht() const
{
    TopTopologicalVariables jets(m_jets);
    return jets.Ht();
}

double LJetsTopoVarsNew::htpluslepton() const
{
    vector<TMBLorentzVector> objects(m_jets);
    objects.push_back(m_lepton);
    TopTopologicalVariables jetsPlusLepton(objects);
    return jetsPlusLepton.Ht();
}

double LJetsTopoVarsNew::methtpluslepton() const
{
    vector<TMBLorentzVector> objects(m_jets);
    objects.push_back(m_lepton);
    objects.push_back(m_met);
    TopTopologicalVariables metjetsPlusLepton(objects);
    return metjetsPlusLepton.Ht();
}

double LJetsTopoVarsNew::h() const
{
    TopTopologicalVariables jets(m_jets);
    return jets.H();
}

double  LJetsTopoVarsNew::H_AllJets_MinusBestJet(){
    //std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
    //GoodJetsMinusBestJet=_GoodJetsMinusBestJet
    //if(debug) cout << "halljGoodJetsMinusBestJet   " << _GoodJetsMinusBestJet.size() << endl;
    TopTopologicalVariables jets(_GoodJetsMinusBestJet);
    return jets.H();
}

double LJetsTopoVarsNew::ktMinPrime() const
{
    TopTopologicalVariables jets(m_jets);
    float ktmin = jets.KtMin();
    float etw = m_met.Pt() + m_lepton.Pt();
    return ktmin/etw;
}

double LJetsTopoVarsNew::dphiLepMet() const
{
    return kinem::delta_phi(m_met.Phi(), m_lepton.Phi());
}

double LJetsTopoVarsNew::dphiLepJ1() 
{
    if (m_jets.size()>0) return kinem::delta_phi(m_lepton.Phi(), m_jets.at(0).Phi());
    else return -100;
}

double LJetsTopoVarsNew::dphiLepJ2() 
{
    if (m_jets.size()>1) return kinem::delta_phi(m_lepton.Phi(), m_jets.at(1).Phi());
    else return -100;
}
double LJetsTopoVarsNew::dphiLepJ3() 
{
    if(m_jets.size()>2) return kinem::delta_phi(m_lepton.Phi(), m_jets.at(2).Phi());
    else return -100;
}

double LJetsTopoVarsNew::dphiLepJ4() 
{
    if(m_jets.size()>3) {
        return kinem::delta_phi(m_lepton.Phi(), m_jets.at(3).Phi());
    } else return -100;
}

double LJetsTopoVarsNew::Jet1Jet2_DeltaPhi() {
    if(m_jets.size()>1) {  
        return TMath::Abs(m_jets.at(0).DeltaPhi(m_jets.at(1)));
    } else return -100;
}

double LJetsTopoVarsNew::dphiLepLeadBTagJet(){
    if (number_of_tagged_jets){
        return kinem::delta_phi(m_lepton.Phi(), m_jets.at(tagged_jet_highpt_index).Phi());
    } else return -100;
 
}

double LJetsTopoVarsNew::dphiLepSecLeadBTagJet(){
    if (number_of_tagged_jets>1){
        return kinem::delta_phi(m_lepton.Phi(), m_jets.at(second_tagged_jet_highpt_index).Phi());
    } else return -100;
 
}

double LJetsTopoVarsNew::dphiLepLightJet(){
    if (number_of_untagged_jets){
        return kinem::delta_phi(m_lepton.Phi(), m_jets.at(untagged_jet_highpt_index).Phi());
    } else return -100;
 
}


double LJetsTopoVarsNew::minDijetMass() const
{
    TopTopologicalVariables jets(m_jets);
    return jets.MinimumPairMass();
}

double LJetsTopoVarsNew::maxJetEta() const 
{
	double jetEta = 0;
	for (unsigned int i=0; i<m_jets.size(); i++) {
        if(TMath::Abs(m_jets.at(i).Eta()) > TMath::Abs(jetEta) ) jetEta = TMath::Abs(m_jets.at(i).Eta());
	}
	return jetEta;
}


double LJetsTopoVarsNew::Et3() const 
{
	double Et3 = 0;
	for (unsigned int i=2; i<m_jets.size(); i++) {
        Et3+=m_jets.at(i).Pt();
	}
	return Et3;
}

double LJetsTopoVarsNew::minDijetDeltaR() const
{

    int nJet = m_jets.size();

    double dRmin = 9999.;
    //double eTmin = 9999.;
    for(int i=0;i<nJet-1;i++){
        for(int j=i+1;j<nJet;j++){
            double dR = m_jets[i].DeltaR(m_jets[j]);
            if(dR<dRmin){
                dRmin = dR;
                //eTmin = std::min(m_jets[i].Pt(),m_jets[j].Pt());
            }
        }
    }
    if(dRmin>100.) {dRmin=-9999.;}
  
    return dRmin;
}


double LJetsTopoVarsNew::Hz() {
	vector<TMBLorentzVector> objects;
	objects.assign(m_jets.begin(), m_jets.end());
	objects.push_back(m_lepton);
	objects.push_back(_neutrino);
	double pz = 0;
	for (vector<TMBLorentzVector>::iterator obj = objects.begin(); obj!=objects.end(); ++obj) pz += abs((*obj).Pz());
	return pz;
}

double LJetsTopoVarsNew::HT2() {
    if (m_jets.size()==0) return 0.;
	vector<TMBLorentzVector> objects;
	objects.assign(++m_jets.begin(), m_jets.end());
	TopTopologicalVariables topo(objects);
	return topo.Ht();
}

double LJetsTopoVarsNew::HT2prime() {
	return HT2()/Hz();
}

double  LJetsTopoVarsNew::HT_AllJets_MinusBestJet(){
    //std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
    //GoodJetsMinusBestJet=_GoodJetsMinusBestJet
    TopTopologicalVariables jets(_GoodJetsMinusBestJet);
    return jets.Ht();
}

double  LJetsTopoVarsNew::AllJets_MinusBestJet_Pt(){
    if(m_jets.size()>1) { 
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        vector<TMBLorentzVector> objects;
        //if(debug) cout << "alljptGoodJetsMinusBestJet   " << GoodJetsMinusBestJet.size() << endl;
        for (unsigned int i=0; i<GoodJetsMinusBestJet.size(); i++) {
            objects.push_back(GoodJetsMinusBestJet.at(i));
        }
        TopTopologicalVariables topo(objects);
        return topo.Pt();
    } else return -100;
  
}

double  LJetsTopoVarsNew::J1_NotBestJet_Pt(){
    vector<TMBLorentzVector> objects;
    //if(debug) cout << "m_jets.size()    " << m_jets.size()<< endl;
    if(m_jets.size()>1) { 
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        //if(debug) cout << "j1ptGoodJetsMinusBestJet   " << GoodJetsMinusBestJet.size() << endl;
        vector<TMBLorentzVector> objects;
        objects.push_back(GoodJetsMinusBestJet.at(0));
        TopTopologicalVariables topo(objects);
        //if(debug) cout << "top.pt    " << topo.Pt() << endl;
    
        return topo.Pt();
    
    } else return -100;
  
}

double  LJetsTopoVarsNew::J1_NotBestJet_Eta(){
    if(m_jets.size()>1) { 
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        if (GoodJetsMinusBestJet.size()){
            return GoodJetsMinusBestJet.at(0).Eta();
        } else return -100;
    } else return -100;
}

double  LJetsTopoVarsNew::J1_NotBestJet_Phi(){
    if(m_jets.size()>1) { 
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        if (GoodJetsMinusBestJet.size()){
            return GoodJetsMinusBestJet.at(0).Phi();
        } else return -100;
    } else return -100;
}



double  LJetsTopoVarsNew::J2_NotBestJet_Pt(){
    vector<TMBLorentzVector> objects;
    //if(debug) cout << "m_jets.size()    " << m_jets.size()<< endl;
    if(m_jets.size()>2) { 
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        vector<TMBLorentzVector> objects;
        objects.push_back(GoodJetsMinusBestJet.at(1));
        TopTopologicalVariables topo(objects);
        return topo.Pt();
    } else return -100;
  
}

double  LJetsTopoVarsNew::J2_NotBestJet_Eta(){
    if(m_jets.size()>2) {
        std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
        GoodJetsMinusBestJet=_GoodJetsMinusBestJet;
        if (GoodJetsMinusBestJet.size()){
            return GoodJetsMinusBestJet.at(1).Eta();
        } else return -100;
    } else return -100;
}


double LJetsTopoVarsNew::W_MT() {
	vector<TMBLorentzVector> objects;
	//objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.TransverseMass();
}

double LJetsTopoVarsNew::W_Pt() {
	vector<TMBLorentzVector> objects;
	//objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.Pt();
}

double LJetsTopoVarsNew::W_M() {
	vector<TMBLorentzVector> objects;
	objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
	//	objects.push_back(m_met);
	objects.push_back(m_lepton);
	TopTopologicalVariables topo(objects);
	return topo.M();
}

double LJetsTopoVarsNew::Jet1Jet2_M() {
    if(m_jets.size()>=2) {
        vector<TMBLorentzVector> objects;
        objects.push_back(m_jets.at(0));
        objects.push_back(m_jets.at(1));
        TopTopologicalVariables topo(objects);
        return topo.M();
    } else return -1;
}

double LJetsTopoVarsNew::Jet1Jet2_Pt() {
    if(m_jets.size()>=2) {	
        vector<TMBLorentzVector> objects;
        objects.push_back(m_jets.at(0));
        objects.push_back(m_jets.at(1));
        TopTopologicalVariables topo(objects);
        return topo.Pt();
    } else return -1;
}

double LJetsTopoVarsNew::Jet1Jet2_DeltaR() {
    if(m_jets.size()>=2) {  
        return m_jets.at(0).DeltaR(m_jets.at(1));
    } else return -1;
}

double LJetsTopoVarsNew::Jet1Jet2W_M() {
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

double LJetsTopoVarsNew::Jet1Jet2W_Pt() {
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

double LJetsTopoVarsNew::DphiJMET() {
    if(m_jets.size()==0) return -9999.;	
	return kinem::delta_phi(m_met.Phi(), m_jets.at(0).Phi());
}

double LJetsTopoVarsNew::LeptonJet_DeltaR() {

    double dR = -1.;
    if (m_jets.size()>=2) {  
        dR = m_lepton.DeltaR(m_jets.at(0))< m_lepton.DeltaR(m_jets.at(1)) ? m_lepton.DeltaR(m_jets.at(0)) : m_lepton.DeltaR(m_jets.at(1));
    } else if (m_jets.size()==1) {
        dR = m_lepton.DeltaR(m_jets.at(0));
    } else if (m_jets.size()==0) {
        dR = -1.;
    }
    return dR;
}

double LJetsTopoVarsNew::Muon_DeltaR() {
	//is this already stored in the muon somewhere?
	double DeltaR = 1e99;
	for (unsigned int i=0; i<m_jets.size(); i++) DeltaR = min(DeltaR, m_lepton.DeltaR(m_jets.at(i)));
	return DeltaR;
}

double LJetsTopoVarsNew::LeadBTagJet_DiscVal(){
    if (number_of_tagged_jets){
        return tagged_jet_highpt_DiscVal;
    } else return -100;
 
}

double LJetsTopoVarsNew::SecLeadBTagJet_DiscVal(){
    if (number_of_tagged_jets>1){
        return second_tagged_jet_highpt_DiscVal;
    } else return -100;
 
}

double LJetsTopoVarsNew::Lead12BTagJet_DiscVal(){
    if (number_of_tagged_jets>1){
        return tagged_jet_highpt_DiscVal + second_tagged_jet_highpt_DiscVal;
    } else return -100;
 
}

double LJetsTopoVarsNew::LightJet_DiscVal(){
    if (number_of_untagged_jets){
        return  untagged_jet_highpt_DiscVal;
    } else return -100;
 
}

double LJetsTopoVarsNew::BestTop() {

  TMBLorentzVector Top;
  //std::cout<<"lepton pt "<<m_lepton.Pt()<<" neutrino pt "<<_neutrino.Pt()<<std::endl;
  TMBLorentzVector W = m_lepton + _neutrino;
  TMBLorentzVector BestTop;
  bool foundindex=false;
  
  double TopMass=0.0;
  double BestTopMass = -9999.0;
  //std::cout<< " Topovar calc TestBestTop njets = " <<m_jets.size() << std::endl;
  vector<TMBLorentzVector> objects;
  SetBestTop_JetIndex(-1);
  std::vector<TMBLorentzVector> GoodJetsMinusBestJet;
  SetGoodJetsMinusBestJet(GoodJetsMinusBestJet);
  for (unsigned int i=0; i< m_jets.size(); i++ ) {
    Top = W + m_jets[i];
    TopMass = Top.M();
    //std::cout << " TopMass ==" << TopMass << "  njet == " << i  << std::endl;
    if ( fabs(172.5-TopMass) <  fabs(172.5-BestTopMass) ) {
        BestTopMass = TopMass;
        BestTop = Top;
      
        TMBLorentzVector TestBestTop = GetBestTop();
        //std::cout << " Topovar == TestBestTop" << TestBestTop.M() << std::endl;
      
        SetBestTop_JetIndex(i);
        foundindex = true;
    }
  } // loop over jets
  
  if ( foundindex){
    SetBestTop(BestTop);
    for (unsigned int i=0; i<m_jets.size(); i++ )
      if (i != _BestTop_JetIndex)
	GoodJetsMinusBestJet.push_back(m_jets[i]);
    SetGoodJetsMinusBestJet(GoodJetsMinusBestJet);
  }
  else {
      std::cout << "In LjetsTopVars \n  Error: No Best Top created!\n" << std::endl;
      return -1;
  }
  return BestTopMass;

}



double  LJetsTopoVarsNew::SecBestTop(){
    double TopMass=0.0;
    TMBLorentzVector Top;
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector SecBestTop;
    unsigned int index = GetBestTop_JetIndex();
    for (unsigned int i=0; i< m_jets.size(); i++ ) {
        if (i != index){
            SecBestTop = W + m_jets[i];
            TopMass = SecBestTop.M();
            break;
        }
    }
    return TopMass;
}

double  LJetsTopoVarsNew::SecBestBTagTop(){
    double TopMass = -10.0;
    TMBLorentzVector Top;
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector SecBestBTagTop;
    unsigned int index = GetBestTop_JetIndex();
    for (unsigned int i=0; i< m_jets.size(); i++ ) {
        if (i != index){
            if (number_of_tagged_jets){
                SecBestBTagTop = W + m_jets[tagged_jet_highpt_index];
                TopMass = SecBestBTagTop.M(); 
	
            }
            if (number_of_tagged_jets>1){
                SecBestBTagTop = W + m_jets[second_tagged_jet_highpt_index];
                TopMass = SecBestBTagTop.M(); 
            }
        }
   
    } return TopMass;
} 

double LJetsTopoVarsNew::BestTop_Pt() {
  
    TMBLorentzVector Top;
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector BestTop;
    bool foundindex=false;
  
    double TopMass=0.0;
    double BestTopMass = 5000.0;
    double BestTopPt = -10.;
    vector<TMBLorentzVector> objects;
    for (unsigned int i=0; i< m_jets.size(); i++ ) {
        Top = W + m_jets[i];
        TopMass = Top.M();
        if ( fabs(172.5-TopMass) <  fabs(172.5-BestTopMass) ) {
            BestTopMass = TopMass;
            BestTop = Top;
            BestTopPt = BestTop.Pt();
            foundindex = true;
        }
    } // loop over jets
  
    if ( !foundindex){
        cout << "In LjetsTopVars \n  Error: No Best Top created!\n" << endl;
        return -1;
    }
    return BestTopPt;
}


double LJetsTopoVarsNew::BTagTopMass()
{
  
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector LeadingbtagTop;
    double mass = -10.;
    if (number_of_tagged_jets){
        LeadingbtagTop = W + m_jets[tagged_jet_highpt_index];
        mass = LeadingbtagTop.M(); 
        SetBTagTop(LeadingbtagTop);
    
    }
    return  mass;
}

double LJetsTopoVarsNew::BTagTop_Pt()
{
  
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector LeadingbtagTop;
    double pt = -10.;
    if (number_of_tagged_jets){
        LeadingbtagTop = W + m_jets[tagged_jet_highpt_index];
        pt = LeadingbtagTop.Pt(); 
    
    }
    return  pt;
}



double LJetsTopoVarsNew::SecBTagTopMass()
{
  
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector SecLeadingbtagTop;
    double mass = -10.;
    if (number_of_tagged_jets>1){
        SecLeadingbtagTop = W + m_jets[second_tagged_jet_highpt_index];
        mass = SecLeadingbtagTop.M(); 
        SetSecBTagTop(SecLeadingbtagTop);
    }
    return  mass;
}

double LJetsTopoVarsNew::SecBTagTop_Pt()
{
    TMBLorentzVector W = m_lepton + _neutrino;
    TMBLorentzVector SecLeadingbtagTop;
    double pt = -10.;
    if (number_of_tagged_jets>1){
        SecLeadingbtagTop = W + m_jets[second_tagged_jet_highpt_index];
        pt = SecLeadingbtagTop.Pt(); 
    }
    return  pt;
}

double LJetsTopoVarsNew::Jet1TagJet2TagW_M(){
    if(m_jets.size()>=2) {
        vector<TMBLorentzVector> objects;
        objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
        //	objects.push_back(m_met);
        objects.push_back(m_lepton);
   
        objects.push_back(m_jets.at(tagged_jet_highpt_index));
        objects.push_back(m_jets.at(second_tagged_jet_highpt_index));
      
        TopTopologicalVariables topo(objects);
        return topo.M();
    } else return -10;
}

double LJetsTopoVarsNew::BestJetJet2W_M() {
    if(m_jets.size()>=2) {

        vector<TMBLorentzVector> objects;
        objects.push_back(_neutrino);  
        objects.push_back(m_lepton);
        unsigned int index = GetBestTop_JetIndex();
        bool got_notbestjet = false;
        for (unsigned int i=0; i< m_jets.size(); i++ ) {
            if (i == index )objects.push_back(m_jets.at(i));
            else if (!got_notbestjet){
                objects.push_back(m_jets.at(i));
                got_notbestjet = true;
            }
        }
        TopTopologicalVariables topo(objects);
        //std::cout<<"mass index "<<index<<" mass "<<topo.M()<<std::endl;
        return topo.M();
    } else return -10;
}

double LJetsTopoVarsNew::BestJet_Pt() {
    if(m_jets.size()>0) {
        unsigned int index = GetBestTop_JetIndex();
        return m_jets.at(index).Pt();
    } else return -10;
}

double LJetsTopoVarsNew::BestJet_Eta() {
    if(m_jets.size()>0) {
        unsigned int index = GetBestTop_JetIndex();
        return m_jets.at(index).Eta();
    } else return -10;
}

double LJetsTopoVarsNew::BestJet_Phi() {
    if(m_jets.size()>0) {
        unsigned int index = GetBestTop_JetIndex();
        return m_jets.at(index).Phi();
    } else return -10;
}


double LJetsTopoVarsNew::AllJets_M() {
    if(m_jets.size()) {
        vector<TMBLorentzVector> objects;
        for (unsigned int i=0; i<m_jets.size(); i++) {
            objects.push_back(m_jets.at(i));
        }
        TopTopologicalVariables topo(objects);
        return topo.M();
    } else return -1;
}

double LJetsTopoVarsNew::AllJetsW_M() {//sqrt_shat
    if(m_jets.size()>=2) {
        vector<TMBLorentzVector> objects;
        objects.push_back(_neutrino);  //_neutrino was made with W mass constraint; use MET instead
        //	objects.push_back(m_met);
        objects.push_back(m_lepton);
        for (unsigned int i=0; i<m_jets.size(); i++) {
            objects.push_back(m_jets.at(i));
        }
        TopTopologicalVariables topo(objects);
        return topo.M();
    } else return -1;
  
}

//---------------angular variables--------------------------
/*
int LJetsTopoVarsNew::LeptonMETxy() { 
    std::cout<<"in leptonmetxy"<<std::endl;
    vector<TMBLorentzVector> LeptonMETxy;
    LeptonMETxy.clear();
    LeptonMETxy.push_back(m_lepton[0]);
    TMBLorentzVector nu;
    nu.SetXYZM(_neutrino[0],_neutrino[1],0.0,0.0);
    LeptonMETxy.push_back(nu);
    SetLeptonMETxy(LeptonMETxy);
    //std::if(debug) cout << " Topovar == _neutrino[0]" << _neutrino[0] << std::endl;
    return 0;
} 
*/

double LJetsTopoVarsNew::Cos_BestJetLepton_BestTop() {

    double Cos_BestJetLepton_Besttop=-10.0;

    vector<TMBLorentzVector> LeptonMETxy;
    LeptonMETxy.clear();
    LeptonMETxy.push_back(m_lepton[0]);
    TMBLorentzVector nu;
    nu.SetXYZM(_neutrino[0],_neutrino[1],0.0,0.0);
    LeptonMETxy.push_back(nu);
    SetLeptonMETxy(LeptonMETxy);

    TMBLorentzVector BestTop = GetBestTop();
    unsigned int index = GetBestTop_JetIndex();

    //std::cout << "  Cos_BestJetLepton_BestTop- jet index == " <<  index <<  std::endl;
    //std::cout << " Topovar Cos_BestJetLepton_Besttop ==" << Cos_BestJetLepton_Besttop << std::endl;
    //std::cout<< " mass "<<BestTop.M()<<std::endl;
    //std::cout << "  Cos_BestJetLepton_BestTop- LeptonMETxy pT == " <<  LeptonMETxy[0].Pt() <<  std::endl;
    //std::cout << "  Cos_BestJetLepton_BestTop- bestTop.M() == " <<  BestTop.M() <<  std::endl;
    
    TopAngleUtils angleutils; 
    if (BestTop.M()) Cos_BestJetLepton_Besttop = angleutils.CosAngle(m_jets.at(index), LeptonMETxy[0], BestTop);
    
    return  Cos_BestJetLepton_Besttop;
  
}
double LJetsTopoVarsNew::Cos_LightjetJetLepton_BestTop() {
   
    double Cos_LightjetJetLepton_BestTop=-10.0;

    vector<TMBLorentzVector> LeptonMETxy;
    LeptonMETxy.clear();
    LeptonMETxy.push_back(m_lepton[0]);
    TMBLorentzVector nu;
    nu.SetXYZM(_neutrino[0],_neutrino[1],0.0,0.0);
    LeptonMETxy.push_back(nu);
    SetLeptonMETxy(LeptonMETxy);

    TMBLorentzVector bestTop = GetBestTop();

    TopAngleUtils angleutils; 
    unsigned int index = untagged_jet_highpt_index;
    if (bestTop.M()&& number_of_untagged_jets) Cos_LightjetJetLepton_BestTop = angleutils.CosAngle(m_jets.at(index), LeptonMETxy[0], bestTop);
    
    return  Cos_LightjetJetLepton_BestTop;
  
}
double LJetsTopoVarsNew::Cos_LightjetJetLepton_BTagTop() {
    
    double Cos_LightjetJetLepton_BTagTop=-10.0;

    vector<TMBLorentzVector> LeptonMETxy;
    LeptonMETxy.clear();
    LeptonMETxy.push_back(m_lepton[0]);
    TMBLorentzVector nu;
    nu.SetXYZM(_neutrino[0],_neutrino[1],0.0,0.0);
    LeptonMETxy.push_back(nu);
    SetLeptonMETxy(LeptonMETxy);

    TMBLorentzVector btagTop = GetBTagTop();
    TopAngleUtils angleutils; 
    unsigned int index = untagged_jet_highpt_index;

    if (btagTop.M() && number_of_untagged_jets) Cos_LightjetJetLepton_BTagTop = angleutils.CosAngle(m_jets.at(index), LeptonMETxy[0], btagTop);
    
    return  Cos_LightjetJetLepton_BTagTop;
  
}

//
//_____________________________________________________________________
void LJetsTopoVarsNew::calcHt(){
  
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
    if (nJet>0) {
        if(m_jets[nJet-1].Pt() > elo){elo=m_jets[nJet-1].Pt();}
        NJW += 0.5*(elo*elo-(15.*15.))*(nJet);
        NJW /= ((55*55)-100.)/2.0;
    }
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
void LJetsTopoVarsNew::calcEvtTopo(){

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
void LJetsTopoVarsNew::calcKt()
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
void LJetsTopoVarsNew::calcMt()
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
