#include "TH1F.h"
#include "TH2F.h"



int maxEvents=-1;

double jet_pT_min=20;
double jet_eta_max=2.5;
double electron_pT_min=1;
double electron_eta_max=2.5;
double muon_pT_min=1;
double muon_eta_max=2.1;

double deltaRMax=0.5;
//--------------------------------------------------------------------------------------------------------
//Histograms

int nEtaBins=75;
float etaBins[]={-4.35, -4.175, -4, -3.825, -3.65, -3.475, -3.3, -3.125, -2.95, -2.868, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.868, 2.95, 3.125, 3.3, 3.475, 3.65, 3.825, 4, 4.175, 4.35, 4.525};
TH1F *h_jet_pT = new TH1F("jet_pt", ";jet p_{T} [GeV]", 60,0,600);
TH1F *h_jet_eta = new TH1F("jet_eta", ";jet #eta", nEtaBins, etaBins);
TH1F *h_jet_mult = new TH1F("jet_mult", ";jet multiplicity", 20, -0.5, 19.5);

TH1F *h_bJet_pT = new TH1F("bJet_pt", ";b-jet p_{T} [GeV]", 60,0,600);
TH1F *h_bJet_eta = new TH1F("bJet_eta", ";b-jet #eta", nEtaBins, etaBins);
TH1F *h_bJet_mult = new TH1F("bJet_mult", ";b-jet multiplicity", 13, -0.5, 12.5);

TH1F *h_electron_pT = new TH1F("electron_pt", ";electron p_{T} [GeV]", 30,0,300);
TH1F *h_electron_eta = new TH1F("electron_eta", ";electron #eta", nEtaBins, etaBins);
TH1F *h_electron_mult = new TH1F("electron_mult", ";electron multiplicity", 6, -0.5, 5.5);

TH1F *h_muon_pT = new TH1F("muon_pt", ";muon p_{T} [GeV]", 30,0,300);
TH1F *h_muon_eta = new TH1F("muon_eta", ";muon #eta", nEtaBins, etaBins);
TH1F *h_muon_mult = new TH1F("muon_mult", ";muon multiplicity", 6, -0.5, 5.5);

TH1F *h_MET=new TH1F("MET",";MET [GeV]",100,0,500);

/*
int nPTBins=5;                                                                                                                                                                                        
float pTBins[]={0, 0.1, 1, 10, 100, 300};                                                                                                                                                             
int nEtaBins=3;                                                                                                                                                                                       
float etaBins[]={0,1.5,2.5,4};                                                                                                                                                                        
*/
int nPTBins=20;
float pTBins[]={0, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 200};
int nEtaBins_2=38;
float etaBins_2[]={0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.868, 2.95, 3.125, 3.3, 3.475, 3.65, 3.825, 4, 4.175, 4.35, 4.525};
TH1F *h_genMuon_pT = new TH1F("genMuon_pt", ";muon p_{T} [GeV]", nPTBins, pTBins);
TH1F *h_genMuon_eta = new TH1F("genMuon_eta", ";muon |#eta|", nEtaBins_2, etaBins_2);
TH2F *h_genMuon_eta_pT = new TH2F("genMuon_eta_pt",";muon |#eta|;muon p_{T} [GeV]",nEtaBins_2, etaBins_2,nPTBins, pTBins);
TH1F *h_genMuon_mult = new TH1F("genMuon_mult", ";muon multiplicity", 6, -0.5, 5.5);
TH1F *h_TMRMuon_pT = new TH1F("TMRMuon_pt", ";muon p_{T} [GeV]", nPTBins, pTBins);
TH1F *h_TMRMuon_eta = new TH1F("TMRMuon_eta", ";muon |#eta|", nEtaBins_2, etaBins_2);
TH1F *h_TMRMuon_mult = new TH1F("TMRMuon_mult", ";muon multiplicity", 6, -0.5, 5.5);
TH2F *h_TMRMuon_eta_pT = new TH2F("TMRMuon_eta_pt",";muon |#eta|;muon p_{T} [GeV]",nEtaBins_2, etaBins_2,nPTBins, pTBins);
TH1F *h_TMRMuon_ID_pT = new TH1F("TMRMuon_ID_pt", ";muon p_{T} [GeV]", nPTBins, pTBins);
TH1F *h_TMRMuon_ID_eta = new TH1F("TMRMuon_ID_eta", ";muon |#eta|", nEtaBins_2, etaBins_2);
TH2F *h_TMRMuon_ID_eta_pT = new TH2F("TMRMuon_ID_eta_pt",";muon |#eta|;muon p_{T} [GeV]",nEtaBins_2, etaBins_2,nPTBins, pTBins);
TH1F *h_TMRMuon_ID_mult = new TH1F("TMRMuon_ID_mult", ";muon multiplicity", 6, -0.5, 5.5);
TH1F *h_TMRMuon_ID_Iso_pT = new TH1F("TMRMuon_ID_Iso_pt", ";muon p_{T} [GeV]", nPTBins, pTBins);
TH1F *h_TMRMuon_ID_Iso_eta = new TH1F("TMRMuon_ID_Iso_eta", ";muon |#eta|", nEtaBins_2, etaBins_2);
TH2F *h_TMRMuon_ID_Iso_eta_pT = new TH2F("TMRMuon_ID_Iso_eta_pt",";muon |#eta|;muon p_{T} [GeV]",nEtaBins_2, etaBins_2,nPTBins, pTBins);
TH1F *h_TMRMuon_ID_Iso_mult = new TH1F("TMRMuon_ID_Iso_mult", ";muon multiplicity", 6, -0.5, 5.5);




void writeOutput(){
  h_MET->Write();

  h_jet_mult->Write();
  h_jet_pT->Write();
  h_jet_eta->Write();

  h_bJet_mult->Write();
  h_bJet_pT->Write();
  h_bJet_eta->Write();

  h_electron_mult->Write();
  h_electron_pT->Write();
  h_electron_eta->Write();

  h_muon_mult->Write();
  h_muon_pT->Write();
  h_muon_eta->Write();

  h_genMuon_pT->Write();
  h_genMuon_eta->Write();
  h_genMuon_eta_pT->Write();
  h_genMuon_mult->Write();
  h_TMRMuon_pT->Write();
  h_TMRMuon_eta->Write();
  h_TMRMuon_eta_pT->Write();
  h_TMRMuon_mult->Write();
  h_TMRMuon_ID_pT->Write();
  h_TMRMuon_ID_eta->Write();
  h_TMRMuon_ID_eta_pT->Write();
  h_TMRMuon_ID_mult->Write();
  h_TMRMuon_ID_Iso_pT->Write();
  h_TMRMuon_ID_Iso_eta->Write();
  h_TMRMuon_ID_Iso_eta_pT->Write();
  h_TMRMuon_ID_Iso_mult->Write();

  TH1F* h_TMRMuon_pT_eff=(TH1F *)h_genMuon_pT->Clone("TMRMuon_pt_eff"); h_TMRMuon_pT_eff->Sumw2(); h_TMRMuon_pT_eff->Divide(h_TMRMuon_pT,h_genMuon_pT,1,1,"B"); h_TMRMuon_pT_eff->Write();
  TH1F* h_TMRMuon_eta_eff=(TH1F *)h_genMuon_eta->Clone("TMRMuon_eta_eff"); h_TMRMuon_eta_eff->Sumw2(); h_TMRMuon_eta_eff->Divide(h_TMRMuon_eta,h_genMuon_eta,1,1,"B"); h_TMRMuon_eta_eff->Write();
  TH1F* h_TMRMuon_eta_pT_eff=(TH1F *)h_genMuon_eta_pT->Clone("TMRMuon_eta_pt_eff"); h_TMRMuon_eta_pT_eff->Sumw2(); h_TMRMuon_eta_pT_eff->Divide(h_TMRMuon_eta_pT,h_genMuon_eta_pT,1,1,"B"); h_TMRMuon_eta_pT_eff->Write();

  TH1F* h_TMRMuon_ID_pT_eff=(TH1F *)h_genMuon_pT->Clone("TMRMuon_ID_pt_eff"); h_TMRMuon_ID_pT_eff->Sumw2(); h_TMRMuon_ID_pT_eff->Divide(h_TMRMuon_ID_pT,h_genMuon_pT,1,1,"B"); h_TMRMuon_ID_pT_eff->Write();
  TH1F* h_TMRMuon_ID_eta_eff=(TH1F *)h_genMuon_eta->Clone("TMRMuon_ID_eta_eff"); h_TMRMuon_ID_eta_eff->Sumw2(); h_TMRMuon_ID_eta_eff->Divide(h_TMRMuon_ID_eta,h_genMuon_eta,1,1,"B"); h_TMRMuon_ID_eta_eff->Write();
  TH1F* h_TMRMuon_ID_eta_pT_eff=(TH1F *)h_genMuon_eta_pT->Clone("TMRMuon_ID_eta_pt_eff"); h_TMRMuon_ID_eta_pT_eff->Sumw2(); h_TMRMuon_ID_eta_pT_eff->Divide(h_TMRMuon_ID_eta_pT,h_genMuon_eta_pT,1,1,"B"); h_TMRMuon_ID_eta_pT_eff->Write();

  TH1F* h_TMRMuon_ID_Iso_pT_eff=(TH1F *)h_genMuon_pT->Clone("TMRMuon_ID_Iso_pt_eff"); h_TMRMuon_ID_Iso_pT_eff->Sumw2(); h_TMRMuon_ID_Iso_pT_eff->Divide(h_TMRMuon_ID_Iso_pT,h_genMuon_pT,1,1,"B"); h_TMRMuon_ID_Iso_pT_eff->Write();
  TH1F* h_TMRMuon_ID_Iso_eta_eff=(TH1F *)h_genMuon_eta->Clone("TMRMuon_ID_Iso_eta_eff"); h_TMRMuon_ID_Iso_eta_eff->Sumw2(); h_TMRMuon_ID_Iso_eta_eff->Divide(h_TMRMuon_ID_Iso_eta,h_genMuon_eta,1,1,"B"); h_TMRMuon_ID_Iso_eta_eff->Write();
  TH1F* h_TMRMuon_ID_Iso_eta_pT_eff=(TH1F *)h_genMuon_eta_pT->Clone("TMRMuon_ID_Iso_eta_pt_eff"); h_TMRMuon_ID_Iso_eta_pT_eff->Sumw2(); h_TMRMuon_ID_Iso_eta_pT_eff->Divide(h_TMRMuon_ID_Iso_eta_pT,h_genMuon_eta_pT,1,1,"B"); h_TMRMuon_ID_Iso_eta_pT_eff->Write();


  return;
}
