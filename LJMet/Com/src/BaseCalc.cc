/*
  Base class for all calculators

   Author: Gena Kukartsev, 2012
*/



#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"



BaseCalc::BaseCalc():
  mName(""),
  mLegend(""){
}



BaseCalc::~BaseCalc(){

}



void BaseCalc::init(){
  //
  // private init method to be called by LjmetFactory
  // when registering the calculator
  //

  mLegend = "["+mName+"]: ";
  std::cout << mLegend << "registering "
	    << mName << std::endl;

  return;
}



std::string BaseCalc::GetName(){
  return mName;
}



int BaseCalc::ProduceEvent(ExRootTreeReader * reader, Int_t entry, BaseEventSelector * selector){
  return 0;
}



int BaseCalc::AnalyzeEvent(ExRootTreeReader * reader, Int_t entry, BaseEventSelector * selector){
  return 0;
}



int BaseCalc::EndJob(){
  return 0;
}



void BaseCalc::setName(std::string name){
  mName = name;
}



void BaseCalc::SetValue(std::string name, int value){
  std::string _name = name+"_"+mName;
  mpEc->SetValue(_name, value);
}



void BaseCalc::SetValue(std::string name, double value){
  std::string _name = name+"_"+mName;
  mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, std::vector<int> value){
  std::string _name = name+"_"+mName;
  mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, std::vector<double> value){
  std::string _name = name+"_"+mName;
  mpEc->SetValue(_name, value);
}



void BaseCalc::SetEventContent(LjmetEventContent * pEc){
  mpEc = pEc;
  return;
}



void BaseCalc::SetPSet(edm::ParameterSet pset){
  mPset = pset;
  return;
}



void BaseCalc::SetHistogram(std::string name,
				     int nbins, 
				     double low, 
				     double high){
  //
  // Declare a new histogram to be created for the module
  //

  mpEc->SetHistogram(mName, name, nbins, low, high);

  return;
}



void BaseCalc::SetHistValue(std::string name, 
				     double value){
  mpEc->SetHistValue(mName, name, value);

  return;
}
