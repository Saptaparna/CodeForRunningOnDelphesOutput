/*
  Singleton manager class for all calculators

   Author: Gena Kukartsev, 2012
*/



#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"



// ensure a single instance
LjmetFactory * LjmetFactory::instance = 0;


LjmetFactory::LjmetFactory():theSelector(0){
  mLegend = "[LjmetFactory]: ";
}



void LjmetFactory::SetExcludedCalcs( std::vector<std::string> vExcl ){

  mvExcludedCalcs = vExcl;

  std::vector<std::string>::const_iterator c;
  for (c=mvExcludedCalcs.begin(); c!=mvExcludedCalcs.end(); ++c){

    std::cout << mLegend << "removing " << *c << std::endl;

    if (mpCalculators.find(*c)!=mpCalculators.end()){

      delete mpCalculators[*c];
      mpCalculators.erase(*c);

    }

  }

  return;
}



int LjmetFactory::Register(BaseCalc * calc, 
			   std::string name){

  std::string _name = name;
  calc->setName(_name);


  if (calc->GetName().size()<1){
    // no name given - need to generate one
    char buf[256];
    int _calc_number = mpCalculators.size()+1;
    sprintf(buf, "calc%d", _calc_number);
    _name.append(buf);
  }


  if (mpCalculators.find(_name)!=mpCalculators.end()){
    std::cout << mLegend << "calculator "
	      << _name << " already registered, rename"
	      << std::endl;
  }
  else{
    calc->init();
    mpCalculators[_name] = calc;
  }


  return 0;
}



int LjmetFactory::Register(BaseEventSelector * object, 
			   std::string name){

  std::string _name = name;
  object->setName(_name);


  if (object->GetName().size()<1){
    // no name given - need to generate one
    char buf[256];
    int _object_number = mpSelectors.size()+1;
    sprintf(buf, "selector%d", _object_number);
    _name.append(buf);
  }


  if (mpSelectors.find(_name)!=mpSelectors.end()){
    std::cout << mLegend << "event selector "
	      << _name << " already registered, rename"
	      << std::endl;
  }
  else{
    object->init();
    mpSelectors[_name] = object;
  }


  return 0;
}



LjmetFactory::~LjmetFactory(){
}



BaseEventSelector * LjmetFactory::GetEventSelector(std::string name){
  //
  // Return pointer to registered event selector
  // Exit if not found
  //

  if (mpSelectors.find(name)==mpSelectors.end()){
    std::cout << mLegend << "event selector "
	      << name << " not registered"
	      << std::endl;
    std::exit(-1);
  }

  // cache the current event selector
  theSelector = mpSelectors[name];

  return theSelector;
}



void LjmetFactory::RunAllProducers(ExRootTreeReader * reader, Int_t entry,
				   BaseEventSelector * selector){
  //
  // Loop over all registered calculators and
  // run all producer methods (comes before selection)
  //

  for (std::map<std::string, BaseCalc * >::const_iterator iCalc = mpCalculators.begin();
       iCalc != mpCalculators.end(); ++iCalc){
    iCalc->second->ProduceEvent(reader, entry, selector);
  }

  return;
}



void LjmetFactory::RunAllCalculators(ExRootTreeReader * reader, Int_t entry,
				     BaseEventSelector * selector,
				     LjmetEventContent & ec){
  //
  // Loop over all registered calculators and compute
  // implemented variables
  //

  for (std::map<std::string, BaseCalc * >::const_iterator iCalc = mpCalculators.begin();
       iCalc != mpCalculators.end(); ++iCalc){
    iCalc->second->SetEventContent(&ec);
    iCalc->second->AnalyzeEvent(reader, entry, selector);
  }

  return;
}



void LjmetFactory::SetAllCalcConfig( std::map<std::string, edm::ParameterSet const> mPar ){
  //
  // Set each calc's parameter set, if present
  //
  for (std::map<std::string, BaseCalc * >::const_iterator iCalc = mpCalculators.begin();
       iCalc != mpCalculators.end(); ++iCalc){
    std::string _name = iCalc->second->GetName();
    if (mPar.find(_name)!=mPar.end()){
      iCalc->second->SetPSet(mPar[_name]);
    }
  }

  return;
}



void LjmetFactory::BeginJobAllCalc(){
  //
  // Run all BeginJob()s
  //
  for (std::map<std::string, BaseCalc * >::const_iterator iCalc = mpCalculators.begin();
       iCalc != mpCalculators.end(); ++iCalc){
    iCalc->second->BeginJob();
  }

  return;
}



void LjmetFactory::EndJobAllCalc(){
  //
  // Run all EndJob()s
  //
  for (std::map<std::string, BaseCalc * >::const_iterator iCalc = mpCalculators.begin();
       iCalc != mpCalculators.end(); ++iCalc){
    iCalc->second->EndJob();
  }

  return;
}



void LjmetFactory::RunEndEvent(ExRootTreeReader * reader, Int_t entry, 
			       LjmetEventContent & ec){
  
  theSelector->EndEvent(reader, entry, ec);

  return;
}



void LjmetFactory::RunBeginEvent(ExRootTreeReader * reader, Int_t entry, 
				 LjmetEventContent & ec){
  
  theSelector->BeginEvent(reader, entry, ec);

  return;
}
