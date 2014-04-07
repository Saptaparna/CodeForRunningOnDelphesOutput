#ifndef LJMet_Com_interface_BaseCalc_h
#define LJMet_Com_interface_BaseCalc_h

/*
  Base class for all calculators

   Author: Gena Kukartsev, 2012
*/



#include <iostream>
#include <vector>

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "LJMet/DelphesData/interface/ExRootTreeReader.h"



class BaseEventSelector;
class LjmetEventContent;

namespace edm{
  class EventBase;
}



class BaseCalc {
  //
  // Base class for all calculators
  //


  friend class LjmetFactory;

  
 public:
  
  BaseCalc();
  virtual ~BaseCalc();
  BaseCalc(const BaseCalc &); // stop default

  std::string GetName();

  virtual int BeginJob() = 0;
  virtual int ProduceEvent(ExRootTreeReader * reader, Int_t entry, BaseEventSelector * selector);
  virtual int AnalyzeEvent(ExRootTreeReader * reader, Int_t entry, BaseEventSelector * selector);
  virtual int EndJob();

  std::string mName;
  std::string mLegend;

  // LJMET event content setters
  void SetHistogram(std::string name, int nbins, double low, double high);
  void SetHistValue(std::string name, double value);
  void SetValue(std::string name, int value);
  void SetValue(std::string name, double value);
  void SetValue(std::string name, std::vector<int> value);
  void SetValue(std::string name, std::vector<double> value);
  
  
 protected:

  edm::ParameterSet mPset;


 private:
  
  virtual void init();
  void setName(std::string name);
  void SetEventContent(LjmetEventContent * pEc);
  void SetPSet(edm::ParameterSet pset);

  LjmetEventContent * mpEc;

};



#endif
