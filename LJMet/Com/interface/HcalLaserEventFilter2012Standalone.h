#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "zlib.h"
#include <stdio.h>
#include <string.h>

class HcalLaserEventFilter2012  {
public:
  HcalLaserEventFilter2012(const std::string & eventFileName);
  ~HcalLaserEventFilter2012() {}
  
  
private:
  bool filter(int run, int lumiSection, int event);

  void readEventListFile(const std::string & eventFileName);
  void addEventString(const std::string & eventString);

 // ----------member data ---------------------------
  typedef std::vector< std::string > strVec;
  typedef std::vector< std::string >::iterator strVecI;
  
  std::vector< std::string > EventList_;  // vector of strings representing bad events, with each string in "run:LS:event" format
  bool verbose_;  // if set to true, then the run:LS:event for any event failing the cut will be printed out
  
  // Set run range of events in the BAD LASER LIST.  
  // The purpose of these values is to shorten the length of the EventList_ vector when running on only a subset of data
  int minrun_;
  int maxrun_;  // if specified (i.e., values > -1), then only events in the given range will be filtered
  int minRunInFile, maxRunInFile;

};

