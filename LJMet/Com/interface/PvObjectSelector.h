#ifndef Analysis_AnalysisFilters_interface_PvObjectSelector_h
#define Analysis_AnalysisFilters_interface_PvObjectSelector_h

#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"

#include <vector>
#include <string>

// make a selector for this selection
class PvObjectSelector : public Selector<reco::Vertex> {
public:

  enum Version_t { DATA2010, N_VERSIONS };
  enum Quality_t { OFF, Top_SelV3, Top_SelV3_data, N_QUALITY };
  
  PvObjectSelector( edm::ParameterSet const & parameters ){

    // legend for identifying output messages
    legend = "[";
    legend.append("PvObjectSelector");
    legend.append("]: ");

    std::string versionStr = parameters.getParameter<std::string>("version");
    std::string qualityStr = parameters.getParameter<std::string>("quality");
       
    Version_t version = N_VERSIONS;
    Quality_t quality = N_QUALITY;

    if ( versionStr == "DATA2010" ) {

      version = DATA2010;

      if      ( qualityStr == "OFF" )   quality = OFF;
      else if ( qualityStr == "Top_SelV3" ) quality = Top_SelV3;
      else if ( qualityStr == "Top_SelV3_data" ) quality = Top_SelV3_data;
                                 
    } else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of DATA2010, ..." << std::endl;
    }

    initialize( version, quality );

  }


  PvObjectSelector( Version_t version, Quality_t quality ) {
    initialize(version, quality);
  }
  
  void initialize( Version_t version, Quality_t quality )
  {

    version_ = version;
    quality_ = quality;

    std::cout << legend << "New selector initialized" << std::endl;
    std::cout << legend << "Version: " << getVersionString(version_) << std::endl;
    std::cout << legend << "Quality: " << getQualityString(quality_) << std::endl;

    if ( version_ == DATA2010){

      push_back("NOT FAKE");
      push_back("PV NDOF");
      push_back("PV Z");
      push_back("PV RHO");

      // set defaults
      set("NOT FAKE", true );
      set("PV NDOF", true);
      set("PV Z", true);
      set("PV RHO", true);

      if ( quality_ == Top_SelV3 ) {

	set("NOT FAKE", true );
	set("PV NDOF", 4.0);
	set("PV Z", 24.0);
	set("PV RHO", 2.0);

      }

      else if ( quality_ == Top_SelV3_data ) {

	set("NOT FAKE", true );
	set("PV NDOF", 4.0);
	set("PV Z", 24.0);
	set("PV RHO", 2.0);

      }

      else if ( quality_ == OFF ) {

	set("NOT FAKE", false );
	set("PV NDOF", false);
	set("PV Z", false);
	set("PV RHO", false);

      }

    }

    else{
      push_back("Undefined", true);
    }

    retInternal_ = getBitTemplate();
  }



  std::string getVersionString(Version_t ver){
    if (ver == DATA2010) return "DATA2010";
    else return "N_VERSIONS";
  }

  std::string getQualityString(Quality_t qual){
    if (qual == OFF) return "OFF";
    else if (qual == Top_SelV3) return "Top_SelV3";
    else if (qual == Top_SelV3_data) return "Top_SelV3_data";
    else return "N_QUALITY";
  }


  
  bool operator() ( reco::Vertex const & pv,  pat::strbitset & ret ) {

    // infinite loop trick to avoid excessive nesting of ifs
    while(1){

      if ( !pv.isFake() || ignoreCut("NOT FAKE") ) passCut(ret, "NOT FAKE");
      else break;
      
      if ( pv.ndof() >= cut("PV NDOF", double() ) || ignoreCut("PV NDOF") ){
	passCut(ret, "PV NDOF" );
      }
      else break;

      if ( fabs(pv.z()) < cut("PV Z", double()) || ignoreCut("PV Z") ){
	passCut(ret, "PV Z" );
      }
      else break;

      if ( fabs(pv.position().Rho()) < cut("PV RHO", double() ) || ignoreCut("PV RHO") ){
	passCut( ret, "PV RHO");
      }
      
      break;

    } // end of while(1) loop

    return (bool)ret;
  }
  
  //using EventSelector::operator();
  
 private:

  std::string legend;
  
  Version_t version_;
  Quality_t quality_;

};

#endif
