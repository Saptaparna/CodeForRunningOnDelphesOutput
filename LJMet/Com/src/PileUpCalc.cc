/*
  Calculator for the pileup weights

   Author: Zaixing Mao, 2012
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"


class LjmetFactory;


class PileUpCalc : public BaseCalc{
  
public:
  
  PileUpCalc():mVerbosity(1){}
  virtual ~PileUpCalc(){}

  virtual int BeginJob(){

    // grab parameter values
    if (mPset.exists("verbosity")){
      mVerbosity = mPset.getParameter<int>("verbosity");
    }


    // Distribution used for Summer2012 MC. 
    // https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
    Double_t MCDist_Summer2012_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
            5.005E-06
    };

    // User generated Data distribution 
    Double_t DataDist_Oct2012[60] = {
      12238.2,
      32262.2,
      88488.8,
      225526,
      487946,
      2.47713e+06,
      1.4766e+07,
      4.44375e+07,
      1.0279e+08,
      1.95543e+08,
      3.33172e+08,
      5.09762e+08,
      6.45795e+08,
      7.16998e+08,
      7.62148e+08,
      8.00239e+08,
      8.22388e+08,
      8.21958e+08,
      8.0435e+08,
      7.75997e+08,
      7.41057e+08,
      7.02468e+08,
      6.61859e+08,
      6.17413e+08,
      5.64036e+08,
      4.97929e+08,
      4.20604e+08,
      3.37939e+08,
      2.56828e+08,
      1.83743e+08,
      1.23728e+08,
      7.88409e+07,
      4.78733e+07,
      2.77934e+07,
      1.53873e+07,
      8.07166e+06,
      3.98747e+06,
      1.85096e+06,
      810013,
      337106,
      135102,
      52853.3,
      20418,
      7846.44,
      3007.4,
      1148.47,
      435.702,
      163.72,
      60.8208,
      22.3317,
      8.11182,
      2.91896,
      1.0414,
      0.368285,
      0.128918,
      0.0445702,
      0.0151809,
      0.00508254,
      0.00166966,
      0.00053755
    };
        
    Double_t DataDist_2012ABC[60] = {
        12261,
        32847.4,
        66317.1,
        299630,
        563890,
        2.53741e+06,
        1.49795e+07,
        4.49597e+07,
        1.03846e+08,
        1.96613e+08,
        3.28962e+08,
        4.9026e+08,
        6.18174e+08,
        6.9113e+08,
        7.39317e+08,
        7.78391e+08,
        7.97071e+08,
        7.91204e+08,
        7.71447e+08,
        7.44401e+08,
        7.12717e+08,
        6.78044e+08,
        6.41833e+08,
        6.0152e+08,
        5.50323e+08,
        4.84419e+08,
        4.06659e+08,
        3.23886e+08,
        2.43456e+08,
        1.7201e+08,
        1.14375e+08,
        7.211e+07,
        4.35065e+07,
        2.52401e+07,
        1.40387e+07,
        7.42295e+06,
        3.69879e+06,
        1.72955e+06,
        760729,
        317511,
        127411,
        49858.3,
        19255.1,
        7394.18,
        2830.89,
        1079.43,
        408.751,
        153.269,
        56.8115,
        20.8129,
        7.54365,
        2.70876,
        0.964359,
        0.340284,
        0.118831,
        0.0409752,
        0.0139165,
        0.00464491,
        0.00152093,
        0.000488015
    };

    Double_t DataDist_2012ABC735[60]={
        11386.9235125 ,
        21331.8482562 ,
        60992.2060452 ,
        210242.004255 ,
        506239.616204 ,
        1156350.29943 ,
        8049161.02903 ,
        27648848.0216 ,
        66662340.5872 ,
        133988893.796 ,
        229691533.32 ,
        361938107.141 ,
        502273693.145 ,
        602777679.846 ,
        660218243.198 ,
        701829201.251 ,
        735635627.431 ,
        751625272.746 ,
        746806664.038 ,
        729917690.281 ,
        706604644.95 ,
        679099154.397 ,
        648752555.998 ,
        616909120.161 ,
        583012825.446 ,
        542340907.569 ,
        490073433.5 ,
        425933276.725 ,
        354100732.412 ,
        280427963.217 ,
        210615056.641 ,
        149646335.793 ,
        100837060.357 ,
        64896905.9795 ,
        40202571.3065 ,
        24062049.6887 ,
        13878092.6988 ,
        7658223.01082 ,
        4012202.87911 ,
        1986776.46766 ,
        930513.498258 ,
        414677.675145 ,
        177621.248632 ,
        74006.9630299 ,
        30333.709191 ,
        12333.4244509 ,
        4996.78541548 ,
        2018.96005364 ,
        812.29590168 ,
        324.590500057 ,
        128.526673404 ,
        50.366683947 ,
        19.533520426 ,
        7.50396813328 ,
        2.85871916828 ,
        1.08077720906 ,
        0.405437305414 ,
        0.150751681511 ,
        0.0554592641238 ,
        0.0201442667595 ,
    };

    Double_t DataDist_2012ABCD[60] = {
        12261.2,
        32854.9,
        90669,
        337108,
        619232,
        3.04778e+06,
        1.75106e+07,
        5.16043e+07,
        1.21968e+08,
        2.46907e+08,
        4.34887e+08,
        6.77153e+08,
        8.76662e+08,
        9.93381e+08,
        1.06899e+09,
        1.12598e+09,
        1.15987e+09,
        1.17009e+09,
        1.16505e+09,
        1.14872e+09,
        1.12265e+09,
        1.08939e+09,
        1.05e+09,
        1.00042e+09,
        9.33497e+08,
        8.45659e+08,
        7.40272e+08,
        6.24608e+08,
        5.06279e+08,
        3.92831e+08,
        2.91344e+08,
        2.06581e+08,
        1.39989e+08,
        9.04621e+07,
        5.55679e+07,
        3.23705e+07,
        1.78808e+07,
        9.3912e+06,
        4.71579e+06,
        2.28234e+06,
        1.07582e+06,
        500387,
        233316,
        111004,
        54799.4,
        28395.5,
        15488.2,
        8844.91,
        5236.19,
        3180.09,
        1964.04,
        1225.15,
        767.779,
        481.279,
        300.644,
        186.558,
        114.687,
        69.6938,
        41.7929,
        24.6979
    };

    Double_t DataDist_2012ABCD735[60] = {
        11387.5,
        21442.9,
        80085.1,
        243858,
        540355,
        1.41391e+06,
        9.45639e+06,
        3.20047e+07,
        7.66546e+07,
        1.62121e+08,
        2.95051e+08,
        4.8596e+08,
        7.00625e+08,
        8.58009e+08,
        9.50564e+08,
        1.015e+09,
        1.06421e+09,
        1.0937e+09,
        1.1032e+09,
        1.09962e+09,
        1.08648e+09,
        1.06484e+09,
        1.03669e+09,
        1.00354e+09,
        9.6361e+08,
        9.11483e+08,
        8.42427e+08,
        7.56475e+08,
        6.58237e+08,
        5.53792e+08,
        4.49088e+08,
        3.50098e+08,
        2.62172e+08,
        1.88645e+08,
        1.30367e+08,
        8.63586e+07,
        5.46877e+07,
        3.30388e+07,
        1.9035e+07,
        1.04783e+07,
        5.53404e+06,
        2.82157e+06,
        1.40006e+06,
        682932,
        331514,
        162442,
        81556,
        42498.7,
        23159.6,
        13203.3,
        7830.22,
        4789.18,
        2994.85,
        1900.82,
        1217.49,
        783.393,
        504.49,
        324.107,
        207.146,
        131.395
    };

    std::vector< float > DataDistOct;
    std::vector< float > DataDistABC;
    std::vector< float > DataDistABC735;
    std::vector< float > DataDistABCD;
    std::vector< float > DataDistABCD735;
    std::vector< float > MCDist;

    for( int i=0; i<60; ++i) {
        DataDistOct.push_back(DataDist_Oct2012[i]);
        DataDistABC.push_back(DataDist_2012ABC[i]);
        DataDistABC735.push_back(DataDist_2012ABC735[i]);
        DataDistABCD.push_back(DataDist_2012ABCD[i]);
        DataDistABCD735.push_back(DataDist_2012ABCD735[i]);
        MCDist.push_back(MCDist_Summer2012_S10[i]);
    }
    LumiWeightsOct_ = edm::LumiReWeighting(MCDist, DataDistOct, mVerbosity);
    LumiWeightsABC_ = edm::LumiReWeighting(MCDist, DataDistABC, mVerbosity);
    LumiWeightsABC735_ = edm::LumiReWeighting(MCDist, DataDistABC735, mVerbosity);
    LumiWeightsABCD_ = edm::LumiReWeighting(MCDist, DataDistABCD, mVerbosity);
    LumiWeightsABCD735_ = edm::LumiReWeighting(MCDist, DataDistABCD735, mVerbosity);


    return 0;
  }
    virtual int ProduceEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}

  
private:

    edm::LumiReWeighting LumiWeightsOct_;
    edm::LumiReWeighting LumiWeightsABC_;
    edm::LumiReWeighting LumiWeightsABC735_;
    edm::LumiReWeighting LumiWeightsABCD_;
    edm::LumiReWeighting LumiWeightsABCD735_;

    int mVerbosity;

};



static int reg = LjmetFactory::GetInstance()->Register(new PileUpCalc(), "PileUpCalc");



int PileUpCalc::ProduceEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){
    //
    // produce and store some new data for other modules
    //

    double _val = 2.34;
    selector->SetTestValue(_val);

    return 0;
}



int PileUpCalc::AnalyzeEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){

    //
    //_____ Pile-up ____________________________________________
    //
    // get pile up handle
    edm::InputTag puInfoSrc("addPileupInfo");
    edm::Handle<std::vector< PileupSummaryInfo > >  hvPuInfo;

    bool isMc = selector->IsMc();
    float nTrue = -1.;
    int nInteractions = -1;
    int bunchXing = -1;
    double MyWeightOct = 1;
    double MyWeightABC = 1;
    double MyWeightABC735 = 1;
    double MyWeightABCD = 1;
    double MyWeightABCD735 = 1;

    if ( isMc ){
        event.getByLabel(puInfoSrc, hvPuInfo);
        for (std::vector<PileupSummaryInfo>::const_iterator iPu=hvPuInfo->begin();
             iPu != hvPuInfo->end();
             ++iPu) {
            if ( iPu->getBunchCrossing() == 0 ){
                bunchXing = iPu->getBunchCrossing();
                nInteractions = iPu->getPU_NumInteractions();
                nTrue = iPu->getTrueNumInteractions();
                //hists["nInteractions"] -> Fill(iPu->getPU_NumInteractions());
                //hists["nTrueInteractions"] -> Fill(iPu->getTrueNumInteractions());
                MyWeightOct = LumiWeightsOct_.weight( nTrue );
                MyWeightABC = LumiWeightsABC_.weight( nTrue );
                MyWeightABC735 = LumiWeightsABC735_.weight( nTrue );
                MyWeightABCD = LumiWeightsABCD_.weight( nTrue );
                MyWeightABCD735 = LumiWeightsABCD735_.weight( nTrue );
                break;
            }
        }
    }
    SetValue("bunchXing", bunchXing);
    SetValue("nInteractions", nInteractions);
    SetValue("nTrueInteractions", nTrue);
    SetValue("weight_PU", MyWeightOct);
    SetValue("weight_PU_ABC", MyWeightABC);
    SetValue("weight_PU_ABC735", MyWeightABC735);
    SetValue("weight_PU_ABCD", MyWeightABCD);
    SetValue("weight_PU_ABCD735", MyWeightABCD735);

  return 0;
}
