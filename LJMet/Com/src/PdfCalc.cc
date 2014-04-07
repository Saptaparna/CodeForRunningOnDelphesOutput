/*
  Calculator for PDF weights and uncertainties

   Author: Gena Kukartsev, 2013
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
//#include "LJMet/Com/interface/LjmetPdfWeightProducer.h"



class LjmetFactory;



class PdfCalc : public BaseCalc{
  
 public:
  
  PdfCalc();
  virtual ~PdfCalc();

  virtual int BeginJob();
  virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
  virtual int EndJob();

  
 private:
  
  //LjmetPdfWeightProducer * pPdfWeights;


};



static int reg = LjmetFactory::GetInstance()->Register(new PdfCalc(), "PdfCalc");



PdfCalc::PdfCalc(){
}



PdfCalc::~PdfCalc(){
}



int PdfCalc::BeginJob(){

  // initialize
  // pPdfWeights = new LjmetPdfWeightProducer(mPset);
  // pPdfWeights->beginJob();

  return 0;
}



int PdfCalc::EndJob(){

  //delete pPdfWeights;

  return 0;
}



int PdfCalc::AnalyzeEvent(edm::EventBase const & event,
			     BaseEventSelector * selector){
  //
  // compute event variables here
  //


  // dummy weights
  std::vector<double> vWeights(45, 1.0);
  double weight_dummy = 1.0;

  std::string _pdf_name = "dummy";

  // for each PDF set, save weights, averages etc.
  SetValue("PdfWeightsVec_"+_pdf_name, vWeights);
  SetValue("PdfWeightPlus_"+_pdf_name, weight_dummy);
  SetValue("PdfWeightMinus_"+_pdf_name, weight_dummy);
  SetValue("PdfWeightAverage_"+_pdf_name, weight_dummy);


    
  return 0;
}
