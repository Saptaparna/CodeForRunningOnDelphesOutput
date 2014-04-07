#include <iostream>

#include "TTree.h"

#include "TFile.h"

#include <sstream>

#include <utility>





void splitFile(std::string infile)

{

  TFile *file=new TFile((infile+".root").c_str());

  TTree *tree=(TTree*)file->Get("ljmet");

  int NEvts=tree->GetEntriesFast();
  cout << "NEvts= " << NEvts << endl;

  int NSplit = 10000;


  int fileNum = 0;



  Int_t   event;

  TBranch *b_event;

  tree->SetBranchAddress("event_number", &event, &b_event);



  b_event->GetEvent(NEvts-1);

  std::cout<<"event_number goes up to "<<event<<"\n ";

  int NFiles = (1+NEvts)/NSplit;
  cout << "NFiles= " << NFiles << endl;

  // Scale NSplit by event/NEvts since the event_number is not continuous

  //NSplit=NSplit*event/NEvts;

  // NSplit=double(NEvts)/double(NFiles);
  
  NSplit=double(event)/double(NFiles);

  cout << "NSplit= " << NSplit << endl;

  std::cout <<"Splitting "<<infile<<".root with "<<NEvts<<" events into "<<NFiles<<" files\n";

  for (fileNum=0;fileNum<NFiles;fileNum++){

    std::ostringstream i_;i_<<(fileNum*NSplit);

    std::ostringstream k_;k_<<(1+fileNum)*NSplit;

    std::ostringstream fileNum_;fileNum_<<fileNum;

    std::string cut="event_number>="+i_.str()+"&&event_number<"+k_.str();

   if (fileNum+1==NFiles){cut="event_number>="+i_.str();}

    std::cout<<"Copying new tree with "<<cut<<" ... ";

    TFile *outfile=new TFile((infile+"_"+fileNum_.str()+"_Skim.root").c_str(), "recreate");

    TTree *tree_new=(TTree*)tree->CopyTree(cut.c_str());

    std::cout<<"done."<<std::endl;



    tree_new->Write();

    outfile->Close();

  }



  return;

}


