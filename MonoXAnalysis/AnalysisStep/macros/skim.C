#include <iostream>
#include <sstream>
#include <map>
#include <utility>
//
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEntryList.h"
//
#include "tdrstyle.h"

Int_t doskim(TString pathIn, vector<TString> sub, TString dirOut, TString tag);

Int_t skim(UInt_t iS)
{

  TString pathIn="/user/ndaci/Data/XMET/MonoJetTrees/V4/test/";
  //TString dirOut="/user/ndaci/Data/XMET/MonoJetTrees/V4/skim/";
  TString dirOut="/user/ndaci/Data/XMET/MonoJetTrees/V4/looserskim/";

  const UInt_t nSub=7;
  vector<TString> sub[nSub];
  TString name[nSub] = {"znn", "wln", "ttbar", "singlet", "qcd", "vv", "zll"};

  sub[0].push_back("znn100to200");
  sub[0].push_back("znn200to400");
  sub[0].push_back("znn400toinf");
  sub[0].push_back("znn50to100");
  //
  sub[1].push_back("wln");
  //sub[0].push_back("w4jets");
  // 
  sub[2].push_back("ttbar"); 
  //
  sub[3].push_back("singletbars");
  sub[3].push_back("singletbart");
  sub[3].push_back("singletbarw");
  sub[3].push_back("singlets");
  sub[3].push_back("singlett");
  sub[3].push_back("singletw");
  //
  sub[4].push_back("qcd1000to1400");
  sub[4].push_back("qcd120to170");
  sub[4].push_back("qcd1400to1800");
  sub[4].push_back("qcd170to300");
  sub[4].push_back("qcd1800toinf");
  sub[4].push_back("qcd300to470");
  sub[4].push_back("qcd470to600");
  sub[4].push_back("qcd600to800");
  sub[4].push_back("qcd800to1000");
  sub[4].push_back("qcd80to120");
  //
  sub[5].push_back("ww"); 
  sub[5].push_back("zz"); 
  sub[5].push_back("wz"); 
  //
  sub[6].push_back("zll");

  //for(UInt_t iS=0 ; iS<nSub ; iS++)
  //doskim(pathIn, sub[iS], dirOut, name[iS]);
    
  doskim(pathIn, sub[iS], dirOut, name[iS]);

  return 0;
}

Int_t doskim(TString pathIn, vector<TString> sub, TString dirOut, TString tag)
{

  // Output file
  TFile* outfile = new TFile(dirOut+"skim_"+tag+".root", "recreate");

  // Input chain
  TChain* chain = new TChain("tree/tree");

  for(UInt_t iD=0 ; iD<sub.size() ; iD++)
    chain->Add(pathIn+"/"+sub[iD]+"/*.root");
  
  // Skim input chain
  TCut trig   = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
  TCut veto   = "(nmuons == 0 && nelectrons == 0 && ntaus == 0)";
  TCut metID  = "(abs(pfmet - calomet) < 2*calomet)";
  TCut metCut = "mumet>200";
  TCut jetKine1 = "signaljetpt > 110 && abs(signaljeteta) < 2.4";
  TCut jetID1 = "(signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2)";

  //TCut select = trig*veto*metID*metCut*jetKine1*jetID1;
  TCut select = trig*veto*metID*jetKine1*jetID1;
  //TTree* outtree = chain->CopyTree(select,"",1);
  TTree* outtree = chain->CopyTree(select);

  //TTree* CopyTree(const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0)

  // Write output
  outtree->Write();
  outfile->Close();

  return 0;
}
