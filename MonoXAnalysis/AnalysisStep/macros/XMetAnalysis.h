#include <iostream>
#include <sstream>
#include <map>
#include <utility>
//
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
//
#include "tdrstyle.h"

using namespace std;

typedef map< TString , pair< TChain* , pair< vector<TString> , Int_t > > > MAP_PROCESS;

class XMetAnalysis {

 public:

  XMetAnalysis(TString tag);
  ~XMetAnalysis();
  
  Int_t StudyQCDKiller();
  Int_t DefineChains();

  Int_t plot(TString select, TString variable, 
	     Int_t nBins, Int_t xFirst, Int_t xLast, 
	     Bool_t stack, Bool_t dolog, Bool_t unity,
	     vector<TString> locProcesses, vector<TString> labelProc);

  Int_t drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut);
  Int_t setStyle(TH1F* h, Int_t color);
  
  TCut  defineCut(TString select);

 private:

  TFile*  _outfile;
  //
  MAP_PROCESS _mapProcess;
  MAP_PROCESS::iterator _itProcess;
  //
  TString  _tag;
  TString  _path;
  Double_t _lumi;
  Double_t _rescale;
  
};
