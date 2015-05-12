#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
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

using namespace std;

typedef map< TString , pair< TChain* , pair< vector<TString> , Int_t > > > MAP_PROCESS;

class XMetAnalysis {

 public:

  XMetAnalysis(TString tag);
  ~XMetAnalysis();

  Int_t Analysis();
  Int_t StudyQCDKiller();
  Int_t DefineChains();

  Int_t plot(TString select, const UInt_t nV, TString* var, 
	     UInt_t* nBins, Float_t* xFirst, Float_t* xLast, 
	     Bool_t stack, Bool_t dolog, Bool_t unity,
	     vector<TString> locProcesses, vector<TString> labelProc);

  Int_t setStyle(TH1F* h, Int_t color);
  
  TCut  defineCut(TString select);

 private:

  TFile*   _outfile;
  ofstream* _outlog;
  //
  MAP_PROCESS _mapProcess;
  MAP_PROCESS::iterator _itProcess;
  //
  TString  _tag;
  TString  _path;
  Double_t _lumi;
  Double_t _rescale;
  
};
