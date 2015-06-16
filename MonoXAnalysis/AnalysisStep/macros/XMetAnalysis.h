#ifndef XMETANALYSIS
#define XMETANALYSIS

#include "myIncludes.h"
#include "XMetProcess.h"

using namespace std;

//typedef map< TString , pair< TChain* , pair< vector<TString> , Int_t > > > MAP_PROCESS;
typedef map< TString , XMetProcess > MAP_PROCESS;

class XMetAnalysis {

 public:

  XMetAnalysis(TString tag);
  ~XMetAnalysis();

  Int_t Analysis();
  Int_t StudyQCDKiller();
  Int_t DefineChains();

  Int_t plot(TString select, const UInt_t nV, TString* var, 
	     UInt_t* nBins, Float_t* xFirst, Float_t* xLast, 
	     vector<TString> locProcesses);

  TCut  defineCut(TString select, TString region);

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

#endif
