#ifndef XMETANALYSIS
#define XMETANALYSIS

#include "myIncludes.h"
#include "XMetProcess.h"

using namespace std;

//typedef map< TString , pair< TChain* , pair< vector<TString> , Int_t > > > MAP_PROCESS;
typedef map< TString , XMetProcess > MAP_PROCESS;

typedef map<TString,TH1F*> M_VAR_H;
typedef map<TString, map<TString,TH1F*> > M_CUT_VAR_H;
typedef map<TString, map<TString, map<TString,TH1F*> > > M_PROCESS_CUT_VAR_H;

class XMetAnalysis {

 public:

  XMetAnalysis(TString tag);
  ~XMetAnalysis();

  Int_t Analysis();
  Int_t StudyQCDKiller();
  Int_t DefineChains();
  
  Int_t plot(TString select, 			 
	     const UInt_t nCut, TString* scanCut, Bool_t* scanReset,
	     const UInt_t nV, TString* var, 
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
  TString  _tag,  _path;
  Double_t _lumi, _rescale, _qcdScale;
  
};

#endif
