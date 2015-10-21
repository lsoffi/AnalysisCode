#ifndef XMETANALYSIS
#define XMETANALYSIS

#include "myIncludes.h"
#include "XMetProcess.h"

using namespace std;

//typedef map< TString , pair< TChain* , pair< vector<TString> , Int_t > > > MAP_PROCESS;
typedef map< TString , XMetProcess > MAP_PROCESS;

typedef map<TString, TH1F*>                M_VAR_H;
typedef map<TString, M_VAR_H >             M_CUT_VAR_H;
typedef map<TString, M_CUT_VAR_H >         M_PROC_CUT_VAR_H;
typedef map<TString, M_PROC_CUT_VAR_H >    M_SEL_PROC_CUT_VAR_H;

// variables and axes
typedef map<TString, TString> M_VAR_AXIS;

class XMetAnalysis {

 public:

  XMetAnalysis(TString tag, TString subdir);
  ~XMetAnalysis();

  // Chains
  Int_t DefineChainsRun1();
  Int_t DefineChainsAN15();

  // Analyses
  Int_t Analysis(Bool_t bProdHistos);
  Int_t AnalysisRun1(Bool_t bProdHistos);
  Int_t AnalysisAN15(Bool_t bProdHistos);
  Int_t StudyQCDKiller(TString signal);
  
  Int_t GetHistos(const UInt_t nS  , TString* select, 
		  const UInt_t nCut, TString* scanCut,
		  const UInt_t nV  , TString* var, 
		  vector<TString> locProcesses);

  Int_t DrawStackPlots(const UInt_t nS  , TString* select, 
		       const UInt_t nCut, TString* scanCut,
		       const UInt_t nV  , TString* var,
		       vector<TString> locProcesses);

  Int_t plot(TString select, 
	     const UInt_t nCut, TString *scanCut, Bool_t *scanReset,
	     const UInt_t nV,    TString *var, 
	     const UInt_t *xbins, Float_t **v_bins,
	     vector<TString> locProcesses);

  TCut  defineCut(TString select, TString region);

 private:

  TFile*   _outfile;
  ofstream* _outlog;
  //
  M_SEL_PROC_CUT_VAR_H _mapHistos;
  M_SEL_PROC_CUT_VAR_H::iterator _itSelProcCutVarH;
  M_PROC_CUT_VAR_H::iterator     _itProcCutVarH;
  M_CUT_VAR_H::iterator          _itCutVarH;
  M_VAR_H::iterator              _itVarH;
  //
  MAP_PROCESS _mapProcess;
  MAP_PROCESS::iterator _itProcess;
  //
  // variables and axes
  M_VAR_AXIS _Axis;
  M_VAR_AXIS _Title;
  //
  TString  _tag,  _pathMC, _pathData, _dirOut;
  Bool_t   _isAN15, _isRun1;
  Double_t _lumi, _rescale, _qcdScale;
  
};

#endif
