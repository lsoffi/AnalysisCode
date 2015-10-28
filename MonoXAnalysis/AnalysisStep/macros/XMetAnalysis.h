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

typedef map<TString, TH2F*>                M_VAR_H2;
typedef map<TString, M_VAR_H2 >             M_CUT_VAR_H2;
typedef map<TString, M_CUT_VAR_H2 >         M_PROC_CUT_VAR_H2;
typedef map<TString, M_PROC_CUT_VAR_H2 >    M_SEL_PROC_CUT_VAR_H2;

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
  Int_t CheckForwardJets(Bool_t bProdHistos);
  Int_t StudyQCDKiller(TString signal);
  Int_t QCDScaleFactor(Bool_t bProdHistos);
  
  Int_t GetHistos(const UInt_t nS  , TString* select, 
		  const UInt_t nCut, TString* scanCut,
		  const UInt_t nV  , TString* var, 
		  vector<TString> locProcesses);

  Int_t DrawStackPlots(const UInt_t nS  , TString* select, 
		       const UInt_t nCut, TString* scanCut,
		       const UInt_t nV  , const UInt_t nV2D,
		       TString* var, vector<TString> locProcesses);

  Int_t Draw2DPlots(const UInt_t nS  , TString* select, 
		    const UInt_t nCut, TString* scanCut,
		    const UInt_t nV  , const UInt_t nV2D,
		    TString* var, vector<TString> locProcesses);

  Int_t plot(TString select, 
	     const UInt_t nCut, TString *scanCut, Bool_t *scanReset,
	     const UInt_t nV,    TString *var, 
	     vector<vector<Float_t>> v_bins, vector<vector<Float_t>> v_bins2,
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

  M_SEL_PROC_CUT_VAR_H2 _mapHistos2D;
  M_SEL_PROC_CUT_VAR_H2::iterator _itSelProcCutVarH2;
  M_PROC_CUT_VAR_H2::iterator     _itProcCutVarH2;
  M_CUT_VAR_H2::iterator          _itCutVarH2;
  M_VAR_H2::iterator              _itVarH2;
  //
  MAP_PROCESS _mapProcess;
  MAP_PROCESS::iterator _itProcess;
  //
  // variables and axes
  M_VAR_AXIS _AxisX;
  M_VAR_AXIS _AxisY;
  M_VAR_AXIS _Title;
  //
  TString  _tag,  _pathMC, _pathData, _dirOut;
  Bool_t   _isAN15, _isRun1, _useLO;
  Double_t _lumi, _rescale, _qcdScale;
  
};

#endif
