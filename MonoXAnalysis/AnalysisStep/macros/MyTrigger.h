#ifndef MYTRIGGER
#define MYTRIGGER

#include "myIncludes.h"
using namespace std;

//           run        lumi      events counter
typedef map< int , map< int , map< int , int> > > MAP_RLE;
typedef map<int, vector<pair<int, int> > > JSON;

// trigger paths and steps
typedef map<TString, STEP> M_STEP;
typedef map<TString, PATH> M_PATH;

// variables and axes
typedef map<TString, TString> M_VAR_AXIS;

// Histogram storage
// < path , < step , < x-axis var , < num/den , histos > > > > 
typedef map< TString, TH1F*       > M_NUM_H;
typedef map< TString, M_NUM_H     > M_VAR_NUM_H;
typedef map< TString, M_VAR_NUM_H > M_STEP_VAR_NUM_H;
typedef map< TString, M_STEP_VAR_NUM_H > M_PATH_STEP_VAR_NUM_H;

// Efficiency storage
typedef map< TString, TEfficiency*     > M_FIT_E;
typedef map< TString, M_FIT_E          > M_VAR_FIT_E;
typedef map< TString, M_VAR_FIT_E      > M_STEP_VAR_FIT_E;
typedef map< TString, M_STEP_VAR_FIT_E > M_PATH_STEP_VAR_FIT_E;

class MyTrigger {

 public:

  // Constructor, destructor
  MyTrigger(TString resultName, TString offlineSel, TString era,
	    TString period, TString seed, TString json, TString field,
	    TString skim, TString HBHECleaning, TString binning);
  ~MyTrigger();

  // Methods //
  
  // json
  Int_t DefineJson();
  // chain
  Int_t InitVar();
  Int_t InitEvent();
  Int_t GetInput();
  Int_t SetBranches();
  // trigger paths
  Int_t DefinePaths();
  // histos
  Int_t ProdHistos();
  // efficiencies
  Int_t ProdEff();
  //Int_t FitEff();

 private:

  // arguments
  TString _resultName; 
  TString _offlineSel; 
  TString _era;
  TString _period; 
  TString _seed; 
  TString _json; 
  TString _field;
  TString _skim; 
  TString _HBHECleaning; 
  TString _binning;

  // I/O
  MAP_RLE   _mapRunLumiEvents;
  Bool_t    _applyJson;
  TString   _dirJson;
  JSON      _jsonMap;
  ofstream* _outIneff;
  TFile*    _outfile;
  TChain*   _ch;

  // trigger
  M_STEP _Steps;
  M_PATH _Paths;

  M_STEP::iterator _itSteps;
  M_PATH::iterator _itPaths;

  // histogram storing
  M_PATH_STEP_VAR_NUM_H _Histos;
  M_NUM_H::iterator               _itNumH;
  M_VAR_NUM_H::iterator           _itVarNumH;
  M_STEP_VAR_NUM_H::iterator      _itStepVarNumH;
  M_PATH_STEP_VAR_NUM_H::iterator _itPathStepVarNumH;

  // variables and axes
  M_VAR_AXIS _Axis;

  // TEfficiency storing
  M_PATH_STEP_VAR_FIT_E _Eff;
  M_FIT_E::iterator               _itFitE;
  M_VAR_FIT_E::iterator           _itVarFitE;
  M_STEP_VAR_FIT_E::iterator      _itStepVarFitE;
  M_PATH_STEP_VAR_FIT_E::iterator _itPathStepVarFitE;

  // chain //
  //
  // trigger objects
  Int_t           _trig_obj_n;
  vector<double> *_trig_obj_pt  ;
  vector<double> *_trig_obj_eta ;
  vector<double> *_trig_obj_phi ;
  vector<string> *_trig_obj_col ;
  //
  // global event info
  int _event, _run, _lumi;
  double   _xsec, _wgt, _kfact, _puwgt;
  int32_t  _puobs, _putrue; 
  //
  // object counters
  uint32_t _nvtx, _nmuons, _nelectrons, _ntaus, _ntightmuons, _ntightelectrons, _nphotons, _njets, _nbjets, _nfatjets;
  //
  // trigger bits
  uint8_t  _hltmet90, _hltmet120, _hltmetwithmu90, _hltmetwithmu120, _hltmetwithmu170, _hltmetwithmu300, _hltjetmet90, _hltjetmet120, _hltphoton165, _hltphoton175, _hltdoublemu, _hltsinglemu, _hltdoubleel, _hltsingleel;
  //
  // offline flags
  uint8_t  _flagcsctight, _flaghbhenoise, _flaghcallaser, _flagecaltrig, _flageebadsc, _flagecallaser, _flagtrkfail, _flagtrkpog, _flaghnoiseloose, _flaghnoisetight, _flaghnoisehilvl;
  //
  // jets/met
  double   _pfmet, _pfmetphi, _t1pfmet, _t1pfmetphi, _pfmupt, _pfmuphi, _mumet, _mumetphi, _phmet, _phmetphi, _t1mumet, _t1mumetphi, _t1phmet, _t1phmetphi;
  double   _signaljetpt, _signaljeteta, _signaljetphi, _signaljetbtag, _signaljetCHfrac, _signaljetNHfrac, _signaljetEMfrac, _signaljetCEMfrac, _signaljetmetdphi;
  double   _secondjetpt, _secondjeteta, _secondjetphi, _secondjetbtag, _secondjetCHfrac, _secondjetNHfrac, _secondjetEMfrac, _secondjetCEMfrac, _secondjetmetdphi;
  double   _thirdjetpt , _thirdjeteta , _thirdjetphi , _thirdjetbtag , _thirdjetCHfrac , _thirdjetNHfrac , _thirdjetEMfrac , _thirdjetCEMfrac , _thirdjetmetdphi ;
  double   _jetjetdphi, _jetmetdphimin, _incjetmetdphimin;
  //
  // leptons/photons
  int32_t _wzid, _l1id, _l2id, _i1id, _i2id, _i3id, _mu1pid, _mu2pid, _mu1id, _mu2id, _el1pid, _el2pid, _el1id, _el2id; 
  double  _mu1pt, _mu1eta, _mu1phi, _mu2pt, _mu2eta, _mu2phi;
  double  _el1pt, _el1eta, _el1phi, _el2pt, _el2eta, _el2phi, _phpt, _pheta, _phphi;
  double  _loosephpt, _loosepheta, _loosephphi, _loosephsieie, _loosephrndiso;

};

#endif