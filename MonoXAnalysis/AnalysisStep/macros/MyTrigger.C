#include "MyTrigger.h"

// file lists
#include "list_SingleMuon_2015D_V2.h"

using namespace std;
Bool_t DEBUG = kFALSE;
Bool_t useCutoff=kFALSE;
UInt_t cutoff=50000; // cut-off

MyTrigger::~MyTrigger()
{
}

MyTrigger::MyTrigger(TString resultName, TString offlineSel, TString era,
		     TString period, TString seed, TString json, TString field,
		     TString skim, TString HBHECleaning, TString binning)
{

  // Record arguments in members
  _resultName = resultName; 
  _offlineSel = offlineSel; 
  _era        = era;
  _period     = period; 
  _seed       = seed; 
  _json       = json; 
  _field      = field;
  _skim       = skim; 
  _binning    = binning;
  _HBHECleaning = HBHECleaning; 

  // Set TDR style for the plots
  gROOT->Reset();
  setTDRStyle();
  gROOT->ForceStyle();

  // Define paths
  DefinePaths();

  // Define log file for inefficiencies
  _outIneff = new ofstream("results/"+_resultName+"/outIneff.txt");

  // Define JSON file
  _dirJson="/user/ndaci/Data/json/13TeV/";
  _applyJson=false;
  DefineJson();

  // Define input TChain
  InitVar();
  GetInput();
  SetBranches();
}

Int_t MyTrigger::ProdHistos()
{
  cout << "- ProdHistos()" << endl;

  // Check and announce #entries
  UInt_t entries=_ch->GetEntries();
  cout << "- Start processing: " << entries << " entries." << endl;

  // Booleans
  bool printOut=false;
  bool jetID1=false;

  // array sizes
  const UInt_t nV=6; // mumet, t1mumet, pfmet, t1pfmet, signaljetpt, signaljetNHfrac
  const UInt_t nF=2; // denom, num
  UInt_t nS;

  // initialize x-axis variables
  double var[nV] = {0, 0, 0, 0, 0, 0};

  // Inefficiency checks
  UInt_t nIneff=0;
  UInt_t nEff=0;

  // Histograms //
  cout << "- Define histograms." << endl;
  TH1F* hTemp;

  // labels
  TString hname, title;
  TString nameV[nV]={"mumet","t1mumet","pfmet","t1pfmet","signaljetpt","signaljetNHfrac"};
  TString nameAxis[nV]={"Reco PFMETNoMu [GeV]",
			"Type1 PFMETNoMu [GeV]",
			"Reco PFMET [GeV]",
			"Type1 PFMET",
			"Leading PFJet p_{T}",
			"Leading PFJet NHEF"};
  for(UInt_t iV=0; iV<nV; iV++) {
    _Axis[nameV[iV]] = nameAxis[iV];
  }

  TString nameF[nF]={"denom","num"};

  // binning
  const UInt_t xbins[nV]    = {21, 21, 21, 21, 21, 27};
  //
  int   xbins_reg[nV] = {40,  40,  40,  40,  40,  50};
  float xlow_reg[ nV] = {100, 100, 100, 100, 100, 0};
  float xup_reg[  nV] = {900, 900, 900, 900, 900, 1};
  //
  float bins_met[]  = {50,  75,  100, 110, 120, 
		       130, 140, 150, 160, 170, 
		       180, 190, 200, 220, 250, 
		       300, 350, 400, 500, 650, 
		       1000};
  //
  float bins_nhef[] = {0.00, 0.02, 0.04, 0.06, 0.08, 
		       0.10, 0.12, 0.14, 0.16, 0.18,
		       0.20, 0.22, 0.24, 0.26, 0.28,
		       0.30, 0.35, 0.40, 0.45, 0.50,
		       0.55, 0.60, 0.65, 0.70, 0.80,
		       0.90, 1.00};
  //  
  float* v_xlow[nV] = {bins_met , bins_met , bins_met , bins_met , bins_met , bins_nhef};

  // Declaration
  cout << "- Declare histograms" << endl;
  for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) { // paths

    _thePath  = &(_itPaths->second);
    _namePath = _thePath->nameP;
    nS       = _thePath->nSteps;
    cout << "-- path: " << _namePath << " ; nS=" << nS << endl;
    
    for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
      cout << "--- declare #" << iS << " : ";
      _nameStep = _thePath->steps[iS].f;

      for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
	for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den

	  hname  = "h_"+nameV[iV]+"_"+nameF[iF]+"_"+_namePath+"_"+_nameStep;
	  title = nameV[iV]+" "+nameF[iF]+" "+_namePath+" "+_nameStep;
	  cout << hname << endl;

	  if(_binning=="regular") {
	    hTemp = new TH1F(hname, title, xbins_reg[iV], xlow_reg[iV], xup_reg[iV]);				    
	  }
	  else {
	    hTemp = new TH1F(hname, title, xbins[iV]-1, v_xlow[iV]);
	  }

	  setStyle(hTemp, kBlack);
	  hTemp->SetXTitle(nameAxis[iV]);

	  _Histos[_namePath][_nameStep][nameV[iV]][nameF[iF]] = hTemp;
	}
      }
    }
  }

  // START LOOP //
  cout << "- Start looping over the chain" << endl;
  UInt_t nProcessed=0;
  UInt_t nHLT90=0;

  for(UInt_t iE=0 ; iE<entries ; iE++) {

    // PRINT OUT //
    printOut = (iE%1000==0);
    if(printOut) {
      cout << "-- processing entry #" 
	   << iE << "/" << entries
	   << endl;
    }

    // INITIALIZE //
    if(DEBUG) cout << "-- InitEvent()" << endl;
    InitEvent();
    
    // GET ENTRY //
    if(DEBUG) cout << "-- GetEntry(" << iE << ")" << endl;
    _ch->GetEntry(iE);

    // json selection
    if(DEBUG) cout << "-- Apply Json" << endl;
    if(_applyJson) {
      if( ! AcceptEventByRunAndLumiSection(_run, _lumi, _jsonMap) ) {
	continue;
      }
    }

    // output json
    if(_mapRunLumiEvents[_run][_lumi][_event]==1)
      continue;
    else _mapRunLumiEvents[_run][_lumi][_event]=1;

    // event selection
    if(_era=="2015C" && _period=="25ns") {
      if(_run==254833) continue;
    }
    
    if(_offlineSel.Contains("TightMuon")) {
      if(_ntightmuons<1) continue;
    }

    if(_offlineSel.Contains("HLTMu")) {
      if(!_hltsinglemu) continue;
    }

    if(_offlineSel.Contains("HLTDiMu")) {
      if(!_hltdoublemu) continue;      
    }

    if(_offlineSel.Contains("Ana")) {
      jetID1 = _signaljetpt>110 && abs(_signaljeteta)<2.5 && _signaljetNHfrac<0.7 && _signaljetEMfrac<0.7 && _signaljetCHfrac>0.2;
      if(!jetID1) continue;
    }

    // Count entries passing offline selection
    nProcessed++ ;
    if(_hltmet90) nHLT90++ ;

    // cut-off
    if(nProcessed>=cutoff && useCutoff) break;

    // print out every 1000 events
    if(printOut) {
      cout << "-- _trig_obj_n=" << _trig_obj_n 
	//<< " _trig_obj_pt->size()=" << _trig_obj_pt->size()
	//<< " _trig_obj_eta->size()=" << _trig_obj_eta->size()
	//<< " _trig_obj_phi->size()=" << _trig_obj_phi->size()
	   << endl;
    }

    // Debug printouts
    if(DEBUG) {
      cout << "-- Run: " << _run
	   << " Lumi: "  << _lumi
	   << " Event: " << _event
	   << endl;
      //if(nProcessed>100) break; // ND debug mode: look only at 100 entries
    }

    // get x-axis variables //
    if(DEBUG) cout << "-- Get x-axis variables" << endl;
    var[0] = _mumet;
    var[1] = _t1mumet;
    var[2] = _pfmet;
    var[3] = _t1pfmet;
    var[4] = _signaljetpt;
    var[5] = _signaljetNHfrac;

    // PROCESS TRIGGER OBJECTS //
    if(DEBUG) cout << "-- initialize trigger output" << endl;
    for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) { // paths

      _thePath = &(_itPaths->second);
      nS = _thePath->nSteps;

      for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	_thePath->steps[iS].pt = _thePath->steps[iS].phi = 0;
	_thePath->steps[iS].pass = false;
      }

    }
    
    // loop: trigger objects
    if(DEBUG) cout << "-- loop: trigger objects" << endl;
    for(UInt_t iObj=0 ; iObj<(UInt_t)_trig_obj_n ; iObj++) {
     
      if(DEBUG) cout << "--- get trigger object's properties" << endl;
      _toPt  = (*_trig_obj_pt )[iObj];
      _toEta = (*_trig_obj_eta)[iObj];
      _toPhi = (*_trig_obj_phi)[iObj];
      _toCol    = (TString)(*_trig_obj_col)[iObj];

      // loop: paths
      if(DEBUG) cout << "--- loop: paths => read trigger objects" << endl;
      for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) { // paths

	if(DEBUG) cout << "---- get _thePath: ";
	_thePath  = &(_itPaths->second);
	_namePath = _thePath->namePath;
	nS = _thePath->nSteps;
	if(DEBUG) cout << _namePath << " (" << nS << " steps)" << endl
		       << "---- loop: steps" << endl;

	for(UInt_t iS=0 ; iS<nS ; iS++) {

	  if(DEBUG) cout << "----- Step: " << _thePath->steps[iS].f ;

	  if(_toCol==_thePath->steps[iS].c) {
	    if( _toPt > _thePath->steps[iS].pt ) {
	      _thePath->steps[iS].pt  = _toPt;
	      _thePath->steps[iS].phi = _toPhi;
	    }
	    if(DEBUG) cout << " => _thePath->steps[iS].pt=" 
			   << _thePath->steps[iS].pt 
			   << " ; _toPt=" << _toPt << endl;
	  } 
	  else if(DEBUG) cout << endl;
	  
	} // end loop: steps
      } // end loop: paths
    } // end loop: trigger objects

    if(DEBUG) cout << "-- loop: paths => fill trigger output logic" << endl; 
    for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) {

      if(DEBUG) cout << "--- get _thePath: ";
      _thePath  = &(_itPaths->second);
      nS        = _thePath->nSteps;
      _namePath = _thePath->nameP;
      if(DEBUG) cout << _namePath << "(" << nS << " steps)" << endl;

      if(DEBUG) cout << "--- loop: steps" << endl;
      for(UInt_t iS=0 ; iS<nS ; iS++) {
	//
	_nameColl=_thePath->steps[iS].c;
	_nameStep=_thePath->steps[iS].n;
	_fired=false;
	if(DEBUG) cout << "---- step #" << iS
		       << " nameStep="  << _nameStep
		       << " _nameColl="  << _nameColl ;
	//
	if(_nameStep=="Full") { // check trigger bit 
	  if(     _nameColl=="hltmet90")        _fired=_hltmet90; 
	  else if(_nameColl=="hltmet120")       _fired=_hltmet120;
	  else if(_nameColl=="hltjetmet90")     _fired=_hltjetmet90;
	  else if(_nameColl=="hltjetmet120")    _fired=_hltjetmet120;
	  else if(_nameColl=="hltmetwithmu90")  _fired=_hltmetwithmu90;
	  else if(_nameColl=="hltmetwithmu120") _fired=_hltmetwithmu120;
	  else if(_nameColl=="hltmetwithmu170") _fired=_hltmetwithmu170;
	  else if(_nameColl=="OR90GeV")         {
	    _fired=_hltmet90 || _hltjetmet90 || _hltmetwithmu170;
	  }
	  else _fired=false;
	}
	//
	else { // check trigger objects
	  _fired = (_thePath->steps[iS].pt>_thePath->steps[iS].T);
	}
	//
	_thePath->steps[iS].pass = _fired;
	if(DEBUG) cout << " : pt="     << _thePath->steps[iS].pt
		       << " : thresh=" << _thePath->steps[iS].T
		       << "   ==>   _fired=" << _fired << endl;
      }
    }
    
    // serial trigger
    if(DEBUG) cout << "-- loop: paths => serial trigger" << endl;
    for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) {

      _thePath = &(_itPaths->second);
      nS = _thePath->nSteps;

      for(UInt_t iS=0 ; iS<nS ; iS++) {
	if(iS==0 || iS==nS-1) {
	  _thePath->steps[iS].serial = _thePath->steps[iS].pass;
	}
	else {
	  _thePath->steps[iS].serial = _thePath->steps[iS-1].serial && _thePath->steps[iS].pass;
	}
      }
    }

    // FILL HISTOGRAMS //
    if(DEBUG) cout << "-- loop: paths => fill histograms" << endl;
    for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) { // paths

      if(DEBUG) cout << "--- get _thePath: " ;
      _thePath  = &(_itPaths->second);
      _namePath = _thePath->nameP;
      nS       = _thePath->nSteps;
      if(DEBUG) cout << _namePath << "(" << nS << " steps)" << endl;

      if(DEBUG) cout << "--- loop: steps" << endl;
      for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths

	if(DEBUG) cout << "---- step: ";
	_nameStep=_thePath->steps[iS].f;
	if(DEBUG) cout << _nameStep << endl;

	if(DEBUG) cout << "---- loop: x-axis var" << endl;
	for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables

	  // forget about energy fractions if MET<=200
	  if(nameV[iV].Contains("frac") && _mumet<=200) continue;

	  if(DEBUG) cout << "----- var: " << nameV[iV] << endl;

	  // conditional denominator (eff1)
	  // if(iS==0 || iS==nS-1) { // L1 or entire path
	  //   h[iV][0][iP][iS]->Fill(var[iV]);
	  // }
	  // else { // intermediate HLT filters
	  //   if(_serial[iP][iS-1]) // events that fire up to step iS-1
	  //     h[iV][0][iP][iS]->Fill(var[iV]);
	  // }

	  // denominator (eff2)
	  if(DEBUG) cout << "----- fill denom" << endl;
	  _Histos[_namePath][_nameStep][nameV[iV]]["denom"]->Fill(var[iV]);
	  
	  // numerator
	  if(DEBUG) cout << "----- fill num: ";
	  if(_thePath->steps[iS].serial) { // event fired step iS of path iP
	    _Histos[_namePath][_nameStep][nameV[iV]]["num"]->Fill(var[iV]);
	    if(DEBUG) cout << "done" << endl;
	  }
	  else {
	    if(DEBUG) cout << "no!" << endl;
	  }

	} // loop:nV
      } // loop:nS
    } // loop:nP
    // end fill histograms //

    // INVESTIGATE INEFFICIENCIES //
    /*
    if(pfmet>400) {
      if(!_hltmet120 || !_hltmetwithmu120 || !_hltmetwithmu170) {
	nIneff++ ;
	FillIneff();
      }
      else nEff++ ;
    }
    */

  } // end loop over chain  

  // Write out histos //
  TFile* outfile = new TFile("results/"+_resultName+"/h_"+_resultName+".root","recreate");
  outfile->cd();

  // Write histograms //
  for(_itPathStepVarNumH=_Histos.begin() ; _itPathStepVarNumH!=_Histos.end() ; _itPathStepVarNumH++)
    for(_itStepVarNumH=_itPathStepVarNumH->second.begin() ; _itStepVarNumH!=_itPathStepVarNumH->second.end() ; _itStepVarNumH++)
      for(_itVarNumH=_itStepVarNumH->second.begin() ; _itVarNumH!=_itStepVarNumH->second.end() ; _itVarNumH++)
	for(_itNumH=_itVarNumH->second.begin() ; _itNumH!=_itVarNumH->second.end() ; _itNumH++)
	  _itNumH->second->Write();

  outfile->Write();
  outfile->Close();

  cout << endl
       << "Processed: " << nProcessed << " events" << endl
       << "HLT 90GeV: " << nHLT90 << endl;

  return 0;
}

Int_t MyTrigger::GetHistos()
{
  cout << "- GetHistos(): start" << endl;

  TString filepath="results/"+_resultName+"/h_"+_resultName+".root";
  TFile* fHistos = new TFile(filepath,"read");
  fHistos->cd();

  cout << "- opened file: " << filepath << endl;

  TH1F* hTemp;
  TString hname;

  // should make these arrays members of the class instead of copy-pastes...
  UInt_t nS=0;
  const UInt_t nF=2;
  const UInt_t nV=6;
  TString nameF[nF]={"denom","num"};
  TString nameV[nV]={"mumet","t1mumet","pfmet","t1pfmet","signaljetpt","signaljetNHfrac"};

  cout << "- loop: paths => get histos" << endl;
  for(_itPaths=_Paths.begin();_itPaths!=_Paths.end();_itPaths++) { // paths

    _thePath = &(_itPaths->second);
    nS = _thePath->nSteps;
    _namePath = _thePath->nameP;
    cout << "-- " << _namePath << endl;

    if(DEBUG) cout << "-- loop: steps" << endl;
    for(UInt_t iS=0; iS<nS; iS++) {

      _nameStep=_thePath->steps[iS].f;
      if(DEBUG) cout << "--- " << _nameStep << endl;

      if(DEBUG) cout << "--- loop: var" << endl;
      for(UInt_t iV=0; iV<nV; iV++) {

	if(DEBUG) cout << "---- " << nameV[iV] << endl;
	for(UInt_t iF=0; iF<nF; iF++) {

	  if(DEBUG) cout << "----- " << nameF[iF] << endl;
	  hname  = "h_"+nameV[iV]+"_"+nameF[iF]+"_"+_namePath+"_"+_nameStep;
	  hTemp = (TH1F*)fHistos->Get(hname);
	  if(hTemp) {
	    cout << "------ got histo: " << hname << endl;
	    _Histos[_namePath][_nameStep][nameV[iV]][nameF[iF]] = hTemp;
	  }
	}
      }
    }
  }

  return 0;
}

Int_t MyTrigger::ProdEff(Bool_t print=kFALSE)
{

  cout << "- ProdEff()" << endl
       << "- using _Histos.size()=" << _Histos.size()
       << endl;

  TFile* outfile = new TFile("results/"+_resultName+"/eff_"+_resultName+".root","recreate");
  outfile->cd();

  const UInt_t nF=2;
  TString nameFunc[nF] = {"sigmoid","cb"};
  //const UInt_t nF=1;
  //TString nameFunc[nF] = {"sigmoid"};

  TH1F *hNum, *hDen;
  TEfficiency *pEff;
  M_NUM_H theMap;
  TString nameTEff;

  // Loop: paths
  for(_itPathStepVarNumH=_Histos.begin() ; _itPathStepVarNumH!=_Histos.end() ; _itPathStepVarNumH++) {

    _namePath = _itPathStepVarNumH->first;
    _namePathFull = _Paths[_namePath].namePath;
    cout << "-- " << _namePath << " : " << _namePathFull << endl;

    // loop: steps
    for(_itStepVarNumH=_itPathStepVarNumH->second.begin() ; _itStepVarNumH!=_itPathStepVarNumH->second.end() ; _itStepVarNumH++) {

      _nameStep = _itStepVarNumH->first;
      _theStep  = &(_Steps[_nameStep]);
      cout << "--- " << _nameStep << endl;

      // loop: var
      for(_itVarNumH=_itStepVarNumH->second.begin() ; _itVarNumH!=_itStepVarNumH->second.end() ; _itVarNumH++) {

	_nameVar = _itVarNumH->first;
	theMap  = _itVarNumH->second;
	hNum    = theMap["num"];
	hDen    = theMap["denom"];

	cout << "---- hNum: " << hNum->GetEntries()
	     << "---- hDen: " << hDen->GetEntries()
	     << endl;

	if(!hNum || !hDen) {
	  cout << "---- ERROR: missing histo" << endl;
	  continue;
	}

	// Produce 1 TEff per fit func
	for(UInt_t iF=0; iF<nF; iF++) {
	  if(hNum && hDen && TEfficiency::CheckConsistency(*hNum, *hDen) ) {
	    pEff = new TEfficiency(*hNum,*hDen);
	    nameTEff = "t_"+_nameVar+"_"+_namePath+"_"+_nameStep+"_"+nameFunc[iF];
	    cout << "----- produced TEff: " << nameTEff << endl;

	    pEff->SetNameTitle(nameTEff, _namePathFull+";"+_Axis[_nameVar]+";Efficiency"); //ND
	    _Eff[_namePath][_nameStep][_nameVar][nameFunc[iF]] = pEff;
	    //_Eff[_namePath][_nameStep][_nameVar][nameFunc[iF]]->Write();
	    pEff->SetDirectory(gDirectory);
	    pEff->Write();
	    pEff->SetDirectory(0);
	    if(print) DrawEff(pEff, nameTEff, _nameVar, "", _theStep->C, _theStep->S);
	  }
	} // end loop: fit func

      } // end loop: var
    } // end loop: steps
  } // end loop: paths

  outfile->cd();
  outfile->Write();
  outfile->Close();

  return 0;
}

Int_t MyTrigger::FitEff()
{

  cout << "- FitEff()" << endl
       << "- using _Eff.size()=" << _Eff.size()
       << endl;

  // Output file
  TFile* outfile = new TFile("results/"+_resultName+"/fits_"+_resultName+".root","recreate");
  outfile->cd();

  UInt_t nS=0;
  const UInt_t nF=2;
  TString nameFunc[nF]  = {"sigmoid","cb"};
  const UInt_t nPar[nF] = {3,5};
  TString nameFuncLoc, nameTEff, s_eff95;

  M_FIT_E theMap;
  double threshold, eff95;
  TEfficiency *pEff;
  TF1 *func, *fitEffTemp;
  TPaveText *pt2;

  // Loop: paths
  for(_itPathStepVarFitE=_Eff.begin() ; _itPathStepVarFitE!=_Eff.end() ; _itPathStepVarFitE++) {

    _namePath = _itPathStepVarFitE->first;
    _thePath  = &(_Paths[_namePath]);
    _namePathFull = _thePath->namePath;
    nS = _thePath->nSteps;
    cout << "-- Path: " << _namePath << " : " << _namePathFull << endl;

    // loop: steps
    for(_itStepVarFitE=_itPathStepVarFitE->second.begin() ; _itStepVarFitE!=_itPathStepVarFitE->second.end() ; _itStepVarFitE++) {

      _nameStep = _itStepVarFitE->first;
      _theStep  = &(_Steps[_nameStep]);

      // in case this step is the trigger bit, use the threshold from the last filter
      if(_nameStep[0]==TString("b")) threshold = _thePath->steps[nS>=2 ? nS-2 : 0].T;
      else                           threshold = _theStep->T;
      cout << "--- Step: " << _nameStep << " (threshold=" << threshold << ")" << endl;

      // loop: var
      for(_itVarFitE=_itStepVarFitE->second.begin() ; _itVarFitE!=_itStepVarFitE->second.end() ; _itVarFitE++) {

	_nameVar = _itVarFitE->first;
	theMap  = _itVarFitE->second;

	// loop: fit functions
	for(UInt_t iF=0; iF<nF; iF++) {

	  pEff = theMap[nameFunc[iF]];
	  nameTEff = pEff->GetName();
	  cout << "----- TEff: " << nameTEff << endl;

	  // energy fractions: store but not fit
	  if(!_nameVar.Contains("frac")) {

	    // prepare fit function
	    nameFuncLoc = "func_"+nameTEff+"_"+nameFunc[iF];
	    if(iF==0) func = new TF1(nameFuncLoc,Sigmoid,0,1000,nPar[iF]);
	    else      func = new TF1(nameFuncLoc,ErfCB,  0,1000,nPar[iF]);
	    prepareFunc(func, nameFunc[iF], threshold);
	    
	    cout << "----- Fitting: " << nameTEff 
		 << " using: " << nameFuncLoc
		 << " (threshold=" << threshold << ")" 
		 << endl;
	    
	    // perform the fit
	    pEff->Fit(func , "R");
	    
	    // work w/ the fit results
	    fitEffTemp = 
	      (TF1*)(pEff->GetListOfFunctions()->FindObject(nameFuncLoc));
	    setStyle(fitEffTemp, _theStep->C, _theStep->S);
	    eff95 = dichotomy(0.95, 0, 1000, 0.0000001, *fitEffTemp, true);
	    s_eff95 = "#epsilon = 95% @ "+TString(Form("%.0f", eff95))+" GeV";
	    
	  } // end if: var is a fraction
	  
	  // draw
	  DrawEff(pEff, "fit_"+nameTEff, _nameVar, s_eff95, _theStep->C, _theStep->S);
	  
	  // write
	  pEff->Write();

	} // end loop: fit func
      } // end loop: var
    } // end loop: steps
  } // end loop: paths

  outfile->Write();
  outfile->Close();

  return 0;
}

Int_t MyTrigger::DrawEff(TEfficiency* pEff, TString nameTEff,
			 TString nameVar  , TString s_eff95,
			 Int_t color      , Int_t style)
{

  TCanvas c("c","c",0,0,600,600);
  TPaveText pt2(0.58,0.15,0.85,0.22,"brNDC");

  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.4);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);

  pEff->SetLineColor(  color);
  pEff->SetMarkerColor(color);
  pEff->SetMarkerStyle(style);

  pEff->Draw("AP");
  
  // stat box (95% eff point)
  if(!nameVar.Contains("frac") && s_eff95!="") {
    pt2.SetLineColor(1);
    pt2.SetTextColor(1);
    pt2.SetTextFont(42);
    pt2.SetTextSize(0.03);
    pt2.SetFillColor(kWhite);
    pt2.SetShadowColor(kWhite);
    pt2.AddText(s_eff95);
    pt2.Draw();
  }
  
  // print	
  c.Print("results/"+_resultName+"/"+nameTEff+".pdf","pdf");

  return 0;
}

Int_t MyTrigger::prepareFunc(TF1* func, TString type, double threshold)
{

  if(type=="cb") {
    func->SetParName(0, "m0");
    func->SetParName(1, "sigma");
    func->SetParName(2, "alpha");
    func->SetParName(3, "n");
    func->SetParName(4, "norm");
    func->SetParameter(0, threshold);
    func->SetParameter(1, 1);
    func->SetParameter(2, 1);
    func->SetParameter(3, 5);
    func->SetParameter(4, 1);
    func->SetParLimits(1, 0.01, 50);
    func->SetParLimits(2, 0.01, 8);
    func->SetParLimits(3, 1.1, 35);
    func->SetParLimits(4, 0.6, 1);
  }
  //
  else {
    func->SetParName(0, "midpoint");
    func->SetParName(1, "steepness");
    func->SetParName(2, "max");
    func->SetParameter(0, threshold);
    func->SetParameter(1, 0.06);
    func->SetParameter(2, 1);
    func->SetParLimits(2, 0.995, 1);
    func->SetLineWidth(2);
  }

  func->SetLineWidth(2);

  return 0;
}

Int_t MyTrigger::DefinePaths()
{

  // STEPS //

  // Level-1 Seeds
  _Steps["L1ETM60"]={.T=60,.c="hltL1extraParticles:MET:HLT",.n="L1",.t="L1",.f="L1ETM60",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};               
  _Steps["L1ETM50"]={.T=50,.c="hltL1extraParticles:MET:HLT",.n="L1",.t="L1",.f="L1ETM50",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};               

  // MonoPFJet
  _Steps["Jet65"]={.T=65,.c="hltAK4CaloJetsCorrectedIDPassed::HLT",.n="Jet",  .t="Jet",  .f="Jet65",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["PFJet80"]={.T=80,.c="hltAK4PFJetsTightIDCorrected::HLT", .n="PFJet",.t="PFJet",.f="PFJet80",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["bJ80MuPFM90"]={.T=0,.c="hltjetmet90",                    .n="Full",.t="Full path",.f="bJ80MuPFM90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["bJ80MuPFM120"]={.T=0,.c="hltjetmet120",                  .n="Full",.t="Full path",.f="bJ80MuPFM120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // OR
  _Steps["bOR90GeV"] = {.T=0,.c="OR90GeV",                  .n="Full",.t="Full path",.f="bOR90GeV",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMNoMu90
  _Steps["MET65"]   ={.T=65,.c="hltMet::HLT",               .n="MET",  .t="MET",.f="MET65",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  _Steps["METC55"]  ={.T=55,.c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.f="METC55",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};      
  _Steps["METJ55"]  ={.T=55,.c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.f="METJ55",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["MHT65"]   ={.T=65,.c="hltMht::HLT",               .n="MHT",  .t="MHT",.f="MHT65",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["MuMHT90"] ={.T=90,.c="hltPFMHTNoMuTightID::HLT",  .n="PFMHT",.t="PFMHTNoMu",.f="MuMHT90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["MuMET90"] ={.T=90,.c="hltPFMETNoMuProducer::HLT", .n="PFMET",.t="PFMETNoMu",.f="MuMET90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["bMuPFM90"]={.T=0, .c="hltmet90",                  .n="Full", .t="Full path",.f="bMuPFM90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMNoMu120
  _Steps["MET80"] ={.T=80, .c="hltMet::HLT",               .n="MET",  .t="MET",.f="MET80",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METC70"]={.T=70, .c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.f="METC70",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METJ70"]={.T=70, .c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.f="METJ70",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["MHT80"]    ={.T=80, .c="hltMht::HLT",               .n="MHT",  .t="MHT",.f="MHT80",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["MuMHT120"] ={.T=120,.c="hltPFMHTNoMuTightID::HLT",  .n="PFMHT",.t="PFMHTNoMu",.f="MuMHT120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["MuMET120"] ={.T=120,.c="hltPFMETNoMuProducer::HLT", .n="PFMET",.t="PFMETNoMu",.f="MuMET120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["bMuPFM120"]={.T=0,  .c="hltmet120",                 .n="Full", .t="Full path",.f="bMuPFM120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFM90
  _Steps["MET70"]  ={.T=70,.c="hltMet::HLT",          .n="MET",  .t="MET",.f="MET70",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  _Steps["MHT70"]  ={.T=70,.c="hltMht::HLT",          .n="MHT",  .t="MHT",.f="MHT70",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["PFMHT90"]={.T=90,.c="hltPFMHTTightID::HLT", .n="PFMHT",.t="PFMHT",.f="PFMHT90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["PFMET90"]={.T=90,.c="hltPFMETProducer::HLT",.n="PFMET",.t="PFMET",.f="PFMET90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["bPFM90"] ={.T=0, .c="hltmetwithmu90",       .n="Full", .t="Full path",.f="bPFM90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFM120
  _Steps["MET90"]   ={.T=90,.c="hltMet::HLT",           .n="MET",  .t="MET",.f="MET90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  _Steps["MHT90"]   ={.T=90,.c="hltMht::HLT",           .n="MHT",  .t="MHT",.f="MHT90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["PFMHT120"]={.T=120,.c="hltPFMHTTightID::HLT", .n="PFMHT",.t="PFMHT",.f="PFMHT120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["PFMET120"]={.T=120,.c="hltPFMETProducer::HLT",.n="PFMET",.t="PFMET",.f="PFMET120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  _Steps["bPFM120"] ={.T=0, .c="hltmetwithmu120",       .n="Full", .t="Full path",.f="bPFM120",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMET170
  //_Steps["MET90"]={.T=90, .c="hltMet::HLT",               .n="MET",  .t="MET",.f="MET90",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METC80"]={.T=80, .c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.f="METC80",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METJ80"]={.T=80, .c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.f="METJ80",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["PFMET170"]={.T=170,.c="hltPFMETProducer::HLT",   .n="PFMET",.t="PFMET",.f="PFMET170",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["bMET170"] ={.T=0,  .c="hltmetwithmu170",         .n="Full", .t="Full path",.f="bMET170",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // CaloMET200
  _Steps["MET210"] ={.T=210,.c="hltMet::HLT",               .n="MET", .t="MET",.f="MET210",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METC200"]={.T=200,.c="hltMetClean::HLT",          .n="METC",.t="Cleaned MET",.f="METC200",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  _Steps["METJ200"]={.T=200,.c="hltMetCleanUsingJetID::HLT",.n="METJ",.t="JetID-cleaned MET",.f="METJ200",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  //_Steps["bMET170"]={.T=0,   .c="hltmetwithmu170",           .n="Full", .t="Full path",.f="bMET170",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // Set STEP style
  pair<Int_t, Int_t> theStyle = make_pair(0,0);
  for(_itSteps=_Steps.begin();_itSteps!=_Steps.end();_itSteps++) {
    _theStep = &(_itSteps->second);
    theStyle = getStyle(_theStep->n);
    _theStep->C = theStyle.first;
    _theStep->S = theStyle.second;
  }

  // PATHS //
  cout << "Define paths" << endl;
  Bool_t useHBHECleaning = (_HBHECleaning!="NoHBHE");
  vector<STEP> vStepEmpty;
  
  _Paths["PFMNoMu90"]={.nameP="PFMNoMu90",.namePath="HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["PFMNoMu90"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFMNoMu90"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFMNoMu90"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFMNoMu90"].steps.push_back(_Steps["MET65"]);
  if(useHBHECleaning) _Paths["PFMNoMu90"].steps.push_back(_Steps["METC55"]);
  _Paths["PFMNoMu90"].steps.push_back(_Steps["METJ55"]);
  _Paths["PFMNoMu90"].steps.push_back(_Steps["MHT65"]);
  _Paths["PFMNoMu90"].steps.push_back(_Steps["MuMHT90"]);
  _Paths["PFMNoMu90"].steps.push_back(_Steps["MuMET90"]);
  _Paths["PFMNoMu90"].steps.push_back(_Steps["bMuPFM90"]);
  _Paths["PFMNoMu90"].nSteps = _Paths["PFMNoMu90"].steps.size();

  _Paths["PFMNoMu120"]={.nameP="PFMNoMu120",.namePath="HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["PFMNoMu120"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFMNoMu120"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFMNoMu120"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFMNoMu120"].steps.push_back(_Steps["MET80"]);
  if(useHBHECleaning) _Paths["PFMNoMu120"].steps.push_back(_Steps["METC70"]);
  _Paths["PFMNoMu120"].steps.push_back(_Steps["METJ70"]);
  _Paths["PFMNoMu120"].steps.push_back(_Steps["MHT80"]);
  _Paths["PFMNoMu120"].steps.push_back(_Steps["MuMHT120"]);
  _Paths["PFMNoMu120"].steps.push_back(_Steps["MuMET120"]);
  _Paths["PFMNoMu120"].steps.push_back(_Steps["bMuPFM120"]);
  _Paths["PFMNoMu120"].nSteps = _Paths["PFMNoMu120"].steps.size();

  _Paths["Jet80PFMNoMu90"]={.nameP="Jet80PFMNoMu90",.namePath="HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["Jet80PFMNoMu90"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["MET65"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["Jet65"]);
  if(useHBHECleaning) _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["METC55"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["METJ55"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["MHT65"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["PFJet80"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["MuMHT90"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["MuMET90"]);
  _Paths["Jet80PFMNoMu90"].steps.push_back(_Steps["bJ80MuPFM90"]);
  _Paths["Jet80PFMNoMu90"].nSteps = _Paths["Jet80PFMNoMu90"].steps.size();

  _Paths["Jet80PFMNoMu120"]={.nameP="Jet80PFMNoMu120",.namePath="HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["Jet80PFMNoMu120"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["MET80"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["Jet65"]);
  if(useHBHECleaning) _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["METC70"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["METJ70"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["MHT80"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["PFJet80"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["MuMHT120"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["MuMET120"]);
  _Paths["Jet80PFMNoMu120"].steps.push_back(_Steps["bJ80MuPFM120"]);
  _Paths["Jet80PFMNoMu120"].nSteps = _Paths["Jet80PFMNoMu120"].steps.size();

  _Paths["PFM90"]={.nameP="PFM90",.namePath="HLT_PFMET90_PFMHT90_IDTight",.nSteps=6,.steps=vStepEmpty};
  _Paths["PFM90"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFM90"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFM90"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFM90"].steps.push_back(_Steps["MET70"]);
  _Paths["PFM90"].steps.push_back(_Steps["MHT70"]);
  _Paths["PFM90"].steps.push_back(_Steps["PFMHT90"]);
  _Paths["PFM90"].steps.push_back(_Steps["PFMET90"]);
  _Paths["PFM90"].steps.push_back(_Steps["bPFM90"]);
  _Paths["PFM90"].nSteps = _Paths["PFM90"].steps.size();

  _Paths["PFM120"]={.nameP="PFM120",.namePath="HLT_PFMET120_PFMHT120_IDTight",.nSteps=6,.steps=vStepEmpty};
  _Paths["PFM120"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFM120"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFM120"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFM120"].steps.push_back(_Steps["MET90"]);
  _Paths["PFM120"].steps.push_back(_Steps["MHT90"]);
  _Paths["PFM120"].steps.push_back(_Steps["PFMHT120"]);
  _Paths["PFM120"].steps.push_back(_Steps["PFMET120"]);
  _Paths["PFM120"].steps.push_back(_Steps["bPFM120"]);
  _Paths["PFM120"].nSteps = _Paths["PFM120"].steps.size();

  _Paths["PFMET170"]={.nameP="PFMET170",.namePath="HLT_PFMET170_JetIdCleaned",.nSteps=6,.steps=vStepEmpty};
  _Paths["PFMET170"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFMET170"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFMET170"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFMET170"].steps.push_back(_Steps["MET90"]);
  if(useHBHECleaning) _Paths["PFMET170"].steps.push_back(_Steps["METC80"]);
  _Paths["PFMET170"].steps.push_back(_Steps["METJ80"]);
  _Paths["PFMET170"].steps.push_back(_Steps["PFMET170"]);
  _Paths["PFMET170"].steps.push_back(_Steps["bMET170"]);
  _Paths["PFMET170"].nSteps = _Paths["PFMET170"].steps.size();

  _Paths["CaloMET200"]={.nameP="CaloMET200",.namePath="HLT_CaloMET200_JetIdCleaned",.nSteps=3,.steps=vStepEmpty};
  _Paths["CaloMET200"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["CaloMET200"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["CaloMET200"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["CaloMET200"].steps.push_back(_Steps["MET210"]);
  if(useHBHECleaning) _Paths["CaloMET200"].steps.push_back(_Steps["METC200"]);
  _Paths["CaloMET200"].steps.push_back(_Steps["METJ200"]);
  _Paths["CaloMET200"].nSteps = _Paths["CaloMET200"].steps.size();

  _Paths["PFMNoMu90_Up"]={.nameP="PFMNoMu90_Up",.namePath="HLT_CaloUp_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["PFMNoMu90_Up"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["MET80"]);
  if(useHBHECleaning) _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["METC70"]);
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["METJ70"]);
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["MHT80"]);
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["MuMHT90"]);
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["MuMET90"]);
  _Paths["PFMNoMu90_Up"].steps.push_back(_Steps["bMuPFM90"]);
  _Paths["PFMNoMu90_Up"].nSteps = _Paths["PFMNoMu90_Up"].steps.size();

  _Paths["Jet80PFMNoMu90_Up"]={.nameP="Jet80PFMNoMu90_Up",.namePath="HLT_CaloUp_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  _Paths["Jet80PFMNoMu90_Up"].steps.clear();
  //
  if(     _seed=="ETM50") _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["L1ETM50"]);
  else if(_seed=="ETM60") _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["L1ETM60"]);
  //
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["MET80"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["Jet65"]);
  if(useHBHECleaning) _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["METC70"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["METJ70"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["MHT80"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["PFJet80"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["MuMHT90"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["MuMET90"]);
  _Paths["Jet80PFMNoMu90_Up"].steps.push_back(_Steps["bMuPFM90"]);
  _Paths["Jet80PFMNoMu90_Up"].nSteps = _Paths["Jet80PFMNoMu90_Up"].steps.size();

  _Paths["OR90GeV"]={.nameP="OR90GeV",.namePath="HLT_PFMNoMu90_PFMET170",.nSteps=1,.steps=vStepEmpty};
  _Paths["OR90GeV"].steps.clear();
  _Paths["OR90GeV"].steps.push_back(_Steps["bOR90GeV"]);
  _Paths["OR90GeV"].nSteps = _Paths["OR90GeV"].steps.size();

  return 0;
}

Int_t MyTrigger::DefineJson()
{
 
  if(_json=="MC") {
    _applyJson = false;
  }
  else if(_json=="DCS") {
    _applyJson = true;
    if(_field=="38T")
      _jsonMap = readJSONFile(_dirJson+"/DCSOnly/json_DCSONLY.txt");
    else if(_field=="0T") 
      _jsonMap = readJSONFile(_dirJson+"/DCSOnly/json_DCSONLY_0T.txt");
  }
  else if(_json=="Prompt") {
    _applyJson = true;
    //
    if(_period=="25ns") {
      if(_field=="38T")
	_jsonMap = readJSONFile(_dirJson+"/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt");
      else if(_field=="0T") {
	_jsonMap = readJSONFile(_dirJson+"/Cert_246908-256869_13TeV_PromptReco_Collisions15_ZeroTesla_25ns_JSON.txt");
      }
      else {
	cout << "ERROR: Please choose field: 38T, 0T. Exit ==> []" << endl;
	return -1;
      }
    }
    //
    else if(_period=="50ns") {
      if(_field=="38T")
	_jsonMap = readJSONFile(_dirJson+"/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_v2.txt");
      else if(_field=="0T") {
	_jsonMap = readJSONFile(_dirJson+"/Cert_246908-252126_13TeV_PromptReco_Collisions15_ZeroTesla_50ns_JSON.txt");
      }
      else {
	cout << "ERROR: Please choose field: 38T, 0T. Exit ==> []" << endl;
	return -1;
      }
    }
    //
    else {
      cout << "ERROR: Please choose period: 25ns, 50ns. Exit ==> []" << endl;
      return -1;
    }
  }

  return 0;
}

Int_t MyTrigger::InitVar()
{

  _trig_obj_pt  = new vector<double> ();
  _trig_obj_eta = new vector<double> ();
  _trig_obj_phi = new vector<double> ();
  _trig_obj_col = new vector<string> ();

  InitEvent();

  return 0;
}

Int_t MyTrigger::InitEvent()
{

  _trig_obj_n = 0;
  _trig_obj_pt->clear();
  _trig_obj_eta->clear();
  _trig_obj_phi->clear();
  _trig_obj_col->clear();

  _event = _run = _lumi = 0;
  _xsec = _wgt = _kfact = _puwgt = 0;
  _puobs = _putrue = 0; 
  _nvtx = _nmuons = _nelectrons = _ntaus = _ntightmuons = _ntightelectrons = _nphotons = _njets = _nbjets = _nfatjets = 0;
  _hltmet90 = _hltmet120 = _hltmetwithmu90 = _hltmetwithmu120 = _hltmetwithmu170 = _hltmetwithmu300 = _hltjetmet90 = _hltjetmet120 = _hltphoton165 = _hltphoton175 = _hltdoublemu = _hltsinglemu = _hltdoubleel = _hltsingleel = 0;
  // _flagcsctight = _flaghbhenoise = _flaghcallaser = _flagecaltrig = _flageebadsc = _flagecallaser = _flagtrkfail = _flagtrkpog = _flaghnoiseloose = _flaghnoisetight = _flaghnoisehilvl = 0;
  _pfmet = _pfmetphi = _t1pfmet = _t1pfmetphi = _pfmupt = _pfmuphi = _mumet = _mumetphi = _phmet = _phmetphi = _t1mumet = _t1mumetphi = _t1phmet = _t1phmetphi = 0;
  _signaljetpt = _signaljeteta = _signaljetphi = _signaljetbtag = _signaljetCHfrac = _signaljetNHfrac = _signaljetEMfrac = _signaljetCEMfrac = _signaljetmetdphi = 0;
  _secondjetpt = _secondjeteta = _secondjetphi = _secondjetbtag = _secondjetCHfrac = _secondjetNHfrac = _secondjetEMfrac = _secondjetCEMfrac = _secondjetmetdphi = 0;
  _thirdjetpt  = _thirdjeteta  = _thirdjetphi  = _thirdjetbtag  = _thirdjetCHfrac  = _thirdjetNHfrac  = _thirdjetEMfrac  = _thirdjetCEMfrac  = _thirdjetmetdphi  = 0;
  _jetjetdphi = _jetmetdphimin = _incjetmetdphimin = 0;

  // leptons/photons
  _wzid = _l1id = _l2id = _i1id = _i2id = _i3id = _mu1pid = _mu2pid = _mu1id = _mu2id = _el1pid = _el2pid = _el1id = _el2id = 0; 
  _mu1pt =_mu1eta =_mu1phi =_mu2pt =_mu2eta =_mu2phi = 0;
  _el1pt =_el1eta =_el1phi =_el2pt =_el2eta =_el2phi =_phpt =_pheta =_phphi = 0;
  _loosephpt =_loosepheta =_loosephphi =_loosephsieie = _loosephrndiso = 0;
  
  return 0;
}

Int_t MyTrigger::GetInput()
{

  TString nameChain="tree/tree";
  if(_skim=="skim") nameChain="tree";
  _ch = new TChain(nameChain);

  vector<TString> f2015D = list_SingleMuon_2015D_V2();

  if(_era=="2015B")
    _ch->Add("/user/ndaci/Data/XMET/Run2015B/SingleMuon/V3/skim.root");
  else if(_era=="2015C") {
    if(_field=="0T") 
      _ch->Add("/user/ndaci/Data/XMET/Run2015C/SingleMuon/V2/skim_met30.root");
    else if(_field=="38T") 
      _ch->Add("/user/ndaci/Data/XMET/Run2015C/SingleMuon/38T_V5/skim_met30.root");
  }
  else if(_era=="2015D") {
    //_ch->Add("/user/ndaci/Data/XMET/Run2015D/SingleMuon/V1/tree_*.root");
    for(UInt_t iF=0 ; iF<f2015D.size() ; iF++) {
      _ch->Add(f2015D[iF]);
    }
  }
  else if(_era=="MC") {
    _ch->Add("/user/ndaci/Data/XMET/MonteCarloSpring15/V2/tree_*.root");
  }
  else {
    cout << "ERROR: Please specify input source in the output dir name: " 
	 << "2015B, 2015C, 2015C, 2015D. "
	 << "Please specify field: 38T, 0T. Exit ==> []"
	 << endl;
    return -1;
  }

  return 0;
}

Int_t MyTrigger::SetBranches()
{

  _ch->SetBranchStatus("*", 1);  

  // global informations
  _ch->SetBranchAddress("event", &_event); 
  _ch->SetBranchAddress("run", &_run); 
  _ch->SetBranchAddress("lumi", &_lumi); 

  // trigger objects
  _ch->SetBranchAddress("trig_obj_n", &_trig_obj_n); 
  _ch->SetBranchAddress("trig_obj_pt", &_trig_obj_pt);
  _ch->SetBranchAddress("trig_obj_eta", &_trig_obj_eta);
  _ch->SetBranchAddress("trig_obj_phi", &_trig_obj_phi);
  _ch->SetBranchAddress("trig_obj_col", &_trig_obj_col);

  // trigger bits
  _ch->SetBranchAddress("hltmet90", &_hltmet90); 
  _ch->SetBranchAddress("hltmet120", &_hltmet120); 
  _ch->SetBranchAddress("hltmetwithmu90", &_hltmetwithmu90); 
  _ch->SetBranchAddress("hltmetwithmu120", &_hltmetwithmu120); 
  _ch->SetBranchAddress("hltmetwithmu170", &_hltmetwithmu170); 
  _ch->SetBranchAddress("hltmetwithmu300", &_hltmetwithmu300); 
  _ch->SetBranchAddress("hltjetmet90", &_hltjetmet90); 
  _ch->SetBranchAddress("hltjetmet120", &_hltjetmet120); 
  _ch->SetBranchAddress("hltphoton165", &_hltphoton165); 
  _ch->SetBranchAddress("hltphoton175", &_hltphoton175); 
  _ch->SetBranchAddress("hltdoublemu", &_hltdoublemu); 
  _ch->SetBranchAddress("hltsinglemu", &_hltsinglemu); 
  _ch->SetBranchAddress("hltdoubleel", &_hltdoubleel); 
  _ch->SetBranchAddress("hltsingleel", &_hltsingleel); 

  // pile-up
  _ch->SetBranchAddress("puwgt", &_puwgt); 
  _ch->SetBranchAddress("puobs", &_puobs); 
  _ch->SetBranchAddress("putrue", &_putrue); 
  _ch->SetBranchAddress("nvtx", &_nvtx); 

  // mets
  _ch->SetBranchAddress("pfmet", &_pfmet); 
  _ch->SetBranchAddress("pfmetphi", &_pfmetphi); 
  _ch->SetBranchAddress("t1pfmet", &_t1pfmet); 
  _ch->SetBranchAddress("t1pfmetphi", &_t1pfmetphi);
  _ch->SetBranchAddress("pfmupt", &_pfmupt);
  _ch->SetBranchAddress("pfmuphi", &_pfmuphi);
  _ch->SetBranchAddress("mumet", &_mumet);
  _ch->SetBranchAddress("mumetphi", &_mumetphi);
  //_ch->SetBranchAddress("phmet", &_phmet); 
  //_ch->SetBranchAddress("phmetphi", &_phmetphi); 
  _ch->SetBranchAddress("t1mumet", &_t1mumet); 
  _ch->SetBranchAddress("t1mumetphi", &_t1mumetphi);
  //_ch->SetBranchAddress("t1phmet", &_t1phmet); 
  //_ch->SetBranchAddress("t1phmetphi", &_t1phmetphi);

  // jets
  _ch->SetBranchAddress("signaljetpt", &_signaljetpt);
  _ch->SetBranchAddress("signaljeteta", &_signaljeteta);
  _ch->SetBranchAddress("signaljetphi", &_signaljetphi);
  _ch->SetBranchAddress("signaljetbtag", &_signaljetbtag);
  _ch->SetBranchAddress("signaljetCHfrac", &_signaljetCHfrac);
  _ch->SetBranchAddress("signaljetNHfrac", &_signaljetNHfrac);
  _ch->SetBranchAddress("signaljetEMfrac", &_signaljetEMfrac);
  _ch->SetBranchAddress("signaljetCEMfrac", &_signaljetCEMfrac);
  _ch->SetBranchAddress("signaljetmetdphi", &_signaljetmetdphi);
  _ch->SetBranchAddress("secondjetpt", &_secondjetpt);
  _ch->SetBranchAddress("secondjeteta", &_secondjeteta);
  _ch->SetBranchAddress("secondjetphi", &_secondjetphi);
  _ch->SetBranchAddress("secondjetbtag", &_secondjetbtag);
  _ch->SetBranchAddress("secondjetCHfrac", &_secondjetCHfrac);
  _ch->SetBranchAddress("secondjetNHfrac", &_secondjetNHfrac);
  _ch->SetBranchAddress("secondjetEMfrac", &_secondjetEMfrac);
  _ch->SetBranchAddress("secondjetCEMfrac", &_secondjetCEMfrac); 
  _ch->SetBranchAddress("secondjetmetdphi", &_secondjetmetdphi); 
  _ch->SetBranchAddress("thirdjetpt", &_thirdjetpt); 
  _ch->SetBranchAddress("thirdjeteta", &_thirdjeteta);
  _ch->SetBranchAddress("thirdjetphi", &_thirdjetphi);
  _ch->SetBranchAddress("thirdjetbtag", &_thirdjetbtag);
  _ch->SetBranchAddress("thirdjetCHfrac", &_thirdjetCHfrac);
  _ch->SetBranchAddress("thirdjetNHfrac", &_thirdjetNHfrac);
  _ch->SetBranchAddress("thirdjetEMfrac", &_thirdjetEMfrac);
  _ch->SetBranchAddress("thirdjetCEMfrac", &_thirdjetCEMfrac);
  _ch->SetBranchAddress("thirdjetmetdphi", &_thirdjetmetdphi);

  // leptons/photons
  _ch->SetBranchAddress("wzid", &_wzid); 
  _ch->SetBranchAddress("l1id", &_l1id); 
  _ch->SetBranchAddress("l2id", &_l2id); 
  _ch->SetBranchAddress("i1id", &_i1id); 
  _ch->SetBranchAddress("i2id", &_i2id); 
  _ch->SetBranchAddress("i3id", &_i3id); 
  //
  _ch->SetBranchAddress("mu1pid", &_mu1pid); 
  _ch->SetBranchAddress("mu1pt", &_mu1pt); 
  _ch->SetBranchAddress("mu1eta", &_mu1eta);
  _ch->SetBranchAddress("mu1phi", &_mu1phi);
  _ch->SetBranchAddress("mu1id", &_mu1id); 
  _ch->SetBranchAddress("mu2pid", &_mu2pid);
  _ch->SetBranchAddress("mu2pt", &_mu2pt); 
  _ch->SetBranchAddress("mu2eta", &_mu2eta);
  _ch->SetBranchAddress("mu2phi", &_mu2phi);
  _ch->SetBranchAddress("mu2id", &_mu2id); 
  //
  _ch->SetBranchAddress("el1pid", &_el1pid); 
  _ch->SetBranchAddress("el1pt", &_el1pt); 
  _ch->SetBranchAddress("el1eta", &_el1eta);
  _ch->SetBranchAddress("el1phi", &_el1phi);
  _ch->SetBranchAddress("el1id", &_el1id); 
  _ch->SetBranchAddress("el2pid", &_el2pid);
  _ch->SetBranchAddress("el2pt", &_el2pt); 
  _ch->SetBranchAddress("el2eta", &_el2eta);
  _ch->SetBranchAddress("el2phi", &_el2phi);
  _ch->SetBranchAddress("el2id", &_el2id); 
  //
  _ch->SetBranchAddress("phpt", &_phpt); 
  _ch->SetBranchAddress("pheta", &_pheta);
  _ch->SetBranchAddress("phphi", &_phphi);
  _ch->SetBranchAddress("loosephpt", &_loosephpt); 
  _ch->SetBranchAddress("loosepheta", &_loosepheta);
  _ch->SetBranchAddress("loosephphi", &_loosephphi);
  _ch->SetBranchAddress("loosephsieie", &_loosephsieie);
  _ch->SetBranchAddress("loosephrndiso", &_loosephrndiso);

  // kinematics
  _ch->SetBranchAddress("jetjetdphi", &_jetjetdphi); 
  _ch->SetBranchAddress("jetmetdphimin", &_jetmetdphimin);
  _ch->SetBranchAddress("incjetmetdphimin", &_incjetmetdphimin);

  // object counters
  _ch->SetBranchAddress("nmuons", &_nmuons);
  _ch->SetBranchAddress("nelectrons", &_nelectrons);
  _ch->SetBranchAddress("ntightmuons", &_ntightmuons);
  _ch->SetBranchAddress("ntightelectrons", &_ntightelectrons);
  _ch->SetBranchAddress("ntaus", &_ntaus);
  _ch->SetBranchAddress("njets", &_njets);
  _ch->SetBranchAddress("nbjets", &_nbjets);
  _ch->SetBranchAddress("nfatjets", &_nfatjets);
  _ch->SetBranchAddress("nphotons", &_nphotons);

  // offline flags
  _ch->SetBranchAddress("flagcsctight", &_flagcsctight); 
  _ch->SetBranchAddress("flaghbhenoise", &_flaghbhenoise);
  _ch->SetBranchAddress("flaghcallaser", &_flaghcallaser);
  _ch->SetBranchAddress("flagecaltrig", &_flagecaltrig); 
  _ch->SetBranchAddress("flageebadsc", &_flageebadsc); 
  _ch->SetBranchAddress("flagecallaser", &_flagecallaser);
  _ch->SetBranchAddress("flagtrkfail", &_flagtrkfail); 
  _ch->SetBranchAddress("flagtrkpog", &_flagtrkpog); 
  _ch->SetBranchAddress("flaghnoiseloose", &_flaghnoiseloose);
  _ch->SetBranchAddress("flaghnoisetight", &_flaghnoisetight);
  _ch->SetBranchAddress("flaghnoisehilvl", &_flaghnoisehilvl);

  return 0;
}

Int_t MyTrigger::FillIneff()
{
  // this method is a bit risky because it calls steps[i] 
  // without knowing if there exists an element #i
  // => I should add a size checker

  (*_outIneff) << "INEFF: run " << _run 
	       << " lumi "      << _lumi
	       << " event "     << _event
	       << endl
	       << "hltmet120:"  << (_hltmet120 ? 1 : 0)
	       << " hltmetwithmu120:" << (_hltmetwithmu120 ? 1 : 0)
	       << " hltmetwithmu170:" << (_hltmetwithmu170 ? 1 : 0)
	       << endl << endl
	       << "1) Trigger objects" << endl
	       << "L1ETM="  << _Paths["PFMNoMu90"].steps[0].pt
	       << "   phi=" << _Paths["PFMNoMu90"].steps[0].phi 
	       << endl
	       << "MET="    << _Paths["PFMNoMu90"].steps[1].pt
	       << "   phi=" << _Paths["PFMNoMu90"].steps[1].phi 
	       << endl
	       << "HBHE-Cleaned MET="  << _Paths["PFMNoMu90"].steps[2].pt
	       << "   phi="            << _Paths["PFMNoMu90"].steps[2].phi 
	       << endl
	       << "JetID-Cleaned MET=" << _Paths["PFMNoMu90"].steps[3].pt
	       << "   phi="            << _Paths["PFMNoMu90"].steps[3].phi 
	       << endl
	       << "MHT="       << _Paths["PFMNoMu90"].steps[4].pt
	       << "   phi="    << _Paths["PFMNoMu90"].steps[4].phi 
	       << endl
	       << "PFMHTNoMu=" << _Paths["PFMNoMu90"].steps[5].pt
	       << "   phi="    << _Paths["PFMNoMu90"].steps[5].phi 
	       << endl
	       << "PFMETNoMu=" << _Paths["PFMNoMu90"].steps[6].pt
	       << "   phi="    << _Paths["PFMNoMu90"].steps[6].phi 
	       << endl
	       << "PFMHT="     << _Paths["PFMNoMu120"].steps[3].pt
	       << "   phi="    << _Paths["PFMNoMu120"].steps[3].phi 
	       << endl
	       << "PFMET="     << _Paths["PFMNoMu120"].steps[4].pt
	       << "   phi="    << _Paths["PFMNoMu120"].steps[4].phi 
	       << endl << endl
	       << "2) Offline MET objects"
	       << endl
	       << "mumet="   << _mumet
	       << "   phi="  << _mumetphi   
	       << endl
	       << "t1mumet=" << _t1mumet
	       << "   phi="  << _t1mumetphi 
	       << endl
	       << "pfmet="   << _pfmet
	       << "   phi="  << _pfmetphi   
	       << endl
	       << "t1pfmet=" << _t1pfmet
	       << "   phi="  << _t1pfmetphi 
	       << endl << endl
	       << "3) Offline objects"
	       << endl
	       << "Jet1(pt,eta,phi)="
	       << "(" << _signaljetpt 
	       << "," << _signaljeteta 
	       << "," << _signaljetphi 
	       << ")" << endl
	       << "Jet2(pt,eta,phi)="
	       << "(" << _secondjetpt 
	       << "," << _secondjeteta 
	       << "," << _secondjetphi 
	       << ")" << endl
	       <<"Jet3(pt,eta,phi)="
	       << "(" << _thirdjetpt 
	       << "," << _thirdjeteta 
	       << "," << _thirdjetphi 
	       << ")" << endl
	       << "Mu1(pt,eta,phi)="
	       << "(" << _mu1pt 
	       << "," << _mu1eta 
	       << "," << _mu1phi 
	       << ")" << endl
	       << "Mu2(pt,eta,phi)="
	       << "(" << _mu2pt 
	       << "," << _mu2eta 
	       << "," << _mu2phi 
	       << ")" << endl
	       << "El1(pt,eta,phi)="
	       << "(" << _el1pt 
	       << "," << _el1eta 
	       << "," << _el1phi 
	       << ")" << endl
	       << "El2(pt,eta,phi)="
	       << "(" << _el2pt 
	       << "," << _el2eta 
	       << "," << _el2phi 
	       << ")" << endl
	       << "Ph1(pt,eta,phi)="
	       << "(" << _loosephpt 
	       << "," << _loosepheta 
	       << "," << _loosephphi 
	       << ")" << endl
	       << endl;

  return 0;
}

pair<Int_t, Int_t> MyTrigger::getStyle(TString name)
{

  if(name=="L1")    return make_pair(kBlack , kOpenSquare);
  if(name=="MET")   return make_pair(kBlue+2, kOpenDiamond);
  if(name=="Jet")   return make_pair(kBlue+3, kOpenTriangleUp);
  if(name=="METC")  return make_pair(kBlue  , kOpenTriangleDown);
  if(name=="METJ")  return make_pair(kCyan+2, kFullTriangleUp);
  if(name=="MHT")   return make_pair(kGreen+2,kFullTriangleDown);
  if(name=="PFJET") return make_pair(kRed+3,  kFullDiamond);
  if(name=="PFMHT") return make_pair(kRed+2,  kOpenCircle);
  if(name=="PFMET") return make_pair(kRed,    kFullCircle);
  if(name=="Full")  return make_pair(kRed,    kFullCircle);

  return make_pair(kBlack,kFullCircle);
}
