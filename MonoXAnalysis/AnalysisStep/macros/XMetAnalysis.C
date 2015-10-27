#include "XMetAnalysis.h"

using namespace std;
Int_t verbose = 1000;

XMetAnalysis::XMetAnalysis(TString tag, TString subdir="AN")
{

  _tag = tag;
  _isAN15 = _tag.Contains("AN15");
  _isRun1 = _tag.Contains("Run1");

  _dirOut  = "/user/ndaci/Results/Monojet/"+subdir+"/";
  _outlog  = new ofstream(_dirOut+"/"+_tag+"/yields_"+_tag+".txt",ios::out);

  (*_outlog) << "- ::XMetAnalysis(" << tag << "," << subdir << endl;

  // Define the chains
  if(_isRun1)      DefineChainsRun1();
  else if(_isAN15) DefineChainsAN15();
  else             DefineChainsRun1();

  cout << "Will use trees for: ";
  if(_isRun1) cout << "8TeV analysis"  << endl;
  if(_isAN15) cout << "2015 analysis"  << endl;
  cout << "MC   trees: "  << _pathMC   << endl
       << "Data trees: "  << _pathData << endl
       << "Results in: "  << _dirOut+"/"+_tag+"/"
       << endl;
}

XMetAnalysis::~XMetAnalysis()
{
  (*_outlog) << "- ::~XMetAnalysis()" << endl;
  _outfile->Close();
  delete _outfile;
}

Int_t XMetAnalysis::Analysis(Bool_t bProdHistos)
{

  (*_outlog) << "- ::Analysis()" << endl;
  if(_isRun1)      {
    AnalysisRun1(bProdHistos);
  (*_outlog) << "- ::AnalysisRun1()" << endl;
  }
  else if(_isAN15) {
    AnalysisAN15(bProdHistos);
    (*_outlog) << "- ::AnalysisAN15()" << endl;
  }
  else {
    AnalysisRun1(bProdHistos);
    (*_outlog) << "- ::AnalysisRun1()" << endl;
  }

  return 0;
}

Int_t XMetAnalysis::AnalysisAN15(Bool_t bProdHistos)
{
  (*_outlog) << "- ::AnalysisAN15()" << endl;

  TString mode;
  if(bProdHistos) mode="recreate";
  else            mode="read";
  _outfile = new TFile(_dirOut+"/"+_tag+"/plots_"+_tag+".root",mode);
  _outfile->cd();

  // Processes to use
  vector<TString> locProcesses;
  //
  locProcesses.push_back("data_met");
  /*
  locProcesses.push_back("data_2m");
  locProcesses.push_back("data_1m");
  locProcesses.push_back("data_1ph");
  locProcesses.push_back("data_2e");
  */
  //
  locProcesses.push_back("znn"); 
  locProcesses.push_back("zll"); 
  locProcesses.push_back("wln"); 
  locProcesses.push_back("ttbar"); 
  locProcesses.push_back("top"); 
  locProcesses.push_back("vv"); 
  locProcesses.push_back("qcd"); 

  // Process to put in the stack plot
  vector<TString> stackProcesses;
  //
  stackProcesses.push_back("data_met");
  //
  stackProcesses.push_back("qcd"); 
  stackProcesses.push_back("znn"); 
  stackProcesses.push_back("zll"); 
  stackProcesses.push_back("wln"); 
  stackProcesses.push_back("ttbar"); 
  stackProcesses.push_back("top"); 
  stackProcesses.push_back("vv"); 

  // Selections and variables
  const UInt_t nS=4;
  TString select[nS] = {"1jet","2jet","3jet","4jet"};
  TString selectTitle[nS] = {"1 Jet", "2 Jets", "3 Jets", "4 Jets And More"};
  for(UInt_t iS=0 ; iS<nS ; iS++) _Title[ select[iS] ] = selectTitle[iS];

  /*
  const UInt_t nCut=5;
  TString scanCut[nCut] = {"JetMet0p4","JetMet0p45","JetMet0p5","JetMet0p55","JetMet0p6"};
  Bool_t  scanReset[nCut] = {true,false,false,false,false};
  */

  const UInt_t nCut=1;
  TString scanCut[  nCut] = {"NoJmCut"};
  //TString scanCut[  nCut] = {"JetMet0p4"}; // fixme
  Bool_t  scanReset[nCut] = {true};

  // Variables
  const UInt_t nV=3;
  TString var[nV]     ={"t1mumet", "jetmetdphimin", "incjetmetdphimin"};
  TString nameAxis[nV]={"Type1 PFMETNoMu [GeV]", "Min #Delta#phi(M,J_{i}^{C})", "Min #Delta#phi(M,J_{i})"};
  for(UInt_t iV=0; iV<nV; iV++) {
    _Axis[var[iV]] = nameAxis[iV];
  }

  // Regular Binning
  UInt_t   nBins[nV]={ 100,   64,   64};
  Float_t xFirst[nV]={ 200,    0,    0};
  Float_t xLast[ nV]={1000,  3.2,  3.2};
  Bool_t regular[nV]={true,true,true};

  // Tuned binning
  Float_t* binsTuned[nV];  
  /// MET
  regular[0] = false;
  nBins[0] = 8;
  Float_t bins_met[8] = {200, 250, 300, 350, 400, 500, 600, 1000};
  binsTuned[0]=bins_met;

  vector<vector<Float_t>> v_bins;
  vector<Float_t> theBins;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    theBins.clear();
    for(UInt_t iB=0 ; iB<nBins[iV] ; iB++) {
      if(regular[iV]) {
	if(nBins[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBins[iV]);
	else             binval = 0;
	theBins.push_back(binval);
      }
      else {
	theBins.push_back(binsTuned[iV][iB]);
      }
    }
    v_bins.push_back(theBins);
  }

  // Produce individual histograms
  if(bProdHistos) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      if(verbose>1) cout << "- selection : " << select[iS] << endl;
      plot(select[iS], nCut, scanCut, scanReset, nV, var, v_bins, locProcesses);
    }
  }
  else {
    GetHistos(nS, select, nCut, scanCut, nV, var, locProcesses);
  }

  DrawStackPlots(nS, select, nCut, scanCut, nV, var, stackProcesses);

  return 0;
}

Int_t XMetAnalysis::AnalysisRun1(Bool_t bProdHistos)
{

  (*_outlog) << "- ::AnalysisRun1()" << endl;

  TString mode;
  if(bProdHistos) mode="recreate";
  else            mode="read";
  _outfile = new TFile(_dirOut+"/"+_tag+"/plots_"+_tag+".root",mode);
  _outfile->cd();

  // Processes to use
  vector<TString> locProcesses;

  // Data
  locProcesses.push_back("data");

  // Signal
  const UInt_t nSpin=2;
  const UInt_t nMass=9;
  TString nameSpin[nSpin] = {"V","AV"};
  TString nameMass[nMass] = {"0p1","1","10","100","200","300","400","700","1000"};
  TString nameProcess;
  for(UInt_t iS=0 ; iS<nSpin ; iS++) {
    for(UInt_t iM=0 ; iM<nMass ; iM++) {
      nameProcess = "dm_"+nameSpin[iS]+"_"+nameMass[iM];
      locProcesses.push_back(nameProcess); 
    }
  }

  // MC backgrounds
  locProcesses.push_back("znn"); 
  locProcesses.push_back("zll"); 
  locProcesses.push_back("wln"); 
  locProcesses.push_back("ttbar"); 
  locProcesses.push_back("top"); 
  locProcesses.push_back("vv"); 
  locProcesses.push_back("qcd"); 

  // Data-driven
  /// without acc, eff, BR corrections
  locProcesses.push_back("zn_cr_DA"); 
  locProcesses.push_back("zn_cr_MC"); 
  locProcesses.push_back("zn_sr_DD"); 
  locProcesses.push_back("wj_cr_DA"); 
  locProcesses.push_back("wj_cr_MC"); 
  locProcesses.push_back("wj_sr_DD"); 

  /// with acc, eff, BR corrections
  locProcesses.push_back("zn_cr_corr_DA"); 
  locProcesses.push_back("zn_cr_corr_MC"); 
  locProcesses.push_back("zn_sr_corr_DD"); 
  locProcesses.push_back("wj_cr_corr_DA"); 
  locProcesses.push_back("wj_cr_corr_MC"); 
  locProcesses.push_back("wj_sr_corr_DD"); 

  // Selections and variables
  const UInt_t nS=3;
  TString select[nS] = {"1jet","2jet","3jet"};
  /*
  const UInt_t nS=1;
  TString select[nS] = {"monojet"};
  */

  const UInt_t nCut=5;
  TString scanCut[  nCut] = {"NoJmCut","JetMet0p2","JetMet0p4","JetMet0p6","JetMet0p8"};
  Bool_t  scanReset[nCut] = {true,true,true,true,true};
  /*
  const UInt_t nCut=1;
  TString scanCut[  nCut] = {"NoJmCut"};
  Bool_t  scanReset[nCut] = {true};
  */

  // variables
  const UInt_t nV=1;
  TString var[nV]    = {"t1mumet"};
  TString nameAxis[nV]={"Type1 PFMETNoMu [GeV]"};
  for(UInt_t iV=0; iV<nV; iV++) {
    _Axis[var[iV]] = nameAxis[iV];
  }

  // Binning
  //
  /// regular binning
  UInt_t  nBins[ nV]={50};
  Float_t xFirst[nV]={0};
  Float_t xLast[ nV]={1000};
  Bool_t regular[nV]={1};
  //
  /// tuned binning
  Float_t* binsTuned[nV];  
  //
  vector<vector<Float_t>> v_bins;
  vector<Float_t> theBins;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    theBins.clear();
    for(UInt_t iB=0 ; iB<nBins[iV] ; iB++) {
      if(regular[iV]) {
	if(nBins[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBins[iV]);
	else             binval = 0;
	theBins.push_back(binval);
      }
      else {
	theBins.push_back(binsTuned[iV][iB]);
      }
    }
    v_bins.push_back(theBins);
  }

  if(bProdHistos) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      if(verbose>1) cout << "- selection : " << select[iS] << endl;
      plot(select[iS], nCut, scanCut, scanReset, nV, var, v_bins, locProcesses);
    }
  }
  else {
    GetHistos(nS, select, nCut, scanCut, nV, var, locProcesses);
  }

  return 0;
}

Int_t XMetAnalysis::CheckForwardJets(Bool_t bProdHistos)
{

  (*_outlog) << "- ::AnalysisAN15()" << endl;

  TString mode;
  if(bProdHistos) mode="recreate";
  else            mode="read";
  _outfile = new TFile(_dirOut+"/"+_tag+"/plots_"+_tag+".root",mode);
  _outfile->cd();

  // Processes to use
  vector<TString> locProcesses;
  //
  locProcesses.push_back("data_met");
  /*
  locProcesses.push_back("data_2m");
  locProcesses.push_back("data_1m");
  locProcesses.push_back("data_1ph");
  locProcesses.push_back("data_2e");
  */
  //
  locProcesses.push_back("znn"); 
  locProcesses.push_back("zll"); 
  locProcesses.push_back("wln"); 
  locProcesses.push_back("ttbar"); 
  locProcesses.push_back("top"); 
  locProcesses.push_back("vv"); 
  locProcesses.push_back("qcd"); 

  // Process to put in the stack plot
  vector<TString> stackProcesses;
  //
  stackProcesses.push_back("data_met");
  //
  stackProcesses.push_back("qcd"); 
  stackProcesses.push_back("znn"); 
  stackProcesses.push_back("zll"); 
  stackProcesses.push_back("wln"); 
  stackProcesses.push_back("ttbar"); 
  stackProcesses.push_back("top"); 
  stackProcesses.push_back("vv"); 

  // Selections and variables
  const UInt_t nS=4;
  TString select[nS] = {"1jet","2jet","3jet","4jet"};
  TString selectTitle[nS] = {"1 Jet", "2 Jets", "3 Jets", "4 Jets And More"};
  for(UInt_t iS=0 ; iS<nS ; iS++) _Title[ select[iS] ] = selectTitle[iS];

  const UInt_t nCut=3;
  TString scanCut[  nCut] = {"NoJmCut","JetMetBelow0p4","JetMet0p4"}; 
  Bool_t  scanReset[nCut] = {true,false,true};

  // Variables
  const UInt_t nV=9;
  TString var[nV]={"njets", 
		   //"nsoftjets", "njets80", 
		   //"nsoftfwdjets", "nfwdjets", "nfwdjets80",
		   "leadfwdjet_pt", "leadfwdjet_eta", "leadfwdjet_phi", 
		   "leadfwdjet_CHfrac" , "leadfwdjet_NHfrac", "leadfwdjet_EMfrac",
		   "leadfwdjet_CEMfrac", "leadfwdjet_Mufrac"};

  TString nameAxis[nV]={"# Central Jets (p_{T}>30 GeV ; |#eta|<2.5)", 
			//"# Soft Central Jets (p_{T}<30 GeV ; |#eta|<2.5)", 
			//"# Hard Central Jets (p_{T}>80 GeV ; |#eta|<2.5)", 
			//"# Soft Forward Jets (p_{T}<30 GeV ; |#eta|<2.5)", 
			//"# Forward Jets (p_{T}>30 GeV ; |#eta|<2.5)", 
			//"# Hard Forward Jets (p_{T}>80 GeV ; |#eta|<2.5)", 
			"Leading Forward Jet p_{T} [GeV]",
			"Leading Forward Jet #eta [rad]",
			"Leading Forward Jet #phi [rad]",
			"Leading Forward Jet Charged Hadronic Energy fraction",
			"Leading Forward Jet Neutral Hadronic Energy fraction",
			"Leading Forward Jet Neutral EM Energy fraction",
			"Leading Forward Jet Charged EM Energy fraction",
			"Leading Forward Jet Muon Energy fraction"};

  for(UInt_t iV=0; iV<nV; iV++) {
    _Axis[var[iV]] = nameAxis[iV];
  }

  // Binning
  //
  /// regular binning
  UInt_t  nBins[ nV]={20, /*20, 20, 20, 20, 20,*/  100, 100,  64, 100,100,100,100,100};
  Float_t xFirst[nV]={ 0, /* 0,  0,  0,  0,  0,*/    0,  -5,   0,   0,  0,  0,  0,  0};
  Float_t xLast[ nV]={20, /*20, 20, 20, 20, 20,*/ 1000,   5, 3.2,   1,  1,  1,  1,  1};
  Bool_t regular[nV]={ 1, /*1,1,1,1,1,*/             1,   1,   1,   1,  1,  1,  1,  1};
  //
  /// tuned binning
  Float_t* binsTuned[nV];  
  //
  vector<vector<Float_t>> v_bins;
  vector<Float_t> theBins;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    theBins.clear();
    for(UInt_t iB=0 ; iB<nBins[iV] ; iB++) {
      if(regular[iV]) {
	if(nBins[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBins[iV]);
	else             binval = 0;
	theBins.push_back(binval);
      }
      else {
	theBins.push_back(binsTuned[iV][iB]);
      }
    }
    v_bins.push_back(theBins);
  }

  // Produce individual histograms
  cout << "- Start producing individual hisograms" << endl;
  if(bProdHistos) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      if(verbose>1) cout << "-- selection : " << select[iS] << endl;
      plot(select[iS], nCut, scanCut, scanReset, nV, var, v_bins, locProcesses);
    }
  }
  else {
    GetHistos(nS, select, nCut, scanCut, nV, var, locProcesses);
  }

  DrawStackPlots(nS, select, nCut, scanCut, nV, var, stackProcesses);

  return 0;
}

Int_t XMetAnalysis::StudyQCDKiller(TString signal="znn")
{

  (*_outlog) << "- ::StudyQCDKiller(" << signal << ")" << endl;

  // Processes to use
  vector<TString> locProcesses;
  locProcesses.push_back(signal); 
  locProcesses.push_back("qcd"); 

  // Selections and variables
  const UInt_t nS=6;
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet","4jet"};

  const UInt_t nCut=1;
  TString scanCut[  nCut] = {"NoCut"};
  Bool_t  scanReset[nCut] = {true};

  const UInt_t nV=5;
  TString var[nV]    = {"jetmetdphimin"     , "incjetmetdphimin",
			//"signaljetmetdphi"  , "secondjetmetdphi", 
			//"thirdjetmetdphi"   , 
			"jetjetdphi"      , 
			//"cosjetjetdphiover2", "abscosjetjetdphiover2",
			//"dphiJ1J3"          , "dphiJ2J3",
			"apcjetmetmax"      , "apcjetmetmin"};
			//"alphat"}; // removed from the trees because memory issues

  // Binning
  //
  /// regular binning
  UInt_t  nBins[nV]  = {10000, 10000, 10000, 10000, 10000}; //, 10000, 10000, 10000, 10000, 10000, 10000, 10000};//, 8000};
  Float_t xFirst[nV] = {    0,     0,     0,     0,     0}; //,     0,     0,     0,     0,     0,     0,     0};//,    0};
  Float_t xLast[nV]  = {  3.2,   3.2,   3.2,   3.2,   3.2}; //,   3.2,   3.2,   3.2,   3.2,   3.2,   3.2,   3.2};//,    2};
  Bool_t regular[nV]={1,1,1,1,1};
  //
  /// tuned binning
  Float_t* binsTuned[nV];  
  //
  vector<vector<Float_t>> v_bins;
  vector<Float_t> theBins;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    theBins.clear();
    for(UInt_t iB=0 ; iB<nBins[iV] ; iB++) {
      if(regular[iV]) {
	if(nBins[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBins[iV]);
	else             binval = 0;
	theBins.push_back(binval);
      }
      else {
	theBins.push_back(binsTuned[iV][iB]);
      }
    }
    v_bins.push_back(theBins);
  }

  // Produce 1 plot per {selection ; variable}
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    if(verbose>1) cout << "- selection : " << select[iS] << endl;
    plot(select[iS], nCut, scanCut, scanReset, nV, var, v_bins, locProcesses);
  }

  if(_outfile) _outfile->Close();
  else {
    cout << "ERROR: outfile undefined! Exit ==> []" << endl;
    return -2;
  }

  return 0;
}

Int_t XMetAnalysis::plot(TString select, 
			 const UInt_t nCut, TString *scanCut, Bool_t *scanReset,
			 const UInt_t nV,    TString *var, 
			 vector<vector<Float_t>> v_bins,
			 vector<TString> locProcesses)
{

  // Write the arguments in the log file /////////////////////////
  (*_outlog) << "- ::plot(" << select << " , " << nCut << " , {" ;
  for(UInt_t iCut=0; iCut<nCut; iCut++) {
    (*_outlog) << scanCut[iCut] << "(" << scanReset[iCut] << ")";
    if(iCut<nCut-1) (*_outlog) << " , ";
  }
  //
  (*_outlog) << "} , " << nV << " , {";
  /*
  for(UInt_t iV=0; iV<nV; iV++) {
    (*_outlog) << var[iV] << "(" << nBins[iV] << "," << xFirst[iV] << "," << xLast[iV] << ")";
    if(iV<nV-1) (*_outlog) << " , ";
  }
  */
  //
  (*_outlog) << "} , {";
  const UInt_t nP = locProcesses.size();
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    (*_outlog) << locProcesses[iP];
    if(iP<nP-1) (*_outlog) << ",";
  }
  //
  (*_outlog) << "} )" << endl << endl;
  //////////////////////////////////////////////////////

  (*_outlog) << "Selection: " << select << endl;

  // Define selections //
  TString selectScan;
  TCut weight="";
  TString region = "signal";
  TCut cut = defineCut(select, region);

  // Declare histograms //
  //M_PROCESS_CUT_VAR_H mapHistos; // made it a member
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax;

  // Loop over chains and generate histograms //
  //
  TString nameDir, locVar, variable;
  Int_t color, style, size;
  TH1F *hTemp, *hSRDD, *hTemp_cr_d, *hTemp_cr_mc,
    *hTemp_num, *hTemp_den;
  pair<Double_t, Double_t> intErr;
  //
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    nameDir = locProcesses[iP];

    // SR data-driven histograms don't have chains
    if(nameDir.Contains("DD")) continue;

    // check if current requested process is available
    if(_mapProcess.find(nameDir)==_mapProcess.end()) { 
      if(verbose>1) {
	cout << "-- ERROR: requested process '" 
	     << locProcesses[iP] << "' unavailable." << endl;
      }
      continue;
    }
    //
    if(verbose>1) cout << "-- process : " << nameDir << endl;

    // Define reweighting
    if( nameDir.Contains("data") ||
	nameDir.Contains("DA") ) {
      weight = "1";
    }
    else { 
      if(_isRun1) weight = "puwgt*wgt"; 
      //else if(_isAN15) weight = "xsec*(wgt/wgtsum)*puwgt";  // FIXME
      //else if(_isAN15) weight = "xsec*puwgt"; // my xsec are buggy => check theXSec in crab
      else if(_isAN15) weight = "editxsec*(wgt/wgtsum)*puwgt"; // using editxsec produced at skim 
    }

    // Choose the phase space region
    region = "signal";
    //
    if(     nameDir.Contains("zn_cr")) {
      if(nameDir.Contains("corr")) {
	region  = "zctrl_corr";
	weight *= "((5.942 * 1.023) / (0.79*(1.0-TMath::Exp(-0.00910276*(t1mumet-36.1669)))))";
      }
      else { region = "zctrl"; }
    }
    //
    else if(nameDir.Contains("wj_cr"))  {
      if(nameDir.Contains("corr")) {
	region  = "wctrl_corr";
	weight *= "(38.7823/TMath::Power(-85.7023 + t1mumet, 0.667232))";
      }
      else { region = "wctrl"; }
    }
    //
    else {region = "signal";}

    // Apply k-factors
    if(_useLO) {
      if(nameDir.Contains("znn") || nameDir.Contains("zll") ) {
	weight *= "1.35";
      }
      if(nameDir.Contains("wln") ) {
	weight *= "1.2";
      }
    }
    ////////////////////////////////////

    // First skim based on select only:
    /////cut = defineCut(select, region);
    /////_mapProcess[nameDir].Skim(select, cut, "reset");

    // Loop over cuts to be scanned
    for(UInt_t iCut=0; iCut<nCut; iCut++) {

      selectScan = select+"_"+scanCut[iCut];
      cut = defineCut(selectScan, region);

      if(verbose>1) cout << "-- cut: " << cut << " selectScan: " 
			 << selectScan << " region: " << region 
			 << endl;
      
      // Skim the chain
      /////_mapProcess[nameDir].Skim(select,     cut, "reuse");
      _mapProcess[nameDir].Skim(selectScan, cut, "reset");

      /*
      if(iCut>0 && scanReset[iCut]) {
	_mapProcess[nameDir].Skim(select,     cut, "reuse");
      }
      _mapProcess[nameDir].Skim(selectScan, cut, "");
      */

      // Loop over requested variables
      for(UInt_t iV=0 ; iV<nV ; iV++) {

	if(verbose>1) cout << "--- variable: " << var[iV] << endl;

	// Re-Build binning 
	const UInt_t nBins = v_bins[iV].size();
	Float_t theBins[nBins];
	for(UInt_t iB=0 ; iB<nBins ; iB++) {
	  theBins[iB] = v_bins[iV][iB];
	}

	// Define histogram and set style
	_mapHistos[select][nameDir][scanCut[iCut]][var[iV]] = 
	  new TH1F("h_"+var[iV]+"_"+nameDir+"_"+selectScan, 
		   var[iV]+" "+nameDir+" "+selectScan,
		   nBins-1, theBins);

	hTemp = _mapHistos[select][nameDir][scanCut[iCut]][var[iV]];
	color = _mapProcess[nameDir].GetColor();
	style = _mapProcess[nameDir].GetStyle();
	size  = _mapProcess[nameDir].GetSize();
	setStyle( hTemp , color , style , size , kTRUE , kTRUE );
	hTemp->SetXTitle(_Axis[var[iV]]);

	// Draw the variable
	locVar = var[iV];
	if(var[iV].Contains("phi"))  locVar = "abs("+var[iV]+")";
	if(     var[iV]=="dphiJ1J3") locVar = "abs(signaljetphi-thirdjetphi)";
	else if(var[iV]=="dphiJ2J3") locVar = "abs(secondjetphi-thirdjetphi)";
	else if(var[iV]=="leadjetmetdphi")        locVar = "abs(signaljetphi-t1mumetphi)";
	else if(var[iV]=="cosjetjetdphiover2")    locVar = "cos(jetjetdphi/2.)";
	else if(var[iV]=="abscosjetjetdphiover2") locVar = "abs(cos(jetjetdphi/2.))";

	// Jet Arrays
	if(var[iV].Contains("leadfwd")) {
	  cut *= "(jet_pt==Max$(jet_pt * (abs(jet_eta)>2.5)) && abs(jet_eta)>2.5)";
	  locVar = var[iV](7,var[iV].Sizeof());
	}

	_mapProcess[nameDir].Draw(hTemp, locVar, cut, weight);
	if(verbose>1) cout << "--- drew variable: " << locVar << endl;

	// Normalize
	Double_t theScale=1.0;
	if( !nameDir.Contains("data") && !nameDir.Contains("DA") ) {
	  theScale = _lumi*_rescale;
	  if(nameDir.Contains("qcd")) theScale *= _qcdScale;
	}
	hTemp->Scale(theScale);

	// Checks
	intErr = Integrate(hTemp);
	if(verbose>1) {
	  cout << "--- check:" 
	       << " nameDir="  << nameDir
	       << " GetName="  << hTemp->GetName()
	       << " Entries="  << hTemp->GetEntries()
	       << " Integral=" << intErr.first
	       << " +/- "      << intErr.second
	       << " Weight="   << weight
	       << endl;
	}

	// Determine extrema
	locMin = hTemp->GetMinimum();
	locMax = hTemp->GetMaximum();
	if(locMin<minPlot) minPlot = locMin;
	if(locMax>maxPlot) maxPlot = locMax;

      } // end loop:variables
    } // end loop:cuts
  } // end loop:processes 
  if(verbose>1) cout << "-- end loop over processes" << endl;
  

  // DATA-DRIVEN //

  // Prepare DD parameters
  if(verbose>1) cout << "-- Ready to produce sr_DD histograms" << endl;
  const UInt_t nDD=2;   // 2 processes to be evaluated: Z(nn), W(ln)
  const UInt_t nCorr=2; // 2 DD methods
  pair<Double_t, Double_t> theNum, theDen, theSF;
  //
  // DD scale-factor-like
  TString myNum[nDD]    = {"znn","wln"};
  TString myDen[nDD]    = {"zn_cr_MC","wj_cr_MC"};
  // DD corrected-like
  TString myDD[nDD]     = {"zn","wj"};
  TString myCorr[nCorr] = {"_","_corr_"};
  //

  // Loop over cuts to be scanned
  for(UInt_t iCut=0; iCut<nCut; iCut++) {
    selectScan = select+"_"+scanCut[iCut];

    for(UInt_t iDD=0 ; iDD<nDD ; iDD++) {
    
      // Check that the current process (znn or wj) has been requested initially
      bool srddIsHere=false;
      for(UInt_t iP=0 ; iP<nP ; iP++) {
	if(locProcesses[iP].Contains(myDD[iDD]+"_cr")) {
	  srddIsHere=true;
	  break; 
	}
      }
      if(!srddIsHere) {
	if(verbose>1) cout << "--- unrequested process: " << myDD[iDD] << endl;
	continue;
      }

      for(UInt_t iCorr=0 ; iCorr<nCorr ; iCorr++) {
	for(UInt_t iV=0 ; iV<nV ; iV++) {
	  //
	  nameDir = "h_"+var[iV]+"_"+myDD[iDD]+"_sr"+myCorr[iCorr]+"DD_"+selectScan;
	  if(verbose>1) cout << "----- produce: " << nameDir << endl;
	  //
	  // Method1: corrected-like
	  if(myCorr[iCorr].Contains("corr")) {
	    hTemp_cr_d  = _mapHistos[select][myDD[iDD]+"_cr"+myCorr[iCorr]+"DA"][scanCut[iCut]][var[iV]];
	    hTemp_cr_mc = _mapHistos[select][myDD[iDD]+"_cr"+myCorr[iCorr]+"MC"][scanCut[iCut]][var[iV]];
	    if(!hTemp_cr_d || !hTemp_cr_mc) continue;
	    //
	    hSRDD = (TH1F*)hTemp_cr_d->Clone(nameDir);
	    if(hTemp_cr_mc) hSRDD->Add( hTemp_cr_mc , -1.0 );
	    _mapHistos[select][myDD[iDD]+"_sr"+myCorr[iCorr]+"DD"][scanCut[iCut]][var[iV]] = hSRDD;
	  }
	  //
	  // Method2: SF-like
	  else {
	    hTemp_cr_d  = _mapHistos[select][myDD[iDD]+"_cr"+myCorr[iCorr]+"DA"][scanCut[iCut]][var[iV]];
	    hTemp_num = _mapHistos[select][myNum[iDD]][scanCut[iCut]][var[iV]];
	    hTemp_den = _mapHistos[select][myDen[iDD]][scanCut[iCut]][var[iV]];
	    if(!hTemp_cr_d || !hTemp_num || !hTemp_den) continue;
	    //
	    hSRDD = (TH1F*)hTemp_cr_d->Clone(nameDir);
	    theSF.first = 1.0;
	    theSF.second= 1.0;
	    if(hTemp_num && hTemp_den) {
	      theNum = Integrate(hTemp_num);
	      theDen = Integrate(hTemp_den);
	      theSF.first  = theDen.first!=0 ? theNum.first / theDen.first : 1.0;
	    }
	    else {
	      theSF.first = 1.0;
	    }
	    hSRDD->Scale( theSF.first );
	    _mapHistos[select][myDD[iDD]+"_sr"+myCorr[iCorr]+"DD"][scanCut[iCut]][var[iV]] = hSRDD;
	    cout << "----- SF(" << nameDir 
		 << ")=" << theNum.first 
		 << "/"  << theDen.first
		 << "="  << theSF.first
		 << endl;
	  }

	  intErr = Integrate(hSRDD);
	  if(verbose>1) {
	    cout << "----- check:" 
		 << " GetName="  << hSRDD->GetName()
		 << " Entries="  << hSRDD->GetEntries()
		 << " Integral=" << intErr.first
		 << " +/- "      << intErr.second;
	      //<< " Weight="   << weight
	    if(!myCorr[iCorr].Contains("corr")) 
	      cout << " SF="  << theSF.first;
	    cout << endl;
	  }

	} // end loop:var
      } // end loop:corr options
    } // end loop:DD processes
  } //end loop:cuts to be scanned
  //////////////////////////

  // Yields outlog //

  // Loop over variables
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    // Print out variable and processes names
    (*_outlog) << "Var: " << var[iV] << endl;

    // Write 1 column title per scanned cut
    (*_outlog) << setw(14) << "Process" ;
    for(UInt_t iCut=0; iCut<nCut; iCut++) {
      (*_outlog) << setw(23) << scanCut[iCut] ;
    }
    (*_outlog) << endl;

    // Loop over processes
    for(UInt_t iP=0 ; iP<nP ; iP++) {

      // Write out process name
      (*_outlog) << setw(14) << locProcesses[iP] ;

      // Get histogram
      for(UInt_t iCut=0; iCut<nCut; iCut++) {
	hTemp = _mapHistos[select][locProcesses[iP]][scanCut[iCut]][var[iV]];
	if(!hTemp) {
	  if(verbose>1) {
	    cout << "ERROR : yield loop did not find histo:"
		 << locProcesses[iP]+" "+var[iV]
		 << endl;
	  }
	  (*_outlog) << setw(23) << "ERROR" ;
	  continue;
	}
	else {
	  intErr = Integrate(hTemp);
	  if(verbose>2) cout << " ==> " << intErr.first 
			     << " +/- " << intErr.second 
			     << endl;
	  (*_outlog) << setw(10) << intErr.first 
		     << setw(5)  << "+/-"
		     << setw(8) << intErr.second;
	}
      } // end loop:cuts

      if(verbose>2) cout << "end loop:cuts" << endl;
      (*_outlog) << endl;

    }   // end loop:processes

    if(verbose>2) cout << "end loop:processes" << endl;
    (*_outlog) << endl;

  } // end loop over var

  if(verbose>2) cout << "end loop:var" << endl;  
  (*_outlog) << endl;

  //////////////////////////

  // Save histograms //
  if(_outfile) _outfile->cd();
  else {
    cout << "ERROR: outfile undefined! Exit ==> []" << endl;
    return -2;
  }

  if(verbose>2) cout << "-- _outfile->cd()... done" << endl;
  TString nameHisto;
  for( _itSelProcCutVarH=_mapHistos.begin() ; 
       _itSelProcCutVarH!=_mapHistos.end() ; 
       _itSelProcCutVarH++) {

    for( _itProcCutVarH=_itSelProcCutVarH->second.begin() ; 
	 _itProcCutVarH!=_itSelProcCutVarH->second.end() ; 
	 _itProcCutVarH++) {

      for( _itCutVarH=_itProcCutVarH->second.begin() ; 
	   _itCutVarH!=_itProcCutVarH->second.end() ; 
	   _itCutVarH++) {

	for( _itVarH=_itCutVarH->second.begin() ; 
	     _itVarH!=_itCutVarH->second.end() ; 
	     _itVarH++) {

	hTemp = _itVarH->second;

	if(hTemp) {
	  nameHisto = hTemp->GetName();
	  hTemp->Write();
	  if(verbose>2) cout << "---- histo written: " << nameHisto << endl;
	  //delete hTemp; // ND
	  //if(verbose>2) cout << "---- histo deleted: " << nameHisto << endl; // ND
	}

	}// end loop:histos
      }  // end loop:variables
    }    // end loop:cuts
  }      // end loop: selections
  //////////////////////////

  // END //
  if(_outfile) _outfile->Write();
  if(verbose>2) cout << "- End plot()" << endl;
  return 1;
}

TCut XMetAnalysis::defineCut(TString select, TString region)
{

  // Trigger
  TCut trig  = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
  if(_isAN15) trig = "hltmet90>0 || hltmet120>0";

  // Lepton selection depending on the region
  TCut leptons = "";
  if(region.Contains("zctrl")) {
    leptons = "(nelectrons==0 && ntaus==0)";
    leptons *= "(zmass > 60 && zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1))";
  }
  else if(region.Contains("wctrl")) {
    leptons = "(nelectrons==0 && ntaus==0 && nmuons==1)";
    leptons *= "(wmt > 50 && wmt < 100 && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1)";
  }
  else if(region.Contains("zfromw")) {
    //leptons = "nelectrons == 0 && ntaus == 0 && nmuons == 1";
    leptons = "nelectrons == 0";
    leptons *= "(abs(l1id) == 13 || abs(l2id) == 13)";
  }
  else if(region.Contains("signal")) {
    leptons = "nelectrons == 0 && ntaus == 0 && nmuons == 0";
  }

  if(_tag.Contains("NoLepVeto")) leptons = "";

  // Photon veto
  TCut photons="";
  if(_isAN15) photons="nphotons==0";

  // MET
  TCut metID = "(abs(pfmet - calomet) < 2*calomet)" ;
  if(_tag.Contains("NoMetClean"))  metID = "";
  if(_isAN15)        metID = "";
  //
  TCut metCut="";
  const UInt_t nIn=2;
  TString inputs[nIn] = {_tag, select};
  for(UInt_t in=0 ; in<nIn ; in++) {
    if(inputs[in].Contains("NoMetCut"))    metCut = "";
    else if(inputs[in].Contains("Met200")) metCut = "t1mumet>200";
    else if(inputs[in].Contains("Met250")) metCut = "t1mumet>250";
    else if(inputs[in].Contains("Met350")) metCut = "t1mumet>350";
    else if(inputs[in].Contains("MetFrom0to200"))   metCut = "t1mumet<=200";
    else if(inputs[in].Contains("MetFrom200to250")) metCut = "t1mumet>200 && t1mumet<=250";
    else if(inputs[in].Contains("MetFrom250to350")) metCut = "t1mumet>250 && t1mumet<=350";
  }

  // JETS
  TCut jetID1, jetID2, jetID3, jetIDMult, jetIDMono;
  //
  if(_isRun1) { // Run1 cuts
    jetID1 = "(signaljetpt > 110 && abs(signaljeteta) < 2.4 && signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2)";
    jetID2 = "(secondjetpt>30 && abs(secondjeteta)<4.5)";
    jetID3 = "(thirdjetpt>30 && abs(thirdjeteta)<4.5)";
    //
    jetIDMult = "(njets==1 || ( (secondjetpt>30 && abs(secondjeteta)<4.5)&&(njets==2 || (njets==3 && thirdjetpt>30 && abs(thirdjeteta)<4.5) ) ) )" ;
    //
    jetIDMono = "(njets==1 || (njets==2 && secondjetpt>30 && abs(secondjeteta)<4.5) )" ;
  }
  else { // Updated cuts
    jetID1 = "(signaljetpt>30 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1)";
    jetID2 = "(secondjetpt>30 && abs(secondjeteta)<2.5 && secondjetCHfrac > 0.1)";
    jetID3 = "( thirdjetpt>30 &&  abs(thirdjeteta)<2.5 &&  thirdjetCHfrac > 0.1)";
    //
    jetIDMult = "(njets==1 || ( (secondjetpt>30 && abs(secondjeteta)<2.5 && secondjetCHfrac > 0.1) && (njets==2 || (njets==3 && thirdjetpt>30 && abs(thirdjeteta)<2.5 && thirdjetCHfrac  > 0.1) ) ) )" ;
    //
    jetIDMono = "(njets==1 || (njets==2 && secondjetCHfrac>0.1 && secondjetpt>30 && abs(secondjeteta)<2.5) )" ;
  }

  // Jet multiplicity, QCD killer
  TCut jetID ="";
  TCut jetBin="";
  TCut dphi  ="";
  TCut alphat="";
  TCut apcjetmetmax="apcjetmetmax>0.55";

  TCut jmdphi="";
  if(     select.Contains("JetMet0p2"))  jmdphi = "abs(jetmetdphimin)>0.2";
  if(select.Contains("JetMet0p4"))       jmdphi = "abs(jetmetdphimin)>0.4";
  if(select.Contains("JetMet0p45"))      jmdphi = "abs(jetmetdphimin)>0.45";
  if(select.Contains("JetMet0p5"))       jmdphi = "abs(jetmetdphimin)>0.5";
  if(select.Contains("JetMet0p55"))      jmdphi = "abs(jetmetdphimin)>0.55";
  if(select.Contains("JetMet0p6"))       jmdphi = "abs(jetmetdphimin)>0.6";
  if(select.Contains("JetMet0p8"))       jmdphi = "abs(jetmetdphimin)>0.8";
  if(select.Contains("JetMetBelow0p4"))  jmdphi = "abs(jetmetdphimin)<0.4";

  if(     select.Contains("alljets")) {
    jetBin = "njets>=1 && njets<=3";
    jetID  = jetID1*jetIDMult;
    dphi   = "njets==1 || abs(jetjetdphi) < 2.5";
    alphat = "njets==1 || alphat>0.55";
  }
  else if(select.Contains("monojet")) {
    jetBin = "njets>=1 && njets<=2";
    jetID  = jetID1*jetIDMono;
    dphi   = "njets==1 || abs(jetjetdphi) < 2.5";
    alphat = "njets==1 || alphat>0.55";
  }
  else if(select.Contains("1jet")) { 
    jetBin = "njets==1"; 
    jetID  = jetID1; 
    dphi   = "1";
    alphat = "1";
  }
  else if(select.Contains("2jet")) { 
    jetBin = "njets==2"; 
    jetID  = jetID1*jetID2; 
    dphi   = "abs(jetjetdphi) < 2.5";
    alphat = "alphat>0.55";
  }
  else if(select.Contains("3jet")) {
    jetBin = "njets==3";
    jetID  = jetID1*jetID2*jetID3;
    dphi   = "abs(jetjetdphi) < 2.5";
    alphat = "alphat>0.55";
  }
  else if(select.Contains("4jet")) {
    jetBin = "njets>=4";
    jetID  = jetID1*jetID2*jetID3;
    dphi   = "abs(jetjetdphi) < 2.5";
    alphat = "alphat>0.55";
  }

  // QCD Killer
  TCut noqcd="";
  //
  for(UInt_t in=0 ; in<nIn ; in++) {
    if(inputs[in].Contains("dphi"))          noqcd *= dphi;
    if(inputs[in].Contains("alphat"))        noqcd *= alphat;
    if(inputs[in].Contains("apcjetmetmax"))  noqcd *= apcjetmetmax;
    if(inputs[in].Contains("JetMet"))        noqcd *= jmdphi;
  }

  // Noise cleaning
  TCut noise="(flaghbheloose>0 && flagcsctight>0 && flageebadsc>0)";

  // Max run
  TCut maxrun="";
  if(_tag.Contains("MaxRun256730")) maxrun="(run<256730)";
  
  return (trig*noise*maxrun*metCut*metID*leptons*photons*jetID*jetBin*noqcd);

  // "(hltmet90>0 || hltmet120>0) && t1mumet>200 && nelectrons == 0 && ntaus == 0 && nmuons == 0 && nphotons==0 signaljetpt>30 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && njets==1"

}

Int_t XMetAnalysis::DefineChainsAN15()
{

  if(verbose>1) cout << "- begin DefineChainsAN15()" << endl;

  //_pathMC   = "/user/ndaci/Data/XMET/AdishTrees_15Oct2015/Spring15MC_25ns/";
  //_pathData = "/user/ndaci/Data/XMET/AdishTrees_15Oct2015/Run2015D/";
  //_pathMC   = "/user/ndaci/Data/XMET/NadirTrees_22Oct2015/Spring15MC_25ns/";
  //_pathData = "/user/ndaci/Data/XMET/NadirTrees_22Oct2015/Run2015D/";
  _pathMC   = "/user/ndaci/Data/XMET/NadirTrees_26Oct2015/Spring15MC_25ns/";
  _pathData = "/user/ndaci/Data/XMET/NadirTrees_26Oct2015/Run2015D/";

  _lumi  = 0.210; // fixme: we have 210 /pb up to run 257599
  //_lumi    = 0.552672886226; // 2015D 05Oct2015 JSON 246908-258750
  _rescale = 1.0; 
  _qcdScale= 1.0; // fixme: will need update
  _useLO   = false; // apply k-factors to LO Z/W samples // fixme

  // // Data
  _mapProcess["data_met"] = XMetProcess("data_met", kBlack, "skimJSON.root");
  /*
  _mapProcess["data_2m" ] = XMetProcess("data_2m",  kBlack, "skimJSON.root");
  _mapProcess["data_1m" ] = XMetProcess("data_1m",  kBlack, "skimJSON.root");
  _mapProcess["data_1ph"] = XMetProcess("data_1ph", kBlack, "skimJSON.root");
  _mapProcess["data_2e" ] = XMetProcess("data_2e",  kBlack, "skimJSON.root");
  */
  
  // MC backgrounds
  _mapProcess["znn"  ] = XMetProcess("znn",   kAzure+7,  "skimMuMet200.root");
  _mapProcess["zll"  ] = XMetProcess("zll",   kPink+9,   "skimMuMet200.root");
  _mapProcess["wln"  ] = XMetProcess("wln",   kGreen+2,  "skimMuMet200.root");
  _mapProcess["ttbar"] = XMetProcess("ttbar", kMagenta+3,"skimMuMet200.root");
  _mapProcess["qcd"  ] = XMetProcess("qcd",   kRed,      "skimMuMet200.root");
  _mapProcess["vv"   ] = XMetProcess("vv",    kBlue+1,   "skimMuMet200.root");
  _mapProcess["top"  ] = XMetProcess("top",   kOrange-3, "skimMuMet200.root");

  // // Data
  // _mapProcess["data_met"] = XMetProcess("data_met", kBlack, "skimMumet100WgtSum.root");
  // _mapProcess["data_2m" ] = XMetProcess("data_2m",  kBlack, "skimMumet100WgtSum.root");
  // _mapProcess["data_1m" ] = XMetProcess("data_1m",  kBlack, "skimMumet100WgtSum.root");
  // _mapProcess["data_1ph"] = XMetProcess("data_1ph", kBlack, "skimMumet100WgtSum.root");
  // _mapProcess["data_2e" ] = XMetProcess("data_2e",  kBlack, "skimMumet100WgtSum.root");
  
  // MC backgrounds
  // _mapProcess["znn"  ] = XMetProcess("znn",   kAzure+7,  "skimMumet100WgtSum.root");
  // _mapProcess["zll"  ] = XMetProcess("zll",   kPink+9,   "skimMumet100WgtSum.root");
  // _mapProcess["wln"  ] = XMetProcess("wln",   kGreen+2,  "skimMumet100WgtSum.root");
  // _mapProcess["ttbar"] = XMetProcess("ttbar", kMagenta+3,"skimMumet100WgtSum.root");
  // _mapProcess["qcd"  ] = XMetProcess("qcd",   kRed,      "skimMumet100WgtSum.root");
  // _mapProcess["vv"   ] = XMetProcess("vv",    kBlue+1,   "skimMumet100WgtSum.root");
  // _mapProcess["top"  ] = XMetProcess("top",   kOrange-3, "skimMumet100WgtSum.root");

  // Sub-directories in _path
  //
  // data
  /*
  _mapProcess["data_met"].AddDir("met");
  _mapProcess["data_1ph"].AddDir("singleph");
  _mapProcess["data_2e" ].AddDir("doubleel");
  _mapProcess["data_1m" ].AddDir("singlemu");
  _mapProcess["data_2m" ].AddDir("doublemu");
  */

  _mapProcess["data_met"].AddDir("met_05Oct2015");
  _mapProcess["data_1ph"].AddDir("singleph_05Oct2015");
  _mapProcess["data_2e" ].AddDir("doubleel_05Oct2015");
  _mapProcess["data_1m" ].AddDir("singlemu_05Oct2015");
  _mapProcess["data_2m" ].AddDir("doublemu_05Oct2015");

  /*
  _mapProcess["data_met"].AddDir("met_PRV4");
  _mapProcess["data_1ph"].AddDir("singleph_PRV4");
  _mapProcess["data_2e" ].AddDir("doubleel_PRV4");
  _mapProcess["data_1m" ].AddDir("singlemu_PRV4");
  _mapProcess["data_2m" ].AddDir("doublemu_PRV4");
  */

  //
  // mc backgrounds
  //
  _mapProcess["znn"].AddDir("znn100to200");
  _mapProcess["znn"].AddDir("znn200to400");
  _mapProcess["znn"].AddDir("znn400to600");
  _mapProcess["znn"].AddDir("znn600toinf");
  //
  _mapProcess["wln"].AddDir("wln100to200");
  _mapProcess["wln"].AddDir("wln200to400");
  _mapProcess["wln"].AddDir("wln400to600");
  _mapProcess["wln"].AddDir("wln600toinf");
  //_mapProcess["wln"].AddDir("wln"); //fixme: NLO sample
  //
  _mapProcess["zll"].AddDir("zll100to200");
  _mapProcess["zll"].AddDir("zll200to400");
  _mapProcess["zll"].AddDir("zll400to600");
  _mapProcess["zll"].AddDir("zll600toinf");
  //_mapProcess["zll"].AddDir("zll"); //fixme: NLO sample
  //
  _mapProcess["ttbar"].AddDir("ttbar");
  //
  _mapProcess["top"].AddDir("singletbart");
  _mapProcess["top"].AddDir("singletbarw");
  _mapProcess["top"].AddDir("singlett");
  _mapProcess["top"].AddDir("singletw");
  //
  _mapProcess["qcd"].AddDir("qcdht1000to1500");
  _mapProcess["qcd"].AddDir("qcdht100to200");
  _mapProcess["qcd"].AddDir("qcdht1500to2000");
  _mapProcess["qcd"].AddDir("qcdht2000toinf");
  _mapProcess["qcd"].AddDir("qcdht200to300");
  _mapProcess["qcd"].AddDir("qcdht300to500");
  _mapProcess["qcd"].AddDir("qcdht500to700");
  _mapProcess["qcd"].AddDir("qcdht700to1000");
  //
  _mapProcess["vv"].AddDir("ww");
  _mapProcess["vv"].AddDir("wz");
  _mapProcess["vv"].AddDir("zz");
  //

  if(verbose>1) cout << "#entries in mapProcess : " << _mapProcess.size() << endl;

  // Add the files to the chains
  for( _itProcess=_mapProcess.begin() ; _itProcess!=_mapProcess.end() ; _itProcess++ ) {

    _itProcess->second.SetNameTree("tree/tree");
    _itProcess->second.SetStyle(kFullCircle);
    _itProcess->second.SetSize(1.0);

    if(_itProcess->first.Contains("data")) 
      _itProcess->second.SetPath(_pathData);
    else                          
      _itProcess->second.SetPath(_pathMC);

    _itProcess->second.AddTrees();
  }

  return 0;
}

Int_t XMetAnalysis::DefineChainsRun1()
{

  if(verbose>1) cout << "- begin DefineChainsRun1()" << endl;

  _pathData= "/user/ndaci/Data/XMET/5Xtrees/";
  _pathMC  = "/user/ndaci/Data/XMET/5Xtrees/";
  _lumi    = 19.7; // Adish
  //_lumi    = 19.5; // Run1 AN
  _rescale = 1.0;
  _qcdScale= 1.6;
  
  // MC backgrounds
  _mapProcess["znn"  ] = XMetProcess("znn",   kAzure+7,  "reducedtree.root");
  _mapProcess["zll"  ] = XMetProcess("zll",   kPink+9,   "reducedtree.root");
  _mapProcess["wln"  ] = XMetProcess("wln", kGreen+2,  "reducedtree.root");
  _mapProcess["ttbar"] = XMetProcess("ttbar", kMagenta+3,"reducedtree.root");
  _mapProcess["qcd"  ] = XMetProcess("qcd",   kRed,      "reducedtree.root");
  _mapProcess["vv"   ] = XMetProcess("vv",    kBlue+1,   "reducedtree.root");
  _mapProcess["top"  ] = XMetProcess("top",   kOrange-3, "reducedtree.root");

  // Data driven backgrounds
  /// no corr
  _mapProcess["zn_cr_DA"] = XMetProcess("zn_cr_DA",kAzure+7,"reducedtree.root");
  _mapProcess["zn_cr_MC"] = XMetProcess("zn_cr_MC",kAzure+7,"reducedtree.root");
  _mapProcess["wj_cr_DA" ] = XMetProcess("wj_cr_DA", kGreen+2, "reducedtree.root");
  _mapProcess["wj_cr_MC" ] = XMetProcess("wj_cr_MC", kGreen+2, "reducedtree.root");
  /// corr
  _mapProcess["zn_cr_corr_DA"] = XMetProcess("zn_cr_corr_DA",kAzure+7,"reducedtree.root");
  _mapProcess["zn_cr_corr_MC"] = XMetProcess("zn_cr_corr_MC",kAzure+7,"reducedtree.root");
  _mapProcess["wj_cr_corr_DA" ] = XMetProcess("wj_cr_corr_DA", kGreen+2, "reducedtree.root");
  _mapProcess["wj_cr_corr_MC" ] = XMetProcess("wj_cr_corr_MC", kGreen+2, "reducedtree.root");

  // Signal
  const UInt_t nSpin=2;
  const UInt_t nMass=9;
  TString nameSpin[nSpin] = {"V","AV"};
  TString nameMass[nMass] = {"0p1","1","10","100","200","300","400","700","1000"};
  // FIXME
  /*
  Double_t dmXSec[nMass][nSpin] 
    = { {,},
	{,},
	{,},
	{,},
	{,},
	{,},
	{,},
	{,},
	{,} };
  //
  Int_t dmNGen[nMass][nSpin] 
    = { {64436,65480},
	{40433,43034},
	{44266,39712},
	{41514,43225},
	{48028,43699},
	{45040,47848},
	{45546,47770},
	{52487,46773},
	{44781,46569} };
  */
  //
  TString nameProcess;
  for(UInt_t iS=0 ; iS<nSpin ; iS++) {
    for(UInt_t iM=0 ; iM<nMass ; iM++) {
      nameProcess = "dm_"+nameSpin[iS]+"_"+nameMass[iM];
      _mapProcess[nameProcess] = XMetProcess(nameProcess, kBlack, "tree_"+nameProcess+".root");
      _mapProcess[nameProcess].AddDir("signal");
      //_mapProcess[nameProcess].SetXSec(dmXSec[iM][iS]); // FIXME
      //_mapProcess[nameProcess].SetNGen(dmNGen[iM][iS]);
    }
  }

  // Data
  _mapProcess["data"]   = XMetProcess("data",   kBlack, "reducedtree.root");

  // Sub-directories in _path
  /// mc backgrounds
  _mapProcess["znn"].AddDir("znn");
  _mapProcess["zll"].AddDir("zll");
  _mapProcess["wln"].AddDir("wln");
  _mapProcess["ttbar"].AddDir("ttbar");
  _mapProcess["top"].AddDir("singletop");
  _mapProcess["qcd"].AddDir("qcd");
  _mapProcess["vv"].AddDir("dibosons");
  /// data-driven (no corr)
  _mapProcess["zn_cr_DA"].AddDir("met");
  _mapProcess["zn_cr_MC"].AddDir("zll");
  _mapProcess["wj_cr_DA"].AddDir("met");
  _mapProcess["wj_cr_MC"].AddDir("wln");
  /// data-driven (corr)
  _mapProcess["zn_cr_corr_DA"].AddDir("met");
  _mapProcess["zn_cr_corr_MC"].AddDir("bkgz");
  _mapProcess["wj_cr_corr_DA"].AddDir("met");
  _mapProcess["wj_cr_corr_MC"].AddDir("bkgw");
  /// signal still missing
  /// data
  _mapProcess["data"].AddDir("met");

  if(verbose>1) cout << "#entries in mapProcess : " << _mapProcess.size() << endl;
  
  // Add the files to the chains
  for( _itProcess=_mapProcess.begin() ; _itProcess!=_mapProcess.end() ; _itProcess++ ) {
    _itProcess->second.SetNameTree("tree/tree");
    if(_itProcess->first=="data") _itProcess->second.SetPath(_pathData);
    else                          _itProcess->second.SetPath(_pathMC);
    _itProcess->second.AddTrees();
  }

  if(verbose>1) cout << "- end DefineChains()" << endl;

  return 0;
}

Int_t XMetAnalysis::GetHistos(const UInt_t nS  , TString *select, 
			      const UInt_t nCut, TString *scanCut,
			      const UInt_t nV  , TString *var, 
			      vector<TString> locProcesses)
{

  // Open histos root file
  if(_outfile) _outfile->cd();
  else {
    cout << "ERROR: outfile undefined! Exit ==> []" << endl;
    return -2;
  }

  const UInt_t nP = locProcesses.size();
  TString name;
  TH1F* hTemp;

  for(UInt_t iS=0 ; iS<nS ; iS++) {       // loop: selections
    for(UInt_t iC=0 ; iC<nCut ; iC++) {   // loop: cuts
      for(UInt_t iV=0 ; iV<nV ; iV++) {   // loop: variables
	for(UInt_t iP=0 ; iP<nP ; iP++) { // loop: processes

	  name = "h_"+var[iV]+"_"+locProcesses[iP]+"_"+select[iS]+"_"+scanCut[iC];
	  hTemp = (TH1F*) _outfile->Get(name);

	  if(verbose>1) {
	    cout << "----- try to retrieve [" << name << "]: ";
	    if(!hTemp) {
	      cout << "FAIL" << endl;
	    }
	    else {
	      cout << "SUCCESS => Integral=" 
		   << hTemp->Integral() << endl;
	    }
	  }

	  _mapHistos[select[iS]][locProcesses[iP]][scanCut[iC]][var[iV]] = hTemp;
	}
      }
    }
  }

  return 0;
}

Int_t XMetAnalysis::DrawStackPlots(const UInt_t nS  , TString* select, 
				   const UInt_t nCut, TString* scanCut,
				   const UInt_t nV  , TString* var,
				   vector<TString> locProcesses)
{

  const UInt_t nP  = locProcesses.size();

  TH1F *hTemp, *hData, *hTotB, *hRatio;
  THStack* stack;
  TLegend* leg;

  TString name, title, modeLeg;
  Int_t color, style, size;
  Bool_t fill;
  
  // Produce stack plots
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iC=0 ; iC<nCut ; iC++) {
      for(UInt_t iV=0 ; iV<nV ; iV++) {

	TCanvas c("c","c",20,20,600,600);
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	name  = "ths_"   +select[iS]+"_"+scanCut[iC]+"_"+var[iV];
	title = "Stack: "+select[iS]+" "+scanCut[iC]+" "+var[iV];
	stack = new THStack(name,title);

	leg = new TLegend(0.7,0.7,0.89,0.89,"","brNDC");
	setStyle(leg);

	// Loop over processes
	for(UInt_t iP=0; iP<nP ; iP++) {

	  hTemp = _mapHistos[select[iS]][locProcesses[iP]][scanCut[iC]][var[iV]];
	  if(!hTemp) continue;

	  hTemp->SetTitle(_Title[select[iS]]);

	  // Decide fill (stack) or point (non stack)
	  if(!locProcesses[iP].Contains("data")) {
	    modeLeg = "F";
	    fill    = kTRUE;
	  }
	  else {
	    modeLeg = "P";
	    fill    = kFALSE;
	  }

	  // Get XMetProcess attributes
	  color = _mapProcess[locProcesses[iP]].GetColor();
	  style = _mapProcess[locProcesses[iP]].GetStyle();
	  size  = _mapProcess[locProcesses[iP]].GetSize();

	  // Apply attributes
	  setStyle(hTemp,color,style,size,kFALSE,fill);
	  leg->AddEntry(hTemp, locProcesses[iP],modeLeg);

	  // Debug printouts
	  if(verbose>2) {
	    cout << "----- stack-get: " << hTemp->GetName() 
		 << "->Integral()="     << hTemp->Integral() 
		 << endl;
	  }

	  // Add to stack
	  if(!locProcesses[iP].Contains("data")) 
	    stack->Add(hTemp);
	  else {
	    hTemp->Draw("P");  // set max with data
	  }
	}

	// Finalize plot
	stack->Draw("HISTSAME"); // plot stack
	//
	for(UInt_t iP=0; iP<nP ; iP++) {
	  hTemp = _mapHistos[select[iS]][locProcesses[iP]][scanCut[iC]][var[iV]];
	  if(!hTemp) continue;
	  if(locProcesses[iP].Contains("data")) hTemp->Draw("PSAME"); // superimpose back data on stack
	}

	// Finalize plot
	leg->Draw();
	c.Print(_dirOut+"/"+_tag+"/"+name+".pdf","pdf");

      }
    }
  }  

  return 0;
}
