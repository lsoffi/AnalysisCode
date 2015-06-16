#include "XMetAnalysis.h"

using namespace std;
Int_t verbose = 2;

XMetAnalysis::XMetAnalysis(TString tag)
{
  _tag = tag;

  _path    = "/user/ndaci/Data/XMET/5Xtrees/";
  //_path    = "/user/ndaci/Data/XMET/7Xtrees/";

  _lumi    = 19.7;
  _rescale = 1.0;
  _outfile = new TFile("plots/"+_tag+"/plots_"+_tag+".root","recreate");
  _outlog  = new ofstream("plots/"+_tag+"/yields_"+_tag+".txt",ios::out);

  DefineChains();
}

XMetAnalysis::~XMetAnalysis()
{
  _outfile->Close();
  delete _outfile;
}

Int_t XMetAnalysis::Analysis()
{

  // Processes to use
  vector<TString> locProcesses;

  // MC backgrounds
  locProcesses.push_back("znn"); //labelProc.push_back("Z(#nu#nu)");
  locProcesses.push_back("zll"); //labelProc.push_back("Z(ll)");
  locProcesses.push_back("wjets"); //labelProc.push_back("W(l#nu)");
  locProcesses.push_back("ttbar"); //labelProc.push_back("t#bar{t}");
  locProcesses.push_back("top"); //labelProc.push_back("t");
  locProcesses.push_back("vv"); //labelProc.push_back("VV");
  locProcesses.push_back("qcd"); //labelProc.push_back("QCD");

  // Data-driven
  /// without acc, eff, BR corrections
  locProcesses.push_back("zn_cr_DA"); //labelProc.push_back("Z(#nu#nu) CR Data");
  locProcesses.push_back("zn_cr_MC"); //labelProc.push_back("Z(#nu#nu) CR MC");
  locProcesses.push_back("zn_sr_DD"); //labelProc.push_back("Z(#nu#nu) SR DD");
  locProcesses.push_back("wj_cr_DA"); //labelProc.push_back("W(l#nu) CR Data");
  locProcesses.push_back("wj_cr_MC"); //labelProc.push_back("W(l#nu) CR MC");
  locProcesses.push_back("wj_sr_DD"); //labelProc.push_back("W(l#nu) SR DD");
  /// with acc, eff, BR corrections
  locProcesses.push_back("zn_cr_corr_DA"); //labelProc.push_back("Z(#nu#nu) CR Data (SF)");
  locProcesses.push_back("zn_cr_corr_MC"); //labelProc.push_back("Z(#nu#nu) CR MC (SF)");
  locProcesses.push_back("zn_sr_corr_DD"); //labelProc.push_back("Z(#nu#nu) SR DD (SF)");
  locProcesses.push_back("wj_cr_corr_DA"); //labelProc.push_back("W(l#nu) CR Data (SF)");
  locProcesses.push_back("wj_cr_corr_MC"); //labelProc.push_back("W(l#nu) CR MC (SF)");
  locProcesses.push_back("wj_sr_corr_DD"); //labelProc.push_back("W(l#nu) SR DD (SF)");

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
      //labelProc.push_back("DM "+nameSpin[iS]+" M="+nameMass[iM]);
    }
  }

  // Data
  locProcesses.push_back("data"); //labelProc.push_back("Data");

  // Selections and variables
  const UInt_t nS=5;
  const UInt_t nV=1;
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"mumet"};

  UInt_t  nBins[nV]  = {50};
  Float_t xFirst[nV] = {0};
  Float_t xLast[nV]  = {1000};

  for(UInt_t iS=0 ; iS<nS ; iS++) {
    //for(UInt_t iS=0 ; iS<1 ; iS++) { // FIXME

    if(verbose>1) cout << "- selection : " << select[iS] << endl;

    //select, nV, var, nBins, xFirst, xLast, vector<TString> locProcesses, vector<TString> labelProc
    plot(select[iS], nV, var, nBins, xFirst, xLast, locProcesses);

  }

  return 0;
}

Int_t XMetAnalysis::StudyQCDKiller()
{

  // STYLE
  //gROOT->Reset();
  //loadPresentationStyle();
  //gROOT->ForceStyle();

  // Processes to use
  vector<TString> locProcesses;
  vector<TString> labelProc;
  locProcesses.push_back("znn"); labelProc.push_back("Z(#nu#nu)");
  locProcesses.push_back("qcd"); labelProc.push_back("QCD");

  // Selections and variables
  const UInt_t nS=5;
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};

  const UInt_t nV=7;
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin",
			"jetjetdphi","jetmetdphimin",
			"dphiJ1J3","dphiJ2J3"};
  //
  //UInt_t  nBins[nV]  = {40, 50, 50, 50,  50, 50, 50};
  //UInt_t  nBins[nV]  = {400, 500, 500, 500, 500, 500, 500};
  //UInt_t  nBins[nV]  = {800, 1000, 1000, 1000,  1000, 1000, 1000};
  //UInt_t  nBins[nV]  = {4000, 5000, 5000, 5000,  5000, 5000, 5000};
  UInt_t  nBins[nV]  = {8000, 10000, 10000, 10000, 10000, 10000, 10000};
  //
  Float_t xFirst[nV] = {0,  0,  0,  0,   0  , 0  , 0};
  Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2, 3.2, 3.2};

  // Produce 1 plot per {selection ; variable}
  for(UInt_t iS=0 ; iS<nS ; iS++) {

    if(verbose>1) cout << "- selection : " << select[iS] << endl;

    plot(select[iS], nV, var, nBins, xFirst, xLast, locProcesses);

  }

  _outfile->Close();
  return 0;
}

Int_t XMetAnalysis::plot(TString select, const UInt_t nV, TString* var,
			 UInt_t* nBins, Float_t* xFirst, Float_t* xLast, 
			 vector<TString> locProcesses)
{

  (*_outlog) << "Selection: " << select << endl;

  // Define selections //
  TCut weight="";
  TString region = "signal";
  TCut cut = defineCut(select, region);

  // Declare histograms //
  map<TString, map<TString,TH1F*> > mapVarHistos;
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax;

  const UInt_t nP = locProcesses.size();
  Float_t integral=0;

  // Loop over chains and generate histograms //
  //
  TString nameDir, locVar, variable;
  Int_t color;
  TH1F *hTemp, *hSRDD, *hTemp_cr_d, *hTemp_cr_mc;
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
    else weight = "puwgt*wgt";

    // Choose the phase space region
    region = "signal";
    if(     nameDir.Contains("zn_cr")) {
      if(nameDir.Contains("corr"))
	region = "zctrl_corr";
      else
	region = "zctrl";
    }
    else if(nameDir.Contains("wj_cr"))  {
      if(nameDir.Contains("corr"))
	region = "wctrl_corr";
      else
	region = "wctrl";
    }
    else {region = "signal";}
    cut = defineCut(select, region);
    
    // Skim the chain
    _mapProcess[nameDir].Skim(select, cut);

    // Loop over requested variables
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      if(verbose>1) cout << "--- variable: " << var[iV] << endl;

      // Define histogram and set style
      mapVarHistos[nameDir][var[iV]] = new TH1F("h_"+var[iV]+"_"+nameDir+"_"+select, 
						var[iV]+" "+nameDir+" "+select,
						nBins[iV], xFirst[iV], xLast[iV]);
      hTemp = mapVarHistos[nameDir][var[iV]];
      color = _mapProcess[nameDir].GetColor();
      setStyle( hTemp , color );

      // Draw the variable
      locVar = var[iV];
      if(var[iV].Contains("phi"))  locVar = "abs("+var[iV]+")";
      if(     var[iV]=="dphiJ1J3") locVar = "abs(signaljetphi-thirdjetphi)";
      else if(var[iV]=="dphiJ2J3") locVar = "abs(secondjetphi-thirdjetphi)";

      _mapProcess[nameDir].Draw(hTemp, locVar, cut, weight);
      if(verbose>1) cout << "--- drew variable: " << locVar << endl;

      // Checks
      if(verbose>1) {
	cout << "--- check:" 
	     << " nameDir="  << nameDir
	     << " GetName="  << hTemp->GetName()
	     << " Entries="  << hTemp->GetEntries()
	     << " Integral=" << hTemp->Integral()
	     << " Weight="   << weight
	     << endl;
      }

      // Normalize
      if( !nameDir.Contains("met") ) hTemp->Scale(_lumi*_rescale);

      // Determine extrema
      locMin = hTemp->GetMinimum();
      locMax = hTemp->GetMaximum();
      if(locMin<minPlot) minPlot = locMin;
      if(locMax>maxPlot) maxPlot = locMax;
    } // end loop:variables
  } // end loop:processes 
  if(verbose>1) cout << "-- end loop over processes" << endl;
  
  // Produce sr_DD histograms
  if(verbose>1) cout << "-- Ready to produce sr_DD histograms" << endl;
  const UInt_t nDD=2;
  const UInt_t nCorr=2;
  TString myDD[nDD]     = {"zn","wj"};
  TString myCorr[nCorr] = {"_","_corr_"};
  //
  for(UInt_t iDD=0 ; iDD<nDD ; iDD++) {
    
    // Check that the current process (znn or wj) has been requested initially
    bool srddIsHere=false;
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      if(locProcesses[iP].Contains(myDD[iDD]+"_cr")) {
	srddIsHere=true;
	break; // hammertime
      }
    }
    if(!srddIsHere) {
      if(verbose>1) cout << "--- unrequested process: " << myDD[iDD] << endl;
      continue;
    }

    for(UInt_t iCorr=0 ; iCorr<nCorr ; iCorr++) {
      for(UInt_t iV=0 ; iV<nV ; iV++) {
	//
	if(verbose>1) cout << "----- produce: " << "h_"+var[iV]+"_"+myDD[iDD]+"_sr"+myCorr[iCorr]+"DD_"+select << endl;
	//
	hTemp_cr_d  = mapVarHistos[myDD[iDD]+"_cr"+myCorr[iCorr]+"DA"][var[iV]];
	hTemp_cr_mc = mapVarHistos[myDD[iDD]+"_cr"+myCorr[iCorr]+"MC"][var[iV]];
	if(!hTemp_cr_d || !hTemp_cr_mc) continue;
	//
	hSRDD = (TH1F*)hTemp_cr_d->Clone("h_"+var[iV]+"_"+myDD[iDD]+"_sr"+myCorr[iCorr]+"DD_"+select);
	if(hTemp_cr_mc) hSRDD->Add( hTemp_cr_mc , -1.0 );
	mapVarHistos[myDD[iDD]+"_sr"+myCorr[iCorr]+"DD"][var[iV]] = hSRDD;
      }
    }
  }

  //////////////////////////

  // Yields outlog //

  // Loop over variables
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    // Print out variable and processes names
    (*_outlog) << "Var: " << var[iV] << endl;
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      (*_outlog) << setw(10) << locProcesses[iP];
    }
    (*_outlog) << endl;
    //  
    // Loop over processes
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      // Get histogram
      hTemp = mapVarHistos[locProcesses[iP]][var[iV]];
      if(!hTemp) {
	if(verbose>1) {
	  cout << "ERROR : yield loop did not find histo:"
	       << locProcesses[iP]+" "+var[iV]
	       << endl;
	}
	(*_outlog) << setw(10) << -999999;
	continue;
      }
      else {
	integral = hTemp->Integral();
	(*_outlog) << setw(10) << integral;
      }
    } // end loop over processes
    //
    (*_outlog) << endl;
  } // end loop over var
  //
  (*_outlog) << endl;

  //////////////////////////

  // Save histograms //
  _outfile->cd();
  if(verbose>1) cout << "-- _outfile->cd()... done" << endl;
  TString nameHisto;
  map<TString, map<TString,TH1F*> >::iterator itVarHistos;
  map<TString,TH1F*>::iterator itHistos;
  //
  for( itVarHistos=mapVarHistos.begin() ; itVarHistos!=mapVarHistos.end() ; itVarHistos++) {
    for( itHistos=itVarHistos->second.begin() ; itHistos!=itVarHistos->second.end() ; itHistos++) {

      hTemp = itHistos->second;

      if(hTemp) {
	nameHisto = hTemp->GetName();
	hTemp->Write();
	if(verbose>1) cout << "---- histo written: " << nameHisto << endl;
	delete hTemp;
	if(verbose>1) cout << "---- histo deleted: " << nameHisto << endl;	
      }
    } // end loop over histos in itVarHistos->second
  } // end loop over maps in mapVarHistos

  //////////////////////////

  // END //
  return 1;
}

TCut XMetAnalysis::defineCut(TString select, TString region)
{

  // Trigger
  TCut trig  = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";

  // Lepton selection depending on the region
  TCut leptons = "";
  if(region.Contains("zctrl")) {
    leptons = "(nelectrons==0 && ntaus==0)";
    leptons *= "(zmass > 60 && zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1))";
    if(region.Contains("corr")) leptons *= "((5.942 * 1.023) / (0.79*(1.0-TMath::Exp(-0.00910276*(mumet-36.1669)))))";
  }
  else if(region.Contains("wctrl")) {
    leptons = "(nelectrons==0 && ntaus==0 && nmuons==1)";
    leptons *= "(wmt > 50 && wmt < 100 && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1)";
    if(region.Contains("corr")) leptons *= "(38.7823/TMath::Power(-85.7023 + mumet, 0.667232))";
  }
  else {
    leptons = "(nmuons == 0 && nelectrons == 0 && ntaus == 0)";
  }

  // MET
  TCut metID = "(abs(pfmet - calomet) < 2*calomet)" ;
  TCut metCut="";
  if(_tag.Contains("NoMetCut"))    metCut = "";
  else if(_tag.Contains("Met200")) metCut = "mumet>200";
  else if(_tag.Contains("Met350")) metCut = "mumet>350";
  else if(_tag.Contains("MetFrom0to200"))   metCut = "mumet<=200";
  else if(_tag.Contains("MetFrom200to250")) metCut = "mumet>200 && mumet<=250";
  else if(_tag.Contains("MetFrom250to350")) metCut = "mumet>250 && mumet<=350";

  // JETS
  TCut jetID1 = "(signaljetpt > 110 && abs(signaljeteta) < 2.4 && signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2)";
  TCut jetID2 = "(secondjetpt>30 && abs(secondjeteta)<2.4 && secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)";
  TCut jetID3 = "(thirdjetpt>30 && abs(thirdjeteta)<4.5 && thirdjetNHfrac  < 0.7 && thirdjetEMfrac  < 0.9)";
  //
  TCut jetIDMult = "(njets==1 || ( (secondjetNHfrac<0.7 && secondjetEMfrac<0.9 && secondjetpt>30 && abs(secondjeteta)<2.5)&&(njets==2 || (njets==3 && thirdjetpt>30 && abs(thirdjeteta)<4.5 && thirdjetNHfrac  < 0.7 && thirdjetEMfrac  < 0.9) ) ) )" ;
  //
  TCut jetIDMono = "(njets==1 || (njets==2 && secondjetNHfrac<0.7 && secondjetEMfrac<0.9 && secondjetpt>30 && abs(secondjeteta)<2.5) )" ;

  // Jet multiplicity, QCD killer
  TCut jetID="";
  TCut jetBin="";
  TCut dphi="";
  TCut alphat="";
  TCut apcjetmetmax="apcjetmetmax>0.55";

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

  TCut noqcd="run>-999";
  if(     select.Contains("dphi"))         noqcd = dphi;
  else if(select.Contains("alphat"))       noqcd = alphat;
  else if(select.Contains("apcjetmetmax")) noqcd = apcjetmetmax;

  return (trig*leptons*metID*metCut*noqcd*jetID*jetBin);

}

Int_t XMetAnalysis::DefineChains()
{
  if(verbose>1) cout << "- begin DefineChains()" << endl;

  // MC backgrounds
  _mapProcess["znn"  ] = XMetProcess("znn",   kAzure+7,  "reducedtree.root");
  _mapProcess["zll"  ] = XMetProcess("zll",   kPink+9,   "reducedtree.root");
  _mapProcess["wjets"] = XMetProcess("wjets", kGreen+2,  "reducedtree.root");
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
  TString nameProcess;
  for(UInt_t iS=0 ; iS<nSpin ; iS++) {
    for(UInt_t iM=0 ; iM<nMass ; iM++) {
      nameProcess = "dm_"+nameSpin[iS]+"_"+nameMass[iM];
      _mapProcess[nameProcess] = XMetProcess(nameProcess, kBlack, "tree_"+nameProcess+".root");
      _mapProcess[nameProcess].AddDir("signal");
    }
  }

  // Data
  _mapProcess["data"]   = XMetProcess("data",   kBlack, "reducedtree.root");

  // Sub-directories in _path
  /// mc backgrounds
  _mapProcess["znn"].AddDir("znn");
  _mapProcess["zll"].AddDir("zll");
  _mapProcess["wjets"].AddDir("wjets");
  _mapProcess["ttbar"].AddDir("ttbar");
  _mapProcess["top"].AddDir("singletop");
  _mapProcess["qcd"].AddDir("qcd");
  _mapProcess["vv"].AddDir("dibosons");
  /// data-driven (no corr)
  _mapProcess["zn_cr_DA"].AddDir("met");
  _mapProcess["zn_cr_MC"].AddDir("bkgz");
  _mapProcess["wj_cr_DA"].AddDir("met");
  _mapProcess["wj_cr_MC"].AddDir("bkgw");
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
    //_itProcess->second.SetNameFile("reducedtree.root");
    _itProcess->second.SetPath(_path);
    _itProcess->second.AddTrees();
  }

  if(verbose>1) cout << "- end DefineChains()" << endl;

  return 0;
}
