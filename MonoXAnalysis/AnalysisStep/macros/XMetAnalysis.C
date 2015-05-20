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
  vector<TString> labelProc;
  locProcesses.push_back("znn"); labelProc.push_back("Z(#nu#nu)");
  locProcesses.push_back("zll"); labelProc.push_back("Z(ll)");
  locProcesses.push_back("wjets"); labelProc.push_back("W(l#nu)");
  locProcesses.push_back("ttbar"); labelProc.push_back("t#bar{t}");
  locProcesses.push_back("top"); labelProc.push_back("t");
  locProcesses.push_back("vv"); labelProc.push_back("VV");
  locProcesses.push_back("qcd"); labelProc.push_back("QCD");

  locProcesses.push_back("znn_cr_DA"); labelProc.push_back("Z(#nu#nu) CR Data");
  locProcesses.push_back("znn_cr_MC"); labelProc.push_back("Z(#nu#nu) CR MC");
  locProcesses.push_back("wj_cr_DA"); labelProc.push_back("W(l#nu) CR Data");
  locProcesses.push_back("wj_cr_MC"); labelProc.push_back("W(l#nu) CR MC");
  locProcesses.push_back("data"); labelProc.push_back("Data");

  // Selections and variables
  const UInt_t nS=5;
  const UInt_t nV=1;
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"mumet"};

  UInt_t  nBins[nV]  = {50};
  Float_t xFirst[nV] = {0};
  Float_t xLast[nV]  = {1000};

  for(UInt_t iS=0 ; iS<nS ; iS++) {

    if(verbose>1) cout << "- selection : " << select[iS] << endl;

    plot(select[iS], nV, var, nBins, xFirst, xLast, false, false, false, locProcesses, labelProc);

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
  const UInt_t nV=5;
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};
  //UInt_t  nBins[nV]  = {40, 50, 50, 50,  50};
  //UInt_t  nBins[nV]  = {400, 500, 500, 500, 500};
  UInt_t  nBins[nV]  = {800, 1000, 1000, 1000,  1000};
  //UInt_t  nBins[nV]  = {4000, 5000, 5000, 5000,  5000};
  //UInt_t  nBins[nV]  = {8000, 10000, 10000, 10000,  10000};
  Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  /*
  const UInt_t nS=1;
  const UInt_t nV=2;
  TString select[nS] = {"3jet"};
  TString var[nV]    = {"alphat","apcjetmetmin"};
  UInt_t  nBins[nV]  = {40,50};
  Float_t xFirst[nV] = {0,0};
  Float_t xLast[nV]  = {2,1};
  */
  // Produce 1 plot per {selection ; variable}
  for(UInt_t iS=0 ; iS<nS ; iS++) {

    if(verbose>1) cout << "- selection : " << select[iS] << endl;

    plot(select[iS], nV, var, nBins, xFirst, xLast, false, false, true, locProcesses, labelProc);

  }

  _outfile->Close();
  return 0;
}

Int_t XMetAnalysis::plot(TString select, const UInt_t nV, TString* var,
			 UInt_t* nBins, Float_t* xFirst, Float_t* xLast, 
			 Bool_t stack, Bool_t dolog, Bool_t unity,
			 vector<TString> locProcesses, vector<TString> labelProc)
{

  (*_outlog) << "Selection: " << select << endl;

  // Define selections //
  TCut weight="";
  TCut cut = defineCut(select, "signal");

  // Declare histograms //
  map<TString, map<TString,TH1F*> > mapVarHistos;
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax;

  const UInt_t nP = locProcesses.size();
  Float_t integral[nP][nV];

  // Loop over chains and generate histograms //
  //
  TString nameDir, locVar, variable;
  Int_t color;
  TH1F* hTemp;
  //
  for(UInt_t iP=0 ; iP<nP ; iP++) {

    // check if current requested process is available
    nameDir = locProcesses[iP];
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

    if(     nameDir.Contains("znn_cr")) {
      cut = defineCut(select, "zctrl");
    }
    else if(nameDir.Contains("wj_cr"))  {
      cut = defineCut(select, "wctrl");
    }

    // Skim the chain
    _mapProcess[nameDir].Skim(select, cut);

    // Loop over requested variables
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      if(verbose>1) cout << "--- variable: " << var[iV] << endl;

      // initialize integrals
      integral[iP][iV] = 0;
      
      // Define histogram and set style
      mapVarHistos[nameDir][var[iV]] = new TH1F("h_"+var[iV]+"_"+nameDir+"_"+select, 
						var[iV]+" "+nameDir+" "+select,
						nBins[iV], xFirst[iV], xLast[iV]);
      hTemp = mapVarHistos[nameDir][var[iV]];
      color   = _mapProcess[nameDir].GetColor();
      setStyle( hTemp , color );

      // Draw the variable
      locVar = var[iV];
      if(var[iV].Contains("phi")) locVar = "abs("+var[iV]+")";
      _mapProcess[nameDir].Draw(hTemp, locVar, cut, weight);
      cout << "--- drew variable: " << locVar << endl;

      // Normalize
      if( !nameDir.Contains("met") ) hTemp->Scale(_lumi*_rescale);
      integral[iP][iV] = hTemp->Integral();
      if(unity && integral[iP][iV]!=0) hTemp->Scale(1/integral[iP][iV]);
      cout << "--- normalized" << endl;

      // Determine extrema
      locMin = hTemp->GetMinimum();
      locMax = hTemp->GetMaximum();
      if(locMin<minPlot) minPlot = locMin;
      if(locMax>maxPlot) maxPlot = locMax;
    } // end loop:variables
  } // end loop:processes 
  cout << "-- end loop over processes" << endl;
  
  // Save histograms //
  _outfile->cd();
  cout << "-- _outfile->cd()... done" << endl;
  map<TString, map<TString,TH1F*> >::iterator itVarHistos;
  map<TString,TH1F*>::iterator itHistos;
  for( itVarHistos=mapVarHistos.begin() ; itVarHistos!=mapVarHistos.end() ; itVarHistos++) {
    for( itHistos=itVarHistos->second.begin() ; itHistos!=itVarHistos->second.end() ; itHistos++) {
      itHistos->second->Write();
      cout << "---- histo written: " 
	   << itHistos->second->GetName() 
	   << endl;
    }
  }

  // Produce the plot for each variable //

  /// Prepare TCanvas and TLegend
  TCanvas* c = new TCanvas("c","c",20,20,600,600);
  TLegend* leg = new TLegend(0.88,0.65,0.98,0.76,"","brNDC");
  setStyle(c);
  setStyle(leg);
  if(dolog) gPad->SetLogy();

  // Loop over variables
  for(UInt_t iV=0 ; iV<nV ; iV++) {

    // Yields outlog
    (*_outlog) << "Var: " << var[iV] << endl;
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      (*_outlog) << setw(10) << locProcesses[iP];
    }
    (*_outlog) << endl;

  
    // LOOP OVER PROCESSES' HISTO //
    Bool_t first=true;
    //
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      //
      // Get histogram
      hTemp = mapVarHistos[locProcesses[iP]][var[iV]];
      //      
      // Yields outlog
      (*_outlog) << setw(10) << integral[iP][iV];
      //
      if(hTemp) {
	if(first) {
	  first=false;
	  //
	  if(!stack) {
	    hTemp->SetMinimum(minPlot);
	    hTemp->SetMaximum(maxPlot);
	    if(!dolog) hTemp->SetMinimum(0.);
	    if(unity) hTemp->SetMaximum(1.1);
	  }
	  else {
	    
	  }
	  //
	  hTemp->Draw("HISTE1");
	}
	//	
	hTemp->Draw("HISTE1SAME");
      }
      //      
      if(iV==0) leg->AddEntry(hTemp,labelProc[iP],"L");
    }
    //
    (*_outlog) << endl;

    // Compute shape compatibility
  

    // Draw and Print
    leg->Draw();
    //
    if(nV>1) {
      if(iV==0)
	c->Print("plots/"+_tag+"/plots_"+select+".pdf(" , "Title:"+var[iV]);
      else if(iV==nV-1)
	c->Print("plots/"+_tag+"/plots_"+select+".pdf)" , "Title:"+var[iV]);
      else
	c->Print("plots/"+_tag+"/plots_"+select+".pdf"  , "Title:"+var[iV]);
    }
    else {
      c->Print("plots/"+_tag+"/plots_"+select+".pdf"  , "Title:"+var[iV]);
    }

  }
  
  (*_outlog) << endl;

  //////////////////////////

  // Delete histo pointers
  for( itVarHistos=mapVarHistos.begin() ; itVarHistos!=mapVarHistos.end() ; itVarHistos++) {
    for( itHistos=itVarHistos->second.begin() ; itHistos!=itVarHistos->second.end() ; itHistos++) {
      delete (itHistos->second);
      cout << "---- histo deleted" << endl;
    }
  }

  delete c;
  delete leg;

  // END
  return 1;
}

TCut XMetAnalysis::defineCut(TString select, TString region)
{

  // Trigger
  TCut trig  = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";

  // Lepton selection depending on the region
  TCut leptons = "";
  if(region=="zctrl") {
    leptons = "(nelectrons==0 && ntaus==0)";
    leptons *= "(zmass > 60 && zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1))";
  }
  else if(region=="wctrl") {
    leptons = "(nelectrons==0 && ntaus==0 && nmuons==1)";
    leptons *= "(wmt > 50 && wmt < 100 && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1)";
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
  /// Z
  _mapProcess["znn_cr_DA"] = XMetProcess("znn_cr_DA", kAzure+7,"reducedtree.root");
  _mapProcess["znn_cr_MC"] = XMetProcess("znn_cr_MC",kAzure+7,"reducedtree.root");
  /// W
  _mapProcess["wj_cr_DA" ] = XMetProcess("wj_cr_DA",  kGreen+2, "reducedtree.root");
  _mapProcess["wj_cr_MC" ] = XMetProcess("wj_cr_MC", kGreen+2, "reducedtree.root");

  // Signal
  //_mapProcess["dm_v_1"]   = XMetProcess("",   kOrange+3, "reducedtree.root");

  // Data
  _mapProcess["data"]   = XMetProcess("data",   kBlack, "reducedtree.root");

  // Sub-directories in _path
  _mapProcess["znn"].AddDir("znn");
  _mapProcess["zll"].AddDir("zll");
  _mapProcess["wjets"].AddDir("wjets");
  _mapProcess["ttbar"].AddDir("ttbar");
  _mapProcess["top"].AddDir("singletop");
  _mapProcess["qcd"].AddDir("qcd");
  _mapProcess["vv"].AddDir("dibosons");
  _mapProcess["znn_cr_DA"].AddDir("met");
  _mapProcess["znn_cr_MC"].AddDir("bkgz");
  _mapProcess["wj_cr_DA"].AddDir("met");
  _mapProcess["wj_cr_MC"].AddDir("bkgw");
  // signal still missing
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
