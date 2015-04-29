#include "XMetAnalysis.h"

using namespace std;
Int_t verbose = 2;

XMetAnalysis::XMetAnalysis(TString tag)
{
  _tag = tag;

  //_path    = "/user/ndaci/Data/XMET/MonoJetTrees/V4/test/";
  _path    = "/user/ndaci/Data/XMET/Phys14/MonoJet7XTrees/";
  _lumi    = 19.7;
  _rescale = 1.0;
  _outfile = new TFile("plots/"+_tag+"/plots_"+_tag+".root","recreate");

  DefineChains();
}

XMetAnalysis::~XMetAnalysis()
{
  // End job //
  //delete _outfile;
  /*
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    delete trees[iP];
    delete mapHistos[nameDir];
    files[iP]->Close();
    //delete files[iP];
    //delete trees[iP];
    //delete mapHistos[nameDir];
  }
  */
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
  //UInt_t  nBins[nV]  = {800, 1000, 1000, 1000,  1000};
  UInt_t  nBins[nV]  = {4000, 5000, 5000, 5000,  5000};
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

  // Define selections //
  TCut weight;
  TCut cut = defineCut(select);

  // Declare histograms //
  map<TString, map<TString,TH1F*> > mapVarHistos;
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax, integral;

  // Loop over chains and generate histograms //
  //
  TString nameDir, locVar, variable;
  Int_t color;
  TChain* chain;
  TH1F* hTemp;
  //
  for(UInt_t iP=0 ; iP<locProcesses.size() ; iP++) {

    // check if current requested process is available
    nameDir = locProcesses[iP];
    if(_mapProcess.find(nameDir)==_mapProcess.end()) { 
      if(verbose>1) cout << "-- ERROR: requested process '" << locProcesses[iP] << "' unavailable." << endl;
      continue;
    }

    // Process and corresponding chain
    chain   = _mapProcess[nameDir].first;
    color   = _mapProcess[nameDir].second.second;
    if(verbose>1) cout << "--- process : " << nameDir << endl;

    // Define reweighting
    if( nameDir.Contains("met") ) weight = "1";
    else weight = "puwgt*wgt";

    // Define skim for current process and requested selection
    //chain->Draw(">>+skim_"+nameDir, cut,"entrylist",100);
    chain->Draw(">>+skim_"+nameDir, cut,"entrylist");
    TEntryList *skim = (TEntryList*)gDirectory->Get("skim_"+nameDir);
    int nEntries = skim->GetN();
    chain->SetEntryList(skim);
    if(verbose>1) cout << "--- produced skim : " << nEntries << " entries" << endl;

    // Loop over requested variables
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      // Define histogram and set style
      mapVarHistos[nameDir][var[iV]] = new TH1F("h_"+var[iV]+"_"+nameDir+"_"+select, 
						var[iV]+" "+nameDir+" "+select,
						nBins[iV], xFirst[iV], xLast[iV]);
      hTemp = mapVarHistos[nameDir][var[iV]];
      setStyle( hTemp , color );

      // Draw the variable
      locVar = var[iV];
      if(var[iV].Contains("phi")) locVar = "abs("+var[iV]+")";
      chain->Draw(locVar+">>"+TString(hTemp->GetName()), cut*weight);
      //chain->Draw(locVar+">>"+TString(hTemp->GetName()), cut*weight, "", 100); // FIXME ND

      // Normalize
      if( !nameDir.Contains("met") ) hTemp->Scale(_lumi*_rescale);
      integral = hTemp->Integral();
      if(unity && integral!=0) hTemp->Scale(1/integral);

      // Determine extrema
      locMin = hTemp->GetMinimum();
      locMax = hTemp->GetMaximum();
      if(locMin<minPlot) minPlot = locMin;
      if(locMax>maxPlot) maxPlot = locMax;
    }

    chain->SetEntryList(0);
    skim->Reset();
  }
  
  // Save histograms //
  _outfile->cd();
  map<TString, map<TString,TH1F*> >::iterator itVarHistos;
  map<TString,TH1F*>::iterator itHistos;
  for( itVarHistos=mapVarHistos.begin() ; itVarHistos!=mapVarHistos.end() ; itVarHistos++) {
    for( itHistos=itVarHistos->second.begin() ; itHistos!=itVarHistos->second.end() ; itHistos++) {
      itHistos->second->Write();
    }
  }

  // Produce the plot for each variable //
  for(UInt_t iV=0 ; iV<nV ; iV++) {

    // Prepare TCanvas /////
    TCanvas c("c","c",20,20,600,600);
    c.SetFillColor(0);
    c.SetBorderMode(0);
    c.SetBorderSize(2);
    c.SetFrameBorderMode(0);
    c.SetFrameBorderMode(0);
    //
    if(dolog) gPad->SetLogy();
    //
    // Add legend
    TLegend* leg = new TLegend(0.88,0.65,0.98,0.76,"","brNDC");
    leg->SetLineColor(1);
    leg->SetTextColor(1);
    leg->SetTextFont(42);
    leg->SetTextSize(0.0244755);
    leg->SetShadowColor(kWhite);
    leg->SetFillColor(kWhite);  
  
    // LOOP OVER PROCESSES' HISTO //
    Bool_t first=true;
    
    for(UInt_t iP=0 ; iP<locProcesses.size() ; iP++) {

      hTemp = mapVarHistos[locProcesses[iP]][var[iV]];

      if(hTemp) {
	if(first) {
	  first=false;
	  
	  if(!stack) {
	    hTemp->SetMinimum(minPlot);
	    hTemp->SetMaximum(maxPlot);
	    if(!dolog) hTemp->SetMinimum(0.);
	    if(unity) hTemp->SetMaximum(1.1);
	  }
	  else {
	    
	  }

	  hTemp->Draw("HISTE1");
	}
	
	hTemp->Draw("HISTE1SAME");
      }
      
      leg->AddEntry(hTemp,labelProc[iP],"L");
    }
  
    // Compute shape compatibility
  

    // Draw and Print
    leg->Draw();
    //
    if(iV==0)
      c.Print("plots/"+_tag+"/plots_"+select+".pdf(" , "Title:"+var[iV]);
    else if(iV==nV-1)
      c.Print("plots/"+_tag+"/plots_"+select+".pdf)" , "Title:"+var[iV]);
    else
      c.Print("plots/"+_tag+"/plots_"+select+".pdf"  , "Title:"+var[iV]);

  }
  //////////////////////////

  return 1;
}

TCut XMetAnalysis::defineCut(TString select)
{

  TCut trig   = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
  TCut veto   = "(nmuons == 0 && nelectrons == 0 && ntaus == 0)";

  TCut metID  = "(abs(pfmet - calomet) < 2*calomet)" ;
  TCut metCut = "" ;

  TCut jetID1 = "(signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2)";
  TCut jetID2 = "(secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9 && secondjetpt>30 && abs(secondjeteta)<2.5)";
  TCut jetID3 = "thirdjetpt>30 && abs(thirdjeteta)<4.5"; // FIXME ND
  //TCut jetID3 = "(thirdjetNHfrac  < 0.7 && thirdjetEMfrac  < 0.9)";
  //TCut jetIDMult = "(njets==1 || ( (secondjetNHfrac<0.7 && secondjetEMfrac<0.9)&&(njets==2 || (njets==3 && thirdjetNHfrac<0.7 && thirdjetEMfrac<0.9) ) ) )" ;
  TCut jetIDMult = "(njets==1 || ( (secondjetNHfrac<0.7 && secondjetEMfrac<0.9 && secondjetpt>30 && abs(secondjeteta)<2.5)&&(njets==2 || (njets==3 && thirdjetpt>30 && abs(thirdjeteta)<2.5) ) ) )" ;
  TCut jetIDMono = "(njets==1 || (njets==2 && secondjetNHfrac<0.7 && secondjetEMfrac<0.9 && secondjetpt>30 && abs(secondjeteta)<2.5) )" ;

  TCut jetKine1 = "signaljetpt > 110 && abs(signaljeteta) < 2.4";

  TCut jetID;
  TCut jetBin;
  TCut dphi;
  TCut alphat;
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
  if(     select.Contains("dphi"))      noqcd = dphi;
  else if(select.Contains("alphat"))    noqcd = alphat;
  else if(select.Contains("apcjetmetmax")) noqcd = apcjetmetmax;

  //cout << trig*veto*metID*noqcd*jetID*jetKine1*jetBin << endl;

  return (trig*veto*metID*noqcd*jetID*jetKine1*jetBin);
}

int XMetAnalysis::setStyle(TH1F* h, Int_t color)
{
  h->Sumw2();

  h->SetMarkerSize(0.5);
  h->SetMarkerStyle(kPlus);

  h->SetLineColor(color);
  h->SetMarkerColor(color);
  //h->SetFillColor(color);

  return 0;
}

Int_t XMetAnalysis::DefineChains()
{
  if(verbose>1) cout << "- begin DefineChains()" << endl;

  // Processes
  //
  // {"bkgnowz","bkgw","bkgz","dibosons","met","qcd","singletop","ttbar","wjets","zjets","znn"};
  // {kGreen, kGreen, kBlue, kRed, kBlack, kYellow, kOrange, kViolet, kGreen, kMagenta-10, kBlue};
  //
  _mapProcess["znn"].second.second       = kBlue;
  _mapProcess["wln"].second.second       = kGreen;
  _mapProcess["ttbar"].second.second     = kViolet;
  _mapProcess["singletop"].second.second = kYellow;
  _mapProcess["qcd"].second.second       = kOrange+2;
  _mapProcess["dibosons"].second.second  = kRed;
  _mapProcess["zll"].second.second       = kMagenta+8;
  //
  _mapProcess["znn"].second.first.push_back("znunu");
  //
  _mapProcess["wln"].second.first.push_back("wjets");
  // 
  _mapProcess["ttbar"].second.first.push_back("ttbar"); 
  //
  //_mapProcess["singletop"].second.first.push_back("");
  //
  _mapProcess["qcd"].second.first.push_back("qcd");
  //
  //_mapProcess["dibosons"].second.first.push_back("ww"); 
  //_mapProcess["dibosons"].second.first.push_back("zz"); 
  //_mapProcess["dibosons"].second.first.push_back("wz"); 
  //
  _mapProcess["zll"].second.first.push_back("zjets");
  //

  if(verbose>1) cout << "#entries in mapProcess : " << _mapProcess.size() << endl;
  
  //
  // Add the files to the chains
  //
  for( _itProcess=_mapProcess.begin() ; _itProcess!=_mapProcess.end() ; _itProcess++ ) {

    _itProcess->second.first = new TChain("tree/tree");
    vector<TString> theDirs = _itProcess->second.second.first;

    for(UInt_t iD=0 ; iD<theDirs.size() ; iD++) {
      _itProcess->second.first->Add(_path+"/"+theDirs[iD]+"/tree_*.root");
    }
    //if(verbose>1) cout << "number of entries : " << _itProcess->second.first->GetEntries() << endl;
  }

  if(verbose>1) cout << "- end DefineChains()" << endl;

  return 0;
}
