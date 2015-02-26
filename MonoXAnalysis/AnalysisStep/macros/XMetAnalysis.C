#include "XMetAnalysis.h"

using namespace std;
Int_t verbose = 2;

XMetAnalysis::XMetAnalysis(TString tag)
{
  _tag = tag;

  _path    = "/user/ndaci/Data/XMET/MonoJetTrees/V4/test/";
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
  const UInt_t nS=4;
  const UInt_t nV=5;
  TString select[nS] = {"alljets","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};
  UInt_t  nBins[nV]  = {40, 50, 50, 50,  50};
  Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  // Produce 1 plot per {selection ; variable}
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    if(verbose>1) cout << "- selection : " << select[iS] << endl;
    for(UInt_t iV=0 ; iV<nV ; iV++) {
      if(verbose>1) cout << "-- variable : " << var[iV] << endl;
      //
      plot(select[iS], var[iV], nBins[iV], xFirst[iV], xLast[iV], false, false, true, locProcesses, labelProc);
    }
  }
  _outfile->Close();

  return 0;
}

Int_t XMetAnalysis::plot(TString select, TString variable, 
			 Int_t nBins, Int_t xFirst, Int_t xLast, 
			 Bool_t stack, Bool_t dolog, Bool_t unity, 
			 vector<TString> locProcesses, vector<TString> labelProc)
{

  // Define selections //
  TCut weight;
  TCut cut = defineCut(select);

  // Declare histograms //
  map<TString, TH1F*> mapHistos;
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax, integral;

  // Loop over chains and generate histograms //
  TString nameDir, locVar;
  Int_t color;
  TChain* chain;
  // loop
  //for( _itProcess=_mapProcess.begin() ; _itProcess!=_mapProcess.end() ; _itProcess++) {
  for(UInt_t iP=0 ; iP<locProcesses.size() ; iP++) {

    // check if current requested process is available
    nameDir = locProcesses[iP];
    if(_mapProcess[nameDir]==_mapProcess.end()) { // FIXME ND
      if(verbose>1) cout << "-- ERROR: requested process '" << locProcesses[iP] << "' unavailable." << endl;
      continue;
    }

    // Process and corresponding chain
    chain   = mapProcess[nameDir].first;
    color   = mapProcess[nameDir].second.second;
    if(verbose>1) cout << "--- process : " << nameDir << endl;

    // Define histogram and set style
    mapHistos[nameDir] = new TH1F("h_"+variable+"_"+nameDir+"_"+select, 
				  variable+" "+nameDir+" "+select,
				  nBins, xFirst, xLast);
    setStyle( mapHistos[nameDir] , color );
    
    // Define reweighting
    if( nameDir.Contains("met") ) weight = "1";
    else weight = "puwgt*wgt";

    locVar = variable;
    if(variable.Contains("phi")) locVar = "abs("+variable+")";
    drawHistogram( chain, mapHistos[nameDir], locVar, cut*weight);

    // Normalize
    if( !nameDir.Contains("met") ) mapHistos[nameDir]->Scale(_lumi*_rescale);
    integral = mapHistos[nameDir]->Integral();
    if(unity && integral!=0) mapHistos[nameDir]->Scale(1/integral);

    // Determine extrema
    locMin = mapHistos[nameDir]->GetMinimum();
    locMax = mapHistos[nameDir]->GetMaximum();
    if(locMin<minPlot) minPlot = locMin;
    if(locMax>maxPlot) maxPlot = locMax;
  }


  // Save histograms //
  _outfile->cd();
  map<TString,TH1F*>::iterator itHistos;
  for( itHistos=mapHistos.begin() ; itHistos!=mapHistos.end() ; itHistos++) {
    itHistos->second->Write();
  }

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
    if(mapHistos[locProcesses[iP]]) {
      if(first) {
	first=false;

	if(!stack) {
	  mapHistos[locProcesses[iP]]->SetMinimum(minPlot);
	  mapHistos[locProcesses[iP]]->SetMaximum(maxPlot);
	  if(!dolog) mapHistos[locProcesses[iP]]->SetMinimum(0.);
	  if(unity) mapHistos[locProcesses[iP]]->SetMaximum(1.1);
	}
	else {

	}

	mapHistos[locProcesses[iP]]->Draw("HISTE1");
      }

      mapHistos[locProcesses[iP]]->Draw("HISTE1SAME");
    }

    leg->AddEntry(locProcesses[iP],labelProc[iP],"P");
  }
  
  // Compute shape compatibility
  

  // Draw and Print
  leg->Draw();
  //
  c.Print("plots/"+_tag+"/plot_"+select+"_"+variable+".png","png");
  c.Print("plots/"+_tag+"/plot_"+select+"_"+variable+".pdf","pdf");
  c.Print("plots/"+_tag+"/plot_"+select+"_"+variable+".eps","eps");
  //////////////////////////

  return 1;
}

TCut XMetAnalysis::defineCut(TString select)
{

  TCut trig   = "(hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
  TCut veto   = "(nmuons == 0 && nelectrons == 0 && ntaus == 0)";

  TCut metID  = "(abs(pfmet - calomet) < 2*calomet)" ;

  TCut jetID1 = "(signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2)";
  TCut jetID2 = "(secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9 && secondjetpt>30 && abs(secondjeteta)<2.5)";
  TCut jetID3 = "thirdjetpt>30 && abs(thirdjeteta)<2.5";
  //TCut jetID3 = "(thirdjetNHfrac  < 0.7 && thirdjetEMfrac  < 0.9)";
  //TCut jetIDMult = "(njets==1 || ( (secondjetNHfrac<0.7 && secondjetEMfrac<0.9)&&(njets==2 || (njets==3 && thirdjetNHfrac<0.7 && thirdjetEMfrac<0.9) ) ) )" ;
  TCut jetIDMult = "(njets==1 || ( (secondjetNHfrac<0.7 && secondjetEMfrac<0.9 && secondjetpt>30 && abs(secondjeteta)<2.5)&&(njets==2 || (njets==3 && thirdjetpt>30 && abs(thirdjeteta)<2.5) ) ) )" ;

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

int XMetAnalysis::drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut) 
{
  tree->Draw(variable+">>"+TString(h->GetName()), cut, "", 1000);
  //tree->Draw(variable+">>"+TString(h->GetName()), cut);
  return 0;
}

int XMetAnalysis::setStyle(TH1F* h, Int_t color)
{
  h->Sumw2();

  //h->SetMarkerSize(0.5);
  //h->SetMarkerStyle(kPlus);

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
  _mapProcess["singletop"].second.second = kOrange;
  _mapProcess["qcd"].second.second       = kYellow;
  _mapProcess["dibosons"].second.second  = kRed;
  _mapProcess["zll"].second.second       = kMagenta+8;
  //
  _mapProcess["znn"].second.first.push_back("znn100to200");
  _mapProcess["znn"].second.first.push_back("znn200to400");
  _mapProcess["znn"].second.first.push_back("znn400toinf");
  _mapProcess["znn"].second.first.push_back("znn50to100");
  //
  _mapProcess["wln"].second.first.push_back("wln");
  //_mapProcess["wln"].second.first.push_back("w4jets");
  // 
  _mapProcess["ttbar"].second.first.push_back("ttbar"); 
  //
  _mapProcess["singletop"].second.first.push_back("singletbars");
  _mapProcess["singletop"].second.first.push_back("singletbart");
  _mapProcess["singletop"].second.first.push_back("singletbarw");
  _mapProcess["singletop"].second.first.push_back("singlets");
  _mapProcess["singletop"].second.first.push_back("singlett");
  _mapProcess["singletop"].second.first.push_back("singletw");
  //
  _mapProcess["qcd"].second.first.push_back("qcd1000to1400");
  _mapProcess["qcd"].second.first.push_back("qcd120to170");
  _mapProcess["qcd"].second.first.push_back("qcd1400to1800");
  _mapProcess["qcd"].second.first.push_back("qcd170to300");
  _mapProcess["qcd"].second.first.push_back("qcd1800toinf");
  _mapProcess["qcd"].second.first.push_back("qcd300to470");
  _mapProcess["qcd"].second.first.push_back("qcd470to600");
  _mapProcess["qcd"].second.first.push_back("qcd600to800");
  _mapProcess["qcd"].second.first.push_back("qcd800to1000");
  _mapProcess["qcd"].second.first.push_back("qcd80to120");
  //
  _mapProcess["dibosons"].second.first.push_back("ww"); 
  _mapProcess["dibosons"].second.first.push_back("zz"); 
  _mapProcess["dibosons"].second.first.push_back("wz"); 
  //
  _mapProcess["zll"].second.first.push_back("zll");
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
