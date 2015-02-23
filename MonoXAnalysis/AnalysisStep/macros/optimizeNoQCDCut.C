#include <iostream>
#include <sstream>
#include <map>
//
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
//
#include "tdrstyle.h"

using namespace std;
int verbose = 2;
bool dolog = false;

// Declare functions ////////
int  setStyle(TH1F* h, Int_t color);
int  drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut);
TCut defineCut(TString select);
int  plot(TString tag, TFile* outfile, double lumi, double rescale, bool unity,
	  map<TString,TChain*> mapChains, map<TString,Int_t> mapColors, 
	  TString select, TString variable, int nBins, int xFirst, int xLast);
int optimizeNoQCDCut(TString tag, bool unity);
/////////////////////////////


int optimizeNoQCDCut(TString tag, bool unity)
{

  // STYLE
  //gROOT->Reset();
  //loadPresentationStyle();
  //gROOT->ForceStyle();

  // Input/output and conditions
  TFile* outfile = new TFile("plots/"+tag+"/plots_"+tag+".root","recreate");
  TString path   ="/user/ndaci/Data/XMET/MonoJetTrees/V4/test/";
  Double_t lumi  = 19.7;

  // Selections and variables
  const UInt_t nS=4;
  const UInt_t nV=5;
  TString select[nS] = {"alljets","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};
  UInt_t  nBins[nV]  = {40, 50, 50, 50,  50};
  Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  // Processes
  map<TString, TChain*> mapChains;
  map<TString, Int_t>   mapColors;
  //
  const UInt_t nP=2;
  TString nameDir[nP]={"qcd","znn"};
  Int_t   color[nP] = {kRed, kBlue};
  //
  vector<TString> dir[nP];
  //
  dir[0].push_back("qcd1000to1400");
  dir[0].push_back("qcd120to170");
  dir[0].push_back("qcd1400to1800");
  dir[0].push_back("qcd170to300");
  dir[0].push_back("qcd1800toinf");
  dir[0].push_back("qcd300to470");
  dir[0].push_back("qcd470to600");
  dir[0].push_back("qcd600to800");
  dir[0].push_back("qcd800to1000");
  dir[0].push_back("qcd80to120");
  //
  dir[1].push_back("znn100to200");
  dir[1].push_back("znn200to400");
  dir[1].push_back("znn400toinf");
  dir[1].push_back("znn50to100");
  //
  // Add the files to the chains
  //
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    mapColors[nameDir[iP]] = color[iP];
    mapChains[nameDir[iP]] = new TChain("tree/tree");
    for(UInt_t iD=0 ; iD<dir[iP].size() ; iD++) {
      mapChains[nameDir[iP]]->Add(path+"/"+dir[iP][iD]+"/tree_*.root");
    }
    //cout << "number of entries : " << trees[iP]->GetEntries() << endl;
  }

  // Produce 1 plot per {selection ; variable}
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    if(verbose>1) cout << "- selection : " << select[iS] << endl;
    for(UInt_t iV=0 ; iV<nV ; iV++) {
      if(verbose>1) cout << "-- variable : " << var[iV] << endl;
      //
      plot(tag, outfile, lumi, 1.0, unity, mapChains, mapColors, select[iS], var[iV], nBins[iV], xFirst[iV], xLast[iV]);
    }
  }
  outfile->Close();

  return 0;
}

int plot(TString tag, TFile* outfile, double lumi, double rescale, bool unity,
	 map<TString,TChain*> mapChains, map<TString,Int_t> mapColors,
	 TString select, TString variable, int nBins, int xFirst, int xLast)
{

  // Define selections //
  TCut weight;
  TCut cut = defineCut(select);

  // Declare histograms
  map<TString, TH1F*> mapHistos;
  Float_t minPlot = 9999999.;
  Float_t maxPlot = -9999999.;
  Float_t locMin, locMax, integral;

  // Loop over chains and generate histograms
  map<TString,TChain*>::iterator itChain;
  TString nameDir, locVar;
  TChain* chain;
  //
  for( itChain=mapChains.begin() ; itChain!=mapChains.end() ; itChain++) {

    // Process and corresponding chain
    nameDir = itChain->first;
    chain   = itChain->second;
    if(verbose>1) cout << "--- process : " << nameDir << endl;

    // Define histogram and set style
    mapHistos[nameDir] = new TH1F("h_"+variable+"_"+nameDir+"_"+select, 
				  variable+" "+nameDir+" "+select,
				  nBins, xFirst, xLast);
    setStyle( mapHistos[nameDir] , mapColors[nameDir] );

    // Define reweighting
    if( nameDir.Contains("met") ) weight = "1";
    else weight = "puwgt*wgt";

    locVar = variable;
    if(variable.Contains("phi")) locVar = "abs("+variable+")";
    drawHistogram( chain, mapHistos[nameDir], locVar, cut*weight);

    // Normalize
    if( !nameDir.Contains("met") ) mapHistos[nameDir]->Scale(lumi*rescale);
    integral = mapHistos[nameDir]->Integral();
    if(unity && integral!=0) mapHistos[nameDir]->Scale(1/integral);

    // Determine extrema
    locMin = mapHistos[nameDir]->GetMinimum();
    locMax = mapHistos[nameDir]->GetMaximum();
    if(locMin<minPlot) minPlot = locMin;
    if(locMax>maxPlot) maxPlot = locMax;
  }

  // Save histograms //
  outfile->cd();
  map<TString,TH1F*>::iterator itHistos;
  for( itHistos=mapHistos.begin() ; itHistos!=mapHistos.end() ; itHistos++) {
    itHistos->second->Write();
  }

  // Produce plot /////
  TCanvas c("c","c",20,20,600,600);
  c.SetFillColor(0);
  c.SetBorderMode(0);
  c.SetBorderSize(2);
  c.SetFrameBorderMode(0);
  c.SetFrameBorderMode(0);
  //
  if(dolog) gPad->SetLogy();
  //
  if(mapHistos["znn"]) {
    mapHistos["znn"]->SetMinimum(minPlot);
    mapHistos["znn"]->SetMaximum(maxPlot);
    if(unity) {
      mapHistos["znn"]->SetMaximum(1.1);
      if(!dolog) mapHistos["znn"]->SetMinimum(0.);
    }
    mapHistos["znn"]->Draw("HISTE1");
  }
  if(mapHistos["qcd"]) mapHistos["qcd"]->Draw("HISTE1SAME");
  //
  // Add legend
  //TLegend* leg = new TLegend(0.78,0.65,0.98,0.76,"","brNDC");
  TLegend* leg = new TLegend(0.88,0.65,0.98,0.76,"","brNDC");
  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0244755);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);  
  //
  leg->AddEntry(mapHistos["znn"],"Z(#nu#nu)");
  leg->AddEntry(mapHistos["qcd"],"QCD");
  leg->Draw();
  //
  c.Print("plots/"+tag+"/plot_"+select+"_"+variable+".png","png");
  c.Print("plots/"+tag+"/plot_"+select+"_"+variable+".pdf","pdf");
  c.Print("plots/"+tag+"/plot_"+select+"_"+variable+".eps","eps");
  //////////////////////////


  // End job //
  //delete outfile;
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

  return 1;
}

TCut defineCut(TString select)
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

int drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut) 
{
  //tree->Draw(variable+">>"+TString(h->GetName()), cut, "", 1000);
  tree->Draw(variable+">>"+TString(h->GetName()), cut);
  return 0;
}

int setStyle(TH1F* h, Int_t color)
{
  h->Sumw2();

  //h->SetMarkerSize(0.5);
  //h->SetMarkerStyle(kPlus);

  h->SetLineColor(color);
  //h->SetMarkerColor(color);
  //h->SetFillColor(color);

  return 0;
}
