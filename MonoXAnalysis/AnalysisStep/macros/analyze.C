#include <iostream>
#include <sstream>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
//#include ".h"
//#include ".h"

using namespace std;

// Declare functions ////////
int  setStyle(TH1F* h, Int_t color);
int  drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut);
TCut defineCut(TString select);
int  plot(TFile* outfile, double lumi, double rescale, TString dirTrees, 
	  TString select, TString variable, int nBins, int xFirst, int xLast);
int analyze();
/////////////////////////////


int analyze()
{

  TFile* outfile = new TFile("plots.root","recreate");
  /*
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "alljets_dphi", "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "1jet_dphi",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "2jet_dphi",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "3jet_dphi",    "mumet", 20, 0, 1000);
  */
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "alljets_alphat", "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "1jet_alphat",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "2jet_alphat",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "3jet_alphat",    "mumet", 20, 0, 1000);

  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "alljets_apcjetmax", "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "1jet_apcjetmax",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "2jet_apcjetmax",    "mumet", 20, 0, 1000);
  plot(outfile, 19.7, 1.0, "/user/ndaci/Data/XMET/MonoJetTrees/V2/", "3jet_apcjetmax",    "mumet", 20, 0, 1000);

  outfile->Close();

  return 0;
}

int plot(TFile* outfile, double lumi, double rescale,TString dirTrees, 
	 TString select, TString variable, int nBins, int xFirst, int xLast)
{

  // Define processes //
  const UInt_t nP=11;
  TString dir[nP] = {"bkgnowz","bkgw","bkgz","dibosons","met","qcd","singletop","ttbar","wjets","zjets","znunu"};
  Int_t color[nP] = {kGreen, kGreen, kBlue, kRed, kBlack, kYellow, kOrange, kViolet, kGreen, kMagenta-10, kBlue};
  map<TString , Int_t> mapNameHisto;

  // Define selections //
  TCut weight;
  TCut cut = defineCut(select);

  TFile* files[nP];
  TTree* trees[nP];
  TH1F*  histo[nP];

  TString treename = "tree/tree";

  for(UInt_t iP=0 ; iP<nP ; iP++) {

    mapNameHisto[dir[iP]] = iP;

    files[iP] = new TFile(dirTrees+"/"+dir[iP]+"/tree.root");
    trees[iP] = (TTree*)files[iP]->Get(treename);
    histo[iP] = new TH1F("h_"+variable+"_"+dir[iP]+"_"+select, 
			 variable+" "+dir[iP]+" "+select,
			 nBins, xFirst, xLast);

    setStyle( histo[iP] , color[iP] );

    if( dir[iP].Contains("met") ) weight = "1";
    else weight = "puwgt*wgt";
    drawHistogram(trees[iP], histo[iP], variable, cut*weight);
    if( !dir[iP].Contains("met") ) histo[iP]->Scale(lumi*rescale);
  }

  // Save histograms //
  outfile->cd();
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    histo[iP]->Write();
  }

  // Produce stack plot /////
  TCanvas c("c","c",20,20,600,600);
  gPad->SetLogy();

  THStack* hStack = new THStack("hStack_"+select,"");
  hStack->Add(histo[mapNameHisto["zjets"]]);
  hStack->Add(histo[mapNameHisto["dibosons"]]);
  hStack->Add(histo[mapNameHisto["qcd"]]);
  hStack->Add(histo[mapNameHisto["singletop"]]);
  hStack->Add(histo[mapNameHisto["ttbar"]]);
  hStack->Add(histo[mapNameHisto["wjets"]]);
  hStack->Add(histo[mapNameHisto["znunu"]]);

  histo[mapNameHisto["met"]]->Draw("P");
  hStack->Draw("HISTSAME");
  histo[mapNameHisto["met"]]->Draw("PSAME");

  c.Print("plots/stack_"+select+".pdf");
  //////////////////////////


  // End job //
  //delete outfile;
  /*
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    delete trees[iP];
    delete histo[iP];
    files[iP]->Close();
    //delete files[iP];
    //delete trees[iP];
    //delete histo[iP];
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
  TCut alphat="alphat>0.55";
  TCut apcjetmax="apcjetmax>0.8";

  if(     select.Contains("alljets")) {
    jetBin = "njets>=1 && njets<=3";
    jetID  = jetID1*jetIDMult;
    dphi   = "njets==1 || abs(jetjetdphi) < 2.5";
  }
  else if(select.Contains("1jet")) { 
    jetBin = "njets==1"; 
    jetID = jetID1; 
    dphi = "1";
  }
  else if(select.Contains("2jet")) { 
    jetBin = "njets==2"; 
    jetID = jetID1*jetID2; 
    dphi = "abs(jetjetdphi) < 2.5";
  }
  else if(select.Contains("3jet")) {
    jetBin = "njets==3";
    jetID = jetID1*jetID2*jetID3;
    dphi = "abs(jetjetdphi) < 2.5";
  }

  TCut noqcd="1";
  if(     select.Contains("dphi"))      noqcd = dphi;
  else if(select.Contains("alphat"))    noqcd = alphat;
  else if(select.Contains("apcjetmax")) noqcd = apcjetmax;

  cout << trig*veto*metID*noqcd*jetID*jetKine1*jetBin << endl;

  return (trig*veto*metID*noqcd*jetID*jetKine1*jetBin);
}

int drawHistogram(TTree* tree, TH1F* h, TString variable, TCut cut) 
{
  tree->Draw(variable+">>"+TString(h->GetName()), cut);
  return 0;
}

int setStyle(TH1F* h, Int_t color)
{
  h->Sumw2();

  h->SetMarkerSize(1.2);
  h->SetMarkerStyle(20);

  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetFillColor(color);

  return 0;
}
