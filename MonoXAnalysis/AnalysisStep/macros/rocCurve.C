#include <iostream>
#include <sstream>
#include <map>
#include <utility>
//
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEntryList.h"
//
#include "tdrstyle.h"

Int_t rocCurve(TString _tag="v11_XMA_QCD_ROC_FinerBinning", bool dolog=false)
{

  // Input file
  TFile* file = new TFile("plots/"+_tag+"/plots_"+_tag+".root","read");
  
  // Output canvas
  TCanvas cRoc("cRoc","cRoc",20,20,600,600);
  cRoc.SetFillColor(0);
  cRoc.SetBorderMode(0);
  cRoc.SetBorderSize(2);
  cRoc.SetFrameBorderMode(0);
  cRoc.SetFrameBorderMode(0);
  //
  if(dolog) gPad->SetLogy();

  // Selections and variables
  const UInt_t nWd=2; // forward/backward cumulative distribution
  const UInt_t nS=5;
  const UInt_t nV=5;

  TGraph* gRoc[nS][nV][nWd];

  TString wd[nWd]    = {"upcut","lowcut"};
  TString wdT[nWd]   = {"(upper cut)","(lower cut)"};
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};

  Int_t colors[nV]   = {kBlack, kBlue, kRed, kGreen+2, kMagenta};

  //UInt_t  nBins[nV]  = {40, 50, 50, 50,  50};
  //UInt_t  nBins[nV]  = {400, 500, 500, 500, 500};
  UInt_t  nBins[nV]  = {800, 1000, 1000, 1000,  1000};
  //UInt_t  nBins[nV]  = {4000, 5000, 5000, 5000,  5000};
  //UInt_t  nBins[nV]  = {8000, 10000, 10000, 10000,  10000};
  
  //Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  //Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      TH1F* h_qcd = (TH1F*) file->Get("h_"+var[iV]+"_qcd_"+select[iS]);
      TH1F* h_znn = (TH1F*) file->Get("h_"+var[iV]+"_znn_"+select[iS]);

      if(!h_qcd) {
	cout << "ERROR: missing histogram" << endl
	     << "h_"+var[iV]+"_qcd_"+select[iS] << endl;
      }
      if(!h_znn) {
	cout << "ERROR: missing histogram" << endl
	     << "h_"+var[iV]+"_znn_"+select[iS] << endl;
      }
      if(!h_qcd || !h_znn) return -1;

      cout << "-- "+var[iV]+" "+select[iS]+" : Integrals : "
	   << "qcd=" << h_qcd->Integral() << " "
	   << "znn=" << h_znn->Integral() << " " << endl;

      const UInt_t nB = nBins[iV];

      for(UInt_t iWd=0 ; iWd<nWd ; iWd++) {
	Double_t cum_qcd=0, cum_znn=0;
	Double_t yCum_qcd[nB];
	Double_t yCum_znn[nB];

	for(UInt_t iB=0 ; iB<nB ; iB++) {

	  cum_qcd += h_qcd->GetBinContent(iB+1);
	  cum_znn += h_znn->GetBinContent(iB+1);
	  //cum_qcd += h_qcd->GetBinContent(nB-iB);
	  //cum_znn += h_znn->GetBinContent(nB-iB);

	  if(iWd==0) { // upper cut : x < cut
	    yCum_qcd[iB] = cum_qcd;
	    yCum_znn[iB] = cum_znn;
	  }
	  else if(iWd==1) { // lower cut : x > cut
	    yCum_qcd[iB] = 1 - cum_qcd;
	    yCum_znn[iB] = 1 - cum_znn;
	  }
	  
	}

	gRoc[iS][iV][iWd] = new TGraph(nB, yCum_znn, yCum_qcd);
	gRoc[iS][iV][iWd]->SetMarkerStyle(kOpenSquare);
	gRoc[iS][iV][iWd]->SetMarkerColor(colors[iV]);
	gRoc[iS][iV][iWd]->SetLineColor(colors[iV]);
	//gRoc[iS][iV][iWd]->SetMarkerSize();
	gRoc[iS][iV][iWd]->SetFillColor(kWhite);
	gRoc[iS][iV][iWd]->GetXaxis()->Set(nB, 0, 1);
	gRoc[iS][iV][iWd]->GetYaxis()->Set(nB, 0, 1);
	gRoc[iS][iV][iWd]->SetMinimum(0);
	gRoc[iS][iV][iWd]->SetMaximum(1.05);
	gRoc[iS][iV][iWd]->GetXaxis()->SetTitle("Z(#nu#nu) efficiency");
	gRoc[iS][iV][iWd]->GetYaxis()->SetTitle("QCD efficiency");
	gRoc[iS][iV][iWd]->SetTitle("ROC Curve : "+select[iS]+" "+var[iV]+" "+wdT[iWd]);
	gRoc[iS][iV][iWd]->Draw("AL");

	//cRoc.Print("plots/"+_tag+"/roc_"+select[iS]+"_"+var[iV]+"_"+wd[iWd]+".pdf","pdf");

	if(iS==0 && iV==0 && iWd==0) 
	  cRoc.Print("plots/"+_tag+"/roc.pdf(","Title:"+select[iS]+" "+var[iV]+" "+wd[iWd]);
	else if(iS==nS-1 && iV==nV-1 && iWd==nWd-1)
	  cRoc.Print("plots/"+_tag+"/roc.pdf)","Title:"+select[iS]+" "+var[iV]+" "+wd[iWd]);
	else
	  cRoc.Print("plots/"+_tag+"/roc.pdf","Title:"+select[iS]+" "+var[iV]+" "+wd[iWd]);

      } // end loop over nWd
    }   // end loop over nV
  }     // end loop over nS

  // Put several killers per plot
  for(UInt_t iS=0 ; iS<nS ; iS++) {

    gRoc[iS][0][1]->SetTitle("Selection: "+select[iS]);

    gRoc[iS][0][1]->Draw("AL"); // alphat
    gRoc[iS][1][1]->Draw("L");  // apcjetmetmax
    gRoc[iS][2][1]->Draw("L");  // apcjetmetmin
    gRoc[iS][3][0]->Draw("L");  // jetjetdphi
    gRoc[iS][4][1]->Draw("L");  // jetmetdphimin

    TLegend* leg = new TLegend(0.38,0.6,0.58,0.8,"","brNDC");
    leg->SetLineColor(1);
    leg->SetTextColor(1);
    leg->SetTextFont(42);
    leg->SetTextSize(0.0244755);
    leg->SetShadowColor(kWhite);
    leg->SetFillColor(kWhite);

    for(UInt_t iV=0 ; iV<nV ; iV++)
      leg->AddEntry(gRoc[iS][iV][1], var[iV], "L");

    leg->Draw();

    TString metCut;
    if(_tag.Contains("NoMetCut"))    metCut = "NoMetCut";
    else if(_tag.Contains("Met200")) metCut = "Met200";
    else if(_tag.Contains("Met350")) metCut = "Met350";
    else if(_tag.Contains("MetFrom0to200"))   metCut = "MetFrom0to200";
    else if(_tag.Contains("MetFrom200to250")) metCut = "MetFrom200to250";
    else if(_tag.Contains("MetFrom250to350")) metCut = "MetFrom250to350";
    TString mysel = select[iS]+"_"+metCut;

    if(iS==0) 
      cRoc.Print("plots/"+_tag+"/allroc.pdf(","Title:"+mysel);
    else if(iS==nS-1)
      cRoc.Print("plots/"+_tag+"/allroc.pdf)","Title:"+mysel);
    else
      cRoc.Print("plots/"+_tag+"/allroc.pdf","Title:"+mysel);

  }

  return 0;
}
