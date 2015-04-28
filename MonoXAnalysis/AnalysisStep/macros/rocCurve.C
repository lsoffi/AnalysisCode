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

Int_t rocCurve(TString _tag="v10_XMA_QCD_ROC", bool dolog=false)
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

  TString wd[nWd]    = {"upcut","lowcut"};
  TString wdT[nWd]   = {"(upper cut)","(lower cut)"};
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};
  UInt_t  nBins[nV]  = {400, 500, 500, 500, 500};
  //Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  //Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      TH1F* h_qcd = (TH1F*) file->Get("h_"+var[iV]+"_qcd_"+select[iS]);
      TH1F* h_znn = (TH1F*) file->Get("h_"+var[iV]+"_znn_"+select[iS]);

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
      
	TGraph gRoc(nB, yCum_znn, yCum_qcd);
	gRoc.SetMarkerStyle(kOpenSquare);
	gRoc.SetMarkerColor(kRed);
	//gRoc.SetMarkerSize();
	gRoc.SetFillColor(kWhite);
	gRoc.GetXaxis()->Set(nB, 0, 1);
	gRoc.GetYaxis()->Set(nB, 0, 1);
	gRoc.GetXaxis()->SetTitle("Z(#nu#nu) efficiency");
	gRoc.GetYaxis()->SetTitle("QCD efficiency");
	gRoc.SetTitle("ROC Curve : "+select[iS]+" "+var[iV]+" "+wdT[iWd]);
	gRoc.Draw("AP");

	cRoc.Print("plots/"+_tag+"/roc_"+select[iS]+"_"+var[iV]+"_"+wd[iWd]+".pdf","pdf");

      } // end loop over nWd
    }   // end loop over nV
  }     // end loop over nS

  return 0;
}
