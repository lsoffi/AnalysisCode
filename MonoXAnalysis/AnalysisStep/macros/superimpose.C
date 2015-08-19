// #include <TTree.h>
// #include <TFile.h>
// #include <iostream>
// #include <stdio.h>
// #include <TEfficiency.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TF1.h>
// #include <TLegend.h>

#include "myIncludes.h"

Int_t superimpose(TString tag="v8_singlemuons_DCS_tunebin_TightMuon")
{

  const UInt_t nF=2;
  TString fitfunc[nF]={"Sigmoid","CB"};
  TString fitname[nF]={"fit2","fit"};

  TFile *file[nF];
  for(UInt_t iF=0 ; iF<nF ; iF++) {
    file[iF] = new TFile("results/"+tag+"_"+fitfunc[iF]+"/f_"+tag+"_"+fitfunc[iF]+".root", "READ");
  }

  const UInt_t nV=4; // mumet, t1mumet, signaljetpt, signaljetNHfrac
  const UInt_t nP=2; // 90GeV, 120GeV
  const UInt_t nS=8; // L1, MET, METClean, METJetID, MHT, PFMHT, PFMET, HLT bit

  Int_t color[nS] = {kBlack, kBlue+2, kBlue, kCyan+2, kGreen+2, kRed+2, kRed, kRed};
  Int_t style[nS] = {kOpenSquare, 
		     kOpenTriangleUp, kOpenTriangleDown, 
		     kFullTriangleUp, kFullTriangleDown, 
		     kOpenCircle, kFullCircle, kFullCircle};

  TString name, title;
  TString nameV[nV]={"mumet","t1mumet","signaljetpt","signaljetNHfrac"};
  TString nameP[nP]={"MET90","MET120"};
  TString nameS[nS]={"L1","MET","METClean","METJetID","MHT","PFMHT","PFMET","bit"};

  TEfficiency *pEff[nF][nV][nP][nS];
  TF1* fitEff[nF][nV][nP][nS];

  for(UInt_t iF=0 ; iF<nF ; iF++) {
    for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	  
	  name  = "t_h_"+nameV[iV]+"_num_"+nameP[iP]+"_"+nameS[iS];

	  pEff[iF][iV][iP][iS] = (TEfficiency*)file[iF]->Get(name);
	  fitEff[iF][iV][iP][iS] = 
	    (TF1*)(pEff[iF][iV][iP][iS]->GetListOfFunctions()->FindObject(fitname[iF]));

	  pEff[iF][iV][iP][iS]->SetLineColor(  color[iS]);
	  pEff[iF][iV][iP][iS]->SetMarkerColor(color[iS]);
	  pEff[iF][iV][iP][iS]->SetMarkerStyle(style[iS]);

	  fitEff[iF][iV][iP][iS]->SetLineColor(  color[iS]);
	  fitEff[iF][iV][iP][iS]->SetMarkerColor(color[iS]);
	  fitEff[iF][iV][iP][iS]->SetMarkerStyle(style[iS]);

	}
      }
    }
  }

  // PRODUCE PLOTS //
  const UInt_t nSPlot = 7;
  UInt_t idxSPlot[  nSPlot] = {0,1,2,3,4,5,7};
  UInt_t idxFitStep[nSPlot] = {0,0,0,0,0,0,0};

  const UInt_t nVPlot = 2;
  UInt_t idxVPlot[nVPlot] = {0,1};

  setTDRStyle();
  TCanvas *c;

  TString namePlot;
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iVPlot=0 ; iVPlot<nVPlot ; iVPlot++) {
      
      UInt_t iV=idxVPlot[iVPlot];
      
      namePlot = "plot_"+nameV[iV]+"_"+nameP[iP];
      c = new TCanvas("c_"+namePlot,"c_"+namePlot,0,0,600,600);
      
      for(UInt_t iSPlot=0 ; iSPlot<nSPlot ; iSPlot++) {
	UInt_t iS = idxSPlot[iSPlot];
	UInt_t iF = idxFitStep[iSPlot];
	pEff[iF][iV][iP][iS]->Draw("SAME");
      }

      c->Print("plots/"+tag+"/"+namePlot+".pdf");
    }
  }

  /*  
  TCanvas *c = new TCanvas();
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.25);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  //gStyle->SetStatTextColor(4);
  gStyle->SetStatTextColor(kRed);
  c->Update();
  eff_HCAL2->Draw("SAME");
  gStyle->SetStatY(0.4);
  //gStyle->SetStatTextColor(3);
  gStyle->SetStatTextColor(kBlue);
  c->Update();
  eff_HCAL3->Draw("SAME");
  gStyle->SetStatY(0.55);
  //gStyle->SetStatTextColor(2);
  gStyle->SetStatTextColor(kGreen+2);
  c->Update();
  
  TLegend *leg = new TLegend(0.58,0.58,0.73,0.68);
  leg->SetFillColor(kWhite);
  leg->AddEntry(eff_HCAL3_fit,"14e33_HCAL3","l");
  leg->AddEntry(eff_HCAL2_fit,"14e33_HCAL2","l");
  leg->AddEntry(eff_file_fit,"14e33_file","l");
  leg->Draw();
  c->Update();
  
  TString output = "results/superimposed_" + cut + ".pdf";
  c->Print(output,"pdf");
  */

  return 0;
}
