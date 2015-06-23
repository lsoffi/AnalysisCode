#ifndef MYINCLUDES
#define MYINCLUDES

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>
//
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEntryList.h"
#include "TMath.h"
//
#include "tdrstyle.h"

using namespace std;

pair<Double_t, Double_t> Integrate(TH1F* h) 
{

  if(!h) {
    cout << "ERROR: called Integrate on a NULL histogram pointer... Abort."
	 << endl;
    return make_pair(-888,-777);
  }

  Double_t norm  = h->Integral(0, h->GetNbinsX() + 1);
  Double_t ent   = h->GetEntries();
  Double_t error = ent>0 ? TMath::Sqrt(ent) * (norm/ent) : -999;

  return make_pair(norm,error);
}

Int_t setStyle(TH1F* h, Int_t color)
{
  h->Sumw2();
  h->SetMarkerSize(0.5);
  h->SetMarkerStyle(kPlus);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  //h->SetFillColor(color);
  return 0;
}

Int_t setStyle(TCanvas* c)
{
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  return 0;
}

Int_t setStyle(TLegend* leg)
{
  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0244755);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);  
  
  return 0;
}


#endif
