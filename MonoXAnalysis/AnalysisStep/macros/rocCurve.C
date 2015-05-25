#include <fstream>
#include <iostream>
#include <iomanip>
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

Int_t setStyle(TGraph* g, Int_t marker, Int_t color, 
	       UInt_t nB, TString title, Bool_t zoom);
Int_t setStyle(TLegend* leg);

Int_t rocCurve(TString _tag="", Float_t wpQCD=0.08, Float_t wpZNN=0.93, bool dolog=false)
{

  // Input file
  TFile* file = new TFile("plots/"+_tag+"/plots_"+_tag+".root","read");

  // Output log
  ofstream outlog("plots/"+_tag+"/wp_"+_tag+".txt",ios::out);
  outlog << setw(10) << "Selection"
	 << setw(14) << "Variable"
	 << setw(12) << "Cut"
	 << setw(7 ) << "wpQCD"
	 << setw(14) << "wpQCD_cut"
	 << setw(14) << "effQCD_wpQCD"
	 << setw(14) << "effZNN_wpQCD"
	 << setw(7 ) << "wpZNN"
	 << setw(14) << "wpZNN_cut"
	 << setw(14) << "effQCD_wpZNN"
	 << setw(14) << "effZNN_wpZNN"
	 << endl;
  
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
  const UInt_t nS=5;  // selections
  const UInt_t nV=5;  // variables
  const UInt_t nGZ=2;  // full/zoomed

  TGraph* gRoc[nS][nV][nWd][nGZ];

  TString tzoom[nGZ] = {"_full", "_zoom"}; 
  TString wd[nWd]    = {"upcut","lowcut"};
  TString wdT[nWd]   = {"(upper cut)","(lower cut)"};
  TString select[nS] = {"alljets","monojet","1jet","2jet","3jet"};
  TString var[nV]    = {"alphat","apcjetmetmax","apcjetmetmin","jetjetdphi","jetmetdphimin"};

  Int_t colors[nV]   = {kBlack, kBlue, kRed, kGreen+2, kMagenta};

  //UInt_t  nBins[nV]  = {40, 50, 50, 50,  50};
  //UInt_t  nBins[nV]  = {400, 500, 500, 500, 500};
  //UInt_t  nBins[nV]  = {800, 1000, 1000, 1000,  1000};
  //UInt_t  nBins[nV]  = {4000, 5000, 5000, 5000,  5000};
  //UInt_t  nBins[nV]  = {8000, 10000, 10000, 10000,  10000};
  
  //Float_t xFirst[nV] = {0,  0,  0,  0,   0  };
  //Float_t xLast[nV]  = {2,  1,  1,  3.2, 3.2};

  Int_t  wpQCD_bin,wpZNN_bin;
  Float_t wpQCD_cut,wpZNN_cut;
  Float_t effQCD_wpQCD, effZNN_wpQCD, effQCD_wpZNN, effZNN_wpZNN;
  wpQCD_bin = wpZNN_bin = wpQCD_cut = wpZNN_cut = -1;
  effQCD_wpQCD = effZNN_wpQCD = effQCD_wpZNN = effZNN_wpZNN = -1;
  
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iV=0 ; iV<nV ; iV++) {

      // Get input QCD and Znn distributions from the file
      TH1F* h_qcd = (TH1F*) file->Get("h_"+var[iV]+"_qcd_"+select[iS]);
      TH1F* h_znn = (TH1F*) file->Get("h_"+var[iV]+"_znn_"+select[iS]);

      // Check that the histograms were found
      if(!h_qcd) {
	cout << "ERROR: missing histogram" << endl
	     << "h_"+var[iV]+"_qcd_"+select[iS] << endl;
      }
      if(!h_znn) {
	cout << "ERROR: missing histogram" << endl
	     << "h_"+var[iV]+"_znn_"+select[iS] << endl;
      }
      if(!h_qcd || !h_znn) return -1;

      // Print the integral of the histograms
      cout << "-- "+var[iV]+" "+select[iS]+" : Integrals : "
	   << "qcd=" << h_qcd->Integral() << " "
	   << "znn=" << h_znn->Integral() << " " << endl;

      // Get the histograms binning
      //const UInt_t nB = nBins[iV];
      const UInt_t nB = h_qcd->GetNbinsX();

      // Build the graph from cumulative distributions
      // for each backward/forward case
      for(UInt_t iWd=0 ; iWd<nWd ; iWd++) {
	Double_t cum_qcd=0, cum_znn=0;
	Double_t yCum_qcd[nB];
	Double_t yCum_znn[nB];

	// Prepare full version of the graph
	/// Loop over the bins of the input distributions
	for(UInt_t iB=0 ; iB<nB ; iB++) {
	  cum_qcd += h_qcd->GetBinContent(iB+1);
	  cum_znn += h_znn->GetBinContent(iB+1);
	  if(iWd==0) { // upper cut : x<cut
	    yCum_qcd[iB] = cum_qcd;
	    yCum_znn[iB] = cum_znn;
	  }
	  else if(iWd==1) { // lower cut : x>cut
	    yCum_qcd[iB] = 1 - cum_qcd;
	    yCum_znn[iB] = 1 - cum_znn;
	  }

	  // Find bin index corresponding to requested QCD/ZNN efficiency (wpQCD/wpZNN)
	  if(iB>0) {
	    if(     iWd==0) { // upper cut x<cut
	      if(yCum_qcd[iB-1]<wpQCD && yCum_qcd[iB]>=wpQCD) wpQCD_bin = iB;
	      if(yCum_znn[iB-1]<wpZNN && yCum_znn[iB]>=wpZNN) wpZNN_bin = iB;
	    }
	    else if(iWd==1) { // lower cut x>cut
	      if(yCum_qcd[iB-1]>=wpQCD && yCum_qcd[iB]<wpQCD) wpQCD_bin = iB;
	      if(yCum_znn[iB-1]>=wpZNN && yCum_znn[iB]<wpZNN) wpZNN_bin = iB;
	    }
	  }
	}

	///////////////////////////////////////////////////////////////////////////////
	// WORKING POINTS /////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	// Find cut value using bin index for requested QCD/ZNN eff
	wpQCD_cut = h_qcd->GetXaxis()->GetBinCenter(wpQCD_bin);
	wpZNN_cut = h_znn->GetXaxis()->GetBinCenter(wpZNN_bin);
	//
	// Determine corresponding efficiencies
	// => precise value, for both QCD and ZNN, for each requested efficiency
	effQCD_wpQCD = (wpQCD_bin>=0 && wpQCD_bin<nB) ? yCum_qcd[wpQCD_bin] : -888888;
	effZNN_wpQCD = (wpQCD_bin>=0 && wpQCD_bin<nB) ? yCum_znn[wpQCD_bin] : -888888;
	effQCD_wpZNN = (wpZNN_bin>=0 && wpZNN_bin<nB) ? yCum_qcd[wpZNN_bin] : -888888;
	effZNN_wpZNN = (wpZNN_bin>=0 && wpZNN_bin<nB) ? yCum_znn[wpZNN_bin] : -888888;
	//
	// Write out values
	outlog << setw(10 ) << select[iS]
	       << setw(14) << var[iV]
	       << setw(12) << wdT[iWd]
	       << setw(7 ) << wpQCD
	       << setw(14) << wpQCD_cut
	       << setw(14) << effQCD_wpQCD
	       << setw(14) << effZNN_wpQCD
	       << setw(7 ) << wpZNN
	       << setw(14) << wpZNN_cut
	       << setw(14) << effQCD_wpZNN
	       << setw(14) << effZNN_wpZNN
	       << endl;
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	// GRAPHS /////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	/// use the arrays to fill the graph and draw it
	gRoc[iS][iV][iWd][0] = new TGraph(nB, yCum_znn, yCum_qcd);
	setStyle(gRoc[iS][iV][iWd][0], kOpenSquare, colors[iV], nB,
		 "ROC Curve : "+select[iS]+" "+var[iV]+" "+wdT[iWd], 
		 kFALSE);
	if(gRoc[iS][iV][iWd][0]) gRoc[iS][iV][iWd][0]->Draw("AL");	
	///
	/// print the graph in a single pdf file
	TString pdftitle = "Title:"+select[iS]+" "+var[iV]+" "+wd[iWd];
	if(iS==0 && iV==0 && iWd==0) 
	  cRoc.Print("plots/"+_tag+"/rocfull.pdf(",pdftitle);
	else if(iS==nS-1 && iV==nV-1 && iWd==nWd-1)
	  cRoc.Print("plots/"+_tag+"/rocfull.pdf)",pdftitle);
	else
	  cRoc.Print("plots/"+_tag+"/rocfull.pdf",pdftitle);
	// end full version

	// Prepare zoomed version by extracting interesting region
	/// determine indices of points where eff(znn)>60%
	vector<UInt_t> idxZoom;
	idxZoom.clear();
	for(UInt_t iB=0 ; iB<nB ; iB++) {
	  if(yCum_znn[iB]>0.60) idxZoom.push_back(iB);
	}
	///
	/// use the indices to build arrays
	const UInt_t nZ = idxZoom.size();
	Double_t yZoo_qcd[nZ];
	Double_t yZoo_znn[nZ];
	for(UInt_t iZ=0 ; iZ<nZ ; iZ++) {	
	  yZoo_qcd[iZ] = yCum_qcd[idxZoom[iZ]];
	  yZoo_znn[iZ] = yCum_znn[idxZoom[iZ]];
	}
	///
	/// use the arrays to fill the graph and draw it
	gRoc[iS][iV][iWd][1] = new TGraph(nZ, yZoo_znn, yZoo_qcd);
	setStyle(gRoc[iS][iV][iWd][1], kOpenSquare, colors[iV], nZ,
		 "ROC Curve : "+select[iS]+" "+var[iV]+" "+wdT[iWd],
		 kTRUE);
	if(gRoc[iS][iV][iWd][1]) gRoc[iS][iV][iWd][1]->Draw("AL");	
	///
	/// print the graph in a single pdf file
	if(iS==0 && iV==0 && iWd==0) 
	  cRoc.Print("plots/"+_tag+"/roczoom.pdf(",pdftitle);
	else if(iS==nS-1 && iV==nV-1 && iWd==nWd-1)
	  cRoc.Print("plots/"+_tag+"/roczoom.pdf)",pdftitle);
	else
	  cRoc.Print("plots/"+_tag+"/roczoom.pdf",pdftitle);
	// end zoomed version
	///////////////////////////////////////////////////////////////////////////////

      } // end loop over nWd
    }   // end loop over nV
  }     // end loop over nS

  // Put several killers per plot
  for(UInt_t iGZ=0 ; iGZ<nGZ ; iGZ++) {
      
    for(UInt_t iS=0 ; iS<nS ; iS++) {

      gRoc[iS][0][1][iGZ]->SetTitle("Selection: "+select[iS]);

      if(gRoc[iS][0][1][iGZ]->GetN()>0) gRoc[iS][0][1][iGZ]->Draw("AL"); // alphat
      else if(gRoc[iS][1][1][iGZ]->GetN()>0) gRoc[iS][1][1][iGZ]->Draw("AL");  // apcjetmetmax
      else if(gRoc[iS][2][1][iGZ]->GetN()>0) gRoc[iS][2][1][iGZ]->Draw("AL");  // apcjetmetmin
      else if(gRoc[iS][3][0][iGZ]->GetN()>0) gRoc[iS][3][0][iGZ]->Draw("AL");  // jetjetdphi
      else if(gRoc[iS][4][1][iGZ]->GetN()>0) gRoc[iS][4][1][iGZ]->Draw("AL");  // jetmetdphimin

      if(gRoc[iS][0][1][iGZ]->GetN()>0) gRoc[iS][0][1][iGZ]->Draw("L"); // alphat
      if(gRoc[iS][1][1][iGZ]->GetN()>0) gRoc[iS][1][1][iGZ]->Draw("L");  // apcjetmetmax
      if(gRoc[iS][2][1][iGZ]->GetN()>0) gRoc[iS][2][1][iGZ]->Draw("L");  // apcjetmetmin
      if(gRoc[iS][3][0][iGZ]->GetN()>0) gRoc[iS][3][0][iGZ]->Draw("L");  // jetjetdphi
      if(gRoc[iS][4][1][iGZ]->GetN()>0) gRoc[iS][4][1][iGZ]->Draw("L");  // jetmetdphimin

      TLegend* leg = new TLegend(0.38,0.6,0.58,0.8,"","brNDC");
      setStyle(leg);
    
      for(UInt_t iV=0 ; iV<nV ; iV++)
	leg->AddEntry(gRoc[iS][iV][1][iGZ], var[iV], "L");

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
	cRoc.Print("plots/"+_tag+"/allroc"+tzoom[iGZ]+".pdf(","Title:"+mysel);
      else if(iS==nS-1)
	cRoc.Print("plots/"+_tag+"/allroc"+tzoom[iGZ]+".pdf)","Title:"+mysel);
      else
	cRoc.Print("plots/"+_tag+"/allroc"+tzoom[iGZ]+".pdf","Title:"+mysel);

    }
  }
  
  return 0;
}

Int_t setStyle(TGraph* g, Int_t marker, Int_t color, 
	       UInt_t nB, TString title, Bool_t zoom)
{

  g->SetMarkerStyle(marker);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  //g->SetMarkerSize();
  g->SetFillColor(kWhite);

  g->GetXaxis()->Set(nB, 0, 1);
  g->GetYaxis()->Set(nB, 0, 1);
  if(zoom) g->GetXaxis()->Set(nB, 0.6, 1);

  g->SetMinimum(0);
  g->SetMaximum(1.05);

  g->GetXaxis()->SetTitle("Z(#nu#nu) efficiency");
  g->GetYaxis()->SetTitle("QCD efficiency");
  g->SetTitle(title);

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
