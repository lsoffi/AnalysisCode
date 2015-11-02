#include "myIncludes.h"

bool FAST=false; // fixme
int  NFAST=1000;

using namespace std;
typedef map<TString, map<TString, map<TString, TH1F*> > > M_SEL_PROC_VAR_H1D;
typedef map<TString, map<TString, map<TString, TH2F*> > > M_SEL_PROC_VAR_H2D;

typedef map<TString,TChain*> M_PROC_CH;

typedef map<TString, TH1F*>                M_VAR_H1D;
typedef map<TString, M_VAR_H1D >           M_CUT_VAR_H1D;
typedef map<TString, M_CUT_VAR_H1D >       M_PROC_CUT_VAR_H1D;
typedef map<TString, M_PROC_CUT_VAR_H1D >  M_SEL_PROC_CUT_VAR_H1D;

typedef map<TString, TH2F*>                M_VAR_H2D;
typedef map<TString, M_VAR_H2D >           M_CUT_VAR_H2D;
typedef map<TString, M_CUT_VAR_H2D >       M_PROC_CUT_VAR_H2D;
typedef map<TString, M_PROC_CUT_VAR_H2D >  M_SEL_PROC_CUT_VAR_H2D;

// Declare functions //
Int_t GetNGen(TString process);
Int_t writeHistos(TString tag, M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D);
Int_t doStackPlots(TString tag, TString region1D, M_SEL_PROC_CUT_VAR_H1D histos1D, const UInt_t nS, const UInt_t nC, TString* selection, TString* extra);
Int_t doTransferPlots(TString tag, TString region2D, M_SEL_PROC_CUT_VAR_H2D histos2D, const UInt_t nS, const UInt_t nC, TString* selection, TString* extra);
Int_t defineChains(M_PROC_CH &chains, TString region1D, TString region2D);
TCut  defineWeight(TString process);
TCut  defineCut(TString sample, TString region);
Int_t defineHistos1D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H1D &histos1D,
		     const UInt_t nS,    const UInt_t nC,
		     TString* selection, TString* extra);
Int_t defineHistos2D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H2D &histos2D,
		     const UInt_t nS,    const UInt_t nC,
		     TString* selection, TString* extra);

Int_t analyze(Float_t lumi1D, Float_t lumi2D, TString tag, TString region1D, TString region2D,
	      M_PROC_CH chains, M_SEL_PROC_CUT_VAR_H1D &histos1D, M_SEL_PROC_CUT_VAR_H2D &histos2D,
	      const UInt_t nS,    const UInt_t nC, TString* selection, TString* extra);

vector<TString> GetVarName(M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D, TString mode);
TGraphErrors TransferFactor(TH2F *hTemp2, Float_t cut);

// MAIN //
Int_t qcd(TString tag, TString region1D="HT", TString region2D="HT")
{

  cout << "- Start qcd study" << endl;

  //Float_t lumi=0.210;  // OUR PREVIOUS HYPOTHESIS
  Float_t lumi   =0.135210; // BRILCALC Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt 256630-257599
  Float_t lumiMET=0.109329; // BRILCALC Same + "HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v*"
  //Float_t lumiHT =0.000172; // BRILCALC 05Oct2015 256630-258158 HLT_PFHT200_v*
  //Float_t lumiHT =0.000470; // BRILCALC 05Oct2015+PRV4 256630-259862 HLT_PFHT200_v*
  Float_t lumiHT = 0.553150; // BRILCALC 05Oct2015 256630-258158 No trigger <=> Nadir_26Oct2015 trees

  Float_t lumi1D, lumi2D;
  if(region1D=="HT") lumi1D=lumiHT;
  else               lumi1D=lumiMET;
  if(region2D=="HT") lumi2D=lumiHT;
  else               lumi2D=lumiMET;

  cout << "- lumi1D=" << lumi1D << " lumi2D=" << lumi2D << endl;

  TString dirOut="/user/ndaci/Results/Monojet/QCD/"+tag+"/";
  ofstream outlog(dirOut+"yields_"+tag+".txt",ios::out);

  // Define Chains
  M_PROC_CH chains;
  cout << "- defineChains(chains," << region1D << "," << region2D << ");... " ;
  defineChains(chains, region1D, region2D);
  cout << " done!" << endl;

  // Define selections to be explored
  const UInt_t nS=4;
  TString selection[nS] = {"inclusive","1jet","2jet","3jet"};
  //const UInt_t nS=2;
  //TString selection[nS] = {"inclusive","1jet"};

  const UInt_t nC=2;
  TString extra[nC]={"None","FwdVeto"};

  // Define histograms
  M_SEL_PROC_CUT_VAR_H1D histos1D;
  if(tag.Contains("1D")) {
    cout << "- defineHistos1D(chains, histos1D, " << nS << "," << nC << ",selection,extra);...";
    defineHistos1D(chains, histos1D, nS, nC, selection, extra);
    cout << " done!" << endl;
  }

  M_SEL_PROC_CUT_VAR_H2D histos2D;
  if(tag.Contains("2D")) {
    cout << "- defineHistos1D(chains, histos2D, " << nS << "," << nC << ",selection,extra);...";
    defineHistos2D(chains, histos2D, nS, nC, selection, extra);
    cout << " done!" << endl;
  }

  // Process input chains
  cout << "- analyze(" << lumi1D << ", " << lumi2D << ", " << tag << ", " << region1D << ", " << region2D 
       << ",chains, histos1D, histos2D, " << nS << ", " << nC << ", selection, extra);...";
  analyze(lumi1D, lumi2D, tag, region1D, region2D, chains, histos1D, histos2D, nS, nC, selection, extra);
  cout << " done!" << endl;

  // Write out the histograms
  cout << "- writeHistos(" << tag << ", histos1D, histos2D);...";
  writeHistos(tag, histos1D, histos2D);
  cout << " done!" << endl;

  // Produce the 1D stack plots
  cout << "- doStackPlots(" << tag << ", " << region1D << ", histos1D, " << nS << ", " << nC << ", selection, extra);" ;
  doStackPlots(   tag, region1D, histos1D, nS, nC, selection, extra);
  cout << " done!" << endl;

  // Produce the 2D scatter plots and transfer plots
  cout << "- doTransferPlots(" << tag << ", " << region2D << ", histos1D, " << nS << ", " << nC << ", selection, extra);" ;
  doTransferPlots(tag, region2D, histos2D, nS, nC, selection, extra);
  cout << " done!" << endl;

  // END //
  cout << "--- THE END ==> [] ---" << endl;
  return 0;
}

Int_t doTransferPlots(TString tag, TString region2D, M_SEL_PROC_CUT_VAR_H2D histos2D, const UInt_t nS, const UInt_t nC, TString* selection, TString* extra)
{

  // Output
  TString dirOut="/user/ndaci/Results/Monojet/QCD/"+tag+"/";
  TFile* outfile = new TFile(dirOut+"/plots_"+tag+".root","recreate");
  outfile->cd();

  // Retrieve 1D variables
  M_SEL_PROC_CUT_VAR_H1D dummy;
  vector<TString> var = GetVarName(dummy, histos2D, "2D");
  const UInt_t nV=var.size();

  // Retrieve histograms  
  TString name,title;
  Int_t color,style; 
  Float_t size;
  TH2F *hData, *hDataClean, *hQCD, *hZLL, *hWLN, *hZNN, *hVV, *hTT, *hTop;
  TLegend* leg;

  // Loop over histos map
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iC=0 ; iC<nC ; iC++) {
      for(UInt_t iV=0 ; iV<nV ; iV++) {

	cout << "---- " << " " << selection[iS] << " " << extra[iC] << " " << var[iV] << endl;

	// Data
	name="h1D_"+selection[iS]+"_dataTotal_"+extra[iC]+"_"+var[iV];

	if(region2D.Contains("MET"))
	  hData = (TH2F*) histos2D[selection[iS]]["data_met"][extra[iC]][var[iV]]->Clone(name);
	else if(region2D.Contains("HT"))
	  hData = (TH2F*) histos2D[selection[iS]]["data_jetht"][extra[iC]][var[iV]]->Clone(name);
	else hData = (TH2F*) histos2D[selection[iS]]["data_jetht"][extra[iC]][var[iV]]->Clone(name);
	
	if(hData) {
	  hData->Write();
	  //setStyle(hData, kBlack, kFullCircle, 1.00, false);
	  cout << "---- hData: Entries=" << hData->GetEntries() << " Integral=" << hData->Integral() << endl;
	}

	// QCD
	name="h1D_"+selection[iS]+"_qcdTotal_"+extra[iC]+"_"+var[iV];
	hQCD=(TH2F*) histos2D[selection[iS]]["qcdht1000to1500"][extra[iC]][var[iV]]->Clone(name);
	if(hQCD) {
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht100to200"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht1500to2000"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht2000toinf"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht200to300"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht300to500"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht500to700"][extra[iC]][var[iV]]);
	  hQCD   ->Add(histos2D[selection[iS]]["qcdht700to1000"][extra[iC]][var[iV]]);
	  hQCD->Write();
	  //setStyle(hQCD, kRed, kOpenSquare, 1.00, true);
	  cout << "---- hQCD: Entries=" << hQCD->GetEntries() << " Integral=" << hQCD->Integral() << endl;
	}

	// ZLL
	name="h1D_"+selection[iS]+"_zllTotal_"+extra[iC]+"_"+var[iV];
	hZLL=(TH2F*) histos2D[selection[iS]]["zll100to200"][extra[iC]][var[iV]]->Clone(name);
	if(hZLL) {
	  hZLL   ->Add(histos2D[selection[iS]]["zll200to400"][extra[iC]][var[iV]]);
	  hZLL   ->Add(histos2D[selection[iS]]["zll400to600"][extra[iC]][var[iV]]);
	  hZLL   ->Add(histos2D[selection[iS]]["zll600toinf"][extra[iC]][var[iV]]);
	  //setStyle(hZLL, kPink+9, kOpenSquare, 2, true);
	  cout << "---- hZLL: Entries=" << hZLL->GetEntries() << " Integral=" << hZLL->Integral() << endl;
	  hZLL->Write();
	}

	// WLN
	name="h1D_"+selection[iS]+"_wlnTotal_"+extra[iC]+"_"+var[iV];
	hWLN=(TH2F*) histos2D[selection[iS]]["wln100to200"][extra[iC]][var[iV]]->Clone(name);
	if(hWLN) {
	  hWLN   ->Add(histos2D[selection[iS]]["wln200to400"][extra[iC]][var[iV]]);
	  hWLN   ->Add(histos2D[selection[iS]]["wln400to600"][extra[iC]][var[iV]]);
	  hWLN   ->Add(histos2D[selection[iS]]["wln600toinf"][extra[iC]][var[iV]]);
	  //setStyle(hWLN, kGreen+2, kOpenSquare, 2, true);
	  cout << "---- hWLN: Entries=" << hWLN->GetEntries() << " Integral=" << hWLN->Integral() << endl;
	  hWLN->Write();
	}

	// ZNN
	name="h1D_"+selection[iS]+"_znnTotal_"+extra[iC]+"_"+var[iV];
	hZNN=(TH2F*) histos2D[selection[iS]]["znn100to200"][extra[iC]][var[iV]]->Clone(name);
	if(hZNN) {
	  hZNN   ->Add(histos2D[selection[iS]]["znn200to400"][extra[iC]][var[iV]]);
	  hZNN   ->Add(histos2D[selection[iS]]["znn400to600"][extra[iC]][var[iV]]);
	  hZNN   ->Add(histos2D[selection[iS]]["znn600toinf"][extra[iC]][var[iV]]);
	  //setStyle(hZNN, kAzure+7, kOpenSquare, 2, true);
	  cout << "---- hZNN: Entries=" << hZNN->GetEntries() << " Integral=" << hZNN->Integral() << endl;
	  hZNN->Write();
	}

	// VV
	name="h1D_"+selection[iS]+"_vvTotal_"+extra[iC]+"_"+var[iV];
	hVV=(TH2F*) histos2D[selection[iS]]["ww"][extra[iC]][var[iV]]->Clone(name);
	if(hVV) {
	  hVV   ->Add(histos2D[selection[iS]]["zz"][extra[iC]][var[iV]]);
	  hVV   ->Add(histos2D[selection[iS]]["wz"][extra[iC]][var[iV]]);
	  //setStyle(hVV, kBlue+1, kOpenSquare, 2, true);
	  cout << "---- hVV: Entries=" << hVV->GetEntries() << " Integral=" << hVV->Integral() << endl;
	  hVV->Write();
	}

	// Single Top AND TTBAR
	name="h1D_"+selection[iS]+"_topTotal_"+extra[iC]+"_"+var[iV];
	hTop=(TH2F*) histos2D[selection[iS]]["singletbart"][extra[iC]][var[iV]]->Clone(name);
	if(hTop) {
	  hTop   ->Add(histos2D[selection[iS]]["singletbarw"][extra[iC]][var[iV]]);
	  hTop   ->Add(histos2D[selection[iS]]["singlett"][extra[iC]][var[iV]]);
	  hTop   ->Add(histos2D[selection[iS]]["singletw"][extra[iC]][var[iV]]);
	  hTop   ->Add(histos2D[selection[iS]]["ttbar"][extra[iC]][var[iV]]);
	  //setStyle(hTop, kOrange-3, kOpenSquare, 2, true);
	  cout << "---- hTop: Entries=" << hTop->GetEntries() << " Integral=" << hTop->Integral() << endl;
	  hTop->Write();
	}

	// Data Cleaned
	TString nameData=TString(hData->GetName());
	hDataClean = (TH2F*) hData->Clone(nameData+"_clean");
	if(hDataClean) {
	  hDataClean->Add(hZNN, -1.0);
	  hDataClean->Add(hWLN, -1.0);
	  hDataClean->Add(hZLL, -1.0);
	  hDataClean->Add(hVV, -1.0);
	  hDataClean->Add(hTop, -1.0);
	  hDataClean->Write();
	}

	// Produce plot //
	//
	/// canvas
	TCanvas c("c","c",20,20,600,600);
	gPad->SetLogy();
	gStyle->SetOptStat(0);

	/// Print
	name = "scat_"+selection[iS]+"_"+extra[iC]+"_"+var[iV];
	if(hData) {
	  hData->Draw("COLZ");
	  c.Print(dirOut+"/"+name+"_data.pdf","pdf");
	}
	if(hDataClean) {
	  hDataClean->Draw("COLZ");
	  c.Print(dirOut+"/"+name+"_dataclean.pdf","pdf");
	}
	if(hQCD) {
	  hQCD->Draw("COLZ");
	  c.Print(dirOut+"/"+name+"_qcd.pdf","pdf");
	}

	// Produce Transfer Factors plots
	TGraphErrors gSF_data, gSF_qcd;
	Float_t cut=0.5;

	/// DATA TF 
	color=kBlack;
	style=kFullCircle;
	size=1.00;
	name  = "tf_"+selection[iS]+"_"+extra[iC]+"_"+var[iV]+"_data";
	gSF_data = TransferFactor(hDataClean, cut);
	gSF_data.SetName(name);
	gSF_data.SetTitle("QCD Transfer Factor (data)");
	gSF_data.SetMinimum(0.0001);
	gSF_data.SetMaximum(1000);
	gSF_data.SetMarkerSize(size);
	gSF_data.SetMarkerStyle(style);
	gSF_data.SetMarkerColor(color);
	gSF_data.SetLineColor(  color);
	/// print data	
	gSF_data.Draw("AP");
	c.Print(dirOut+"/"+name+".pdf","pdf");
	/// print zoome data	
	gSF_data.GetXaxis()->SetRangeUser(40, 200);
	gSF_data.GetXaxis()->SetTitle("PFMETNoMu [GeV]");
	gSF_data.Draw("AP");
	name = "sfzoom_"+selection[iS]+"_"+extra[iC]+"_"+var[iV]+"_data";
	c.Print(dirOut+"/"+name+".pdf","pdf");

	/// QCD TF 
	color=kRed;
	style=kOpenSquare;
	size=1.00;
	name  = "tf_"+selection[iS]+"_"+extra[iC]+"_"+var[iV]+"_qcd";
	gSF_qcd = TransferFactor(hQCD, cut);
	gSF_qcd.SetName(name);
	gSF_qcd.SetTitle("QCD Transfer Factor (QCD MC)");
	gSF_qcd.SetMinimum(0.0001);
	gSF_qcd.SetMaximum(1000);
	gSF_qcd.SetMarkerSize(size);
	gSF_qcd.SetMarkerStyle(style);
	gSF_qcd.SetMarkerColor(color);
	gSF_qcd.SetLineColor(  color);
	/// print qcd	
	gSF_qcd.Draw("AP");
	c.Print(dirOut+"/"+name+".pdf","pdf");
	/// print zoome qcd	
	gSF_qcd.GetXaxis()->SetRangeUser(40, 200);
	gSF_qcd.GetXaxis()->SetTitle("PFMETNoMu [GeV]");
	gSF_qcd.Draw("AP");
	name = "sfzoom_"+selection[iS]+"_"+extra[iC]+"_"+var[iV]+"_qcd";
	c.Print(dirOut+"/"+name+".pdf","pdf");

	// Print both
	gSF_data.SetTitle("QCD Transfer Factor");
	gSF_qcd .SetTitle("QCD Transfer Factor");
	gSF_data.Draw("AP");
	gSF_qcd .Draw("PSAME");
	
	TLegend* leg = new TLegend(0.7,0.79,0.89,0.89,"","brNDC");
	setStyle(leg);
	leg->AddEntry(&gSF_data, "Data (JetHT)", "P");
	leg->AddEntry(&gSF_qcd,  "QCD  (MC)"   , "P");
	leg->Draw();
	
	name = "sfboth_"+selection[iS]+"_"+extra[iC]+"_"+var[iV];
	c.Print(dirOut+"/"+name+".pdf","pdf");

      } // end loop: var
    } // end loop: extra cuts
  } // end loop: selections

  return 0;
}

Int_t doStackPlots(TString tag, TString region1D, M_SEL_PROC_CUT_VAR_H1D histos1D, const UInt_t nS, const UInt_t nC, TString* selection, TString* extra)
{

  // Output
  TString dirOut="/user/ndaci/Results/Monojet/QCD/"+tag+"/";

  // Retrieve 1D variables
  M_SEL_PROC_CUT_VAR_H2D dummy;
  vector<TString> var = GetVarName(histos1D, dummy, "1D");
  const UInt_t nV=var.size();

  // Retrieve histograms  
  TString name,title;
  TH1F *hData, *hQCD, *hZLL, *hWLN, *hZNN, *hVV, *hTT, *hTop;
  THStack* stack;
  TLegend* leg;

  // Loop over histos map
  for(UInt_t iS=0 ; iS<nS ; iS++) {
    for(UInt_t iC=0 ; iC<nC ; iC++) {
      for(UInt_t iV=0 ; iV<nV ; iV++) {

	cout << "---- " << " " << selection[iS] << " " << extra[iC] << " " << var[iV] << endl;

	// Data
	name="h1D_"+selection[iS]+"_dataTotal_"+extra[iC]+"_"+var[iV];

	if(region1D.Contains("MET"))
	  hData = (TH1F*) histos1D[selection[iS]]["data_met"][extra[iC]][var[iV]]->Clone(name);
	else if(region1D.Contains("HT"))
	  hData = (TH1F*) histos1D[selection[iS]]["data_jetht"][extra[iC]][var[iV]]->Clone(name);
	else hData = (TH1F*) histos1D[selection[iS]]["data_jetht"][extra[iC]][var[iV]]->Clone(name);

	setStyle(hData, kBlack, kFullCircle, 1.00, false);
	cout << "---- hData: Entries=" << hData->GetEntries() << " Integral=" << hData->Integral() << endl;

	// QCD
	name="h1D_"+selection[iS]+"_qcdTotal_"+extra[iC]+"_"+var[iV];
	hQCD=(TH1F*) histos1D[selection[iS]]["qcdht1000to1500"][extra[iC]][var[iV]]->Clone(name);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht100to200"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht1500to2000"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht2000toinf"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht200to300"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht300to500"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht500to700"][extra[iC]][var[iV]]);
	hQCD   ->Add(histos1D[selection[iS]]["qcdht700to1000"][extra[iC]][var[iV]]);
	setStyle(hQCD, kRed, kOpenSquare, 1.00, true);
	cout << "---- hQCD: Entries=" << hQCD->GetEntries() << " Integral=" << hQCD->Integral() << endl;

	// ZLL
	name="h1D_"+selection[iS]+"_zllTotal_"+extra[iC]+"_"+var[iV];
	hZLL=(TH1F*) histos1D[selection[iS]]["zll100to200"][extra[iC]][var[iV]]->Clone(name);
	hZLL   ->Add(histos1D[selection[iS]]["zll200to400"][extra[iC]][var[iV]]);
	hZLL   ->Add(histos1D[selection[iS]]["zll400to600"][extra[iC]][var[iV]]);
	hZLL   ->Add(histos1D[selection[iS]]["zll600toinf"][extra[iC]][var[iV]]);
	setStyle(hZLL, kPink+9, kOpenSquare, 2, true);
	cout << "---- hZLL: Entries=" << hZLL->GetEntries() << " Integral=" << hZLL->Integral() << endl;

	// WLN
	name="h1D_"+selection[iS]+"_wlnTotal_"+extra[iC]+"_"+var[iV];
	hWLN=(TH1F*) histos1D[selection[iS]]["wln100to200"][extra[iC]][var[iV]]->Clone(name);
	hWLN   ->Add(histos1D[selection[iS]]["wln200to400"][extra[iC]][var[iV]]);
	hWLN   ->Add(histos1D[selection[iS]]["wln400to600"][extra[iC]][var[iV]]);
	hWLN   ->Add(histos1D[selection[iS]]["wln600toinf"][extra[iC]][var[iV]]);
	setStyle(hWLN, kGreen+2, kOpenSquare, 2, true);
	cout << "---- hWLN: Entries=" << hWLN->GetEntries() << " Integral=" << hWLN->Integral() << endl;

	// ZNN
	name="h1D_"+selection[iS]+"_znnTotal_"+extra[iC]+"_"+var[iV];
	hZNN=(TH1F*) histos1D[selection[iS]]["znn100to200"][extra[iC]][var[iV]]->Clone(name);
	hZNN   ->Add(histos1D[selection[iS]]["znn200to400"][extra[iC]][var[iV]]);
	hZNN   ->Add(histos1D[selection[iS]]["znn400to600"][extra[iC]][var[iV]]);
	hZNN   ->Add(histos1D[selection[iS]]["znn600toinf"][extra[iC]][var[iV]]);
	setStyle(hZNN, kAzure+7, kOpenSquare, 2, true);
	cout << "---- hZNN: Entries=" << hZNN->GetEntries() << " Integral=" << hZNN->Integral() << endl;
	
	// VV
	name="h1D_"+selection[iS]+"_vvTotal_"+extra[iC]+"_"+var[iV];
	hVV=(TH1F*) histos1D[selection[iS]]["ww"][extra[iC]][var[iV]]->Clone(name);
	hVV   ->Add(histos1D[selection[iS]]["zz"][extra[iC]][var[iV]]);
	hVV   ->Add(histos1D[selection[iS]]["wz"][extra[iC]][var[iV]]);
	setStyle(hVV, kBlue+1, kOpenSquare, 2, true);
	cout << "---- hVV: Entries=" << hVV->GetEntries() << " Integral=" << hVV->Integral() << endl;

	// Single Top AND TTBAR
	name="h1D_"+selection[iS]+"_topTotal_"+extra[iC]+"_"+var[iV];
	hTop=(TH1F*) histos1D[selection[iS]]["singletbart"][extra[iC]][var[iV]]->Clone(name);
	hTop   ->Add(histos1D[selection[iS]]["singletbarw"][extra[iC]][var[iV]]);
	hTop   ->Add(histos1D[selection[iS]]["singlett"][extra[iC]][var[iV]]);
	hTop   ->Add(histos1D[selection[iS]]["singletw"][extra[iC]][var[iV]]);
	hTop   ->Add(histos1D[selection[iS]]["ttbar"][extra[iC]][var[iV]]);
	setStyle(hTop, kOrange-3, kOpenSquare, 2, true);
	cout << "---- hTop: Entries=" << hTop->GetEntries() << " Integral=" << hTop->Integral() << endl;

	// Produce plot //
	//
	/// canvas
	TCanvas c("c","c",20,20,600,600);
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	leg = new TLegend(0.7,0.7,0.89,0.89,"","brNDC");
	setStyle(leg);
	leg->AddEntry(hZNN, "Z(#nu#nu)", "F");
	leg->AddEntry(hWLN, "W(L#nu)", "F");
	leg->AddEntry(hZLL, "Z(LL)", "F");
	leg->AddEntry(hTop, "Top", "F");
	leg->AddEntry(hVV,  "Diboson", "F");
	leg->AddEntry(hQCD, "QCD", "F");
	//
	/// MC Stack
	name  = "ths_"   +selection[iS]+"_"+extra[iC]+"_"+var[iV];
	title = "Stack: "+selection[iS]+" "+extra[iC]+" "+var[iV];
	stack = new THStack(name,title);
	stack->Add(hQCD);
	stack->Add(hVV);
	stack->Add(hTop);
	stack->Add(hZLL);
	stack->Add(hWLN);
	stack->Add(hZNN);
	//
	/// Print
	hData->SetMinimum(1.0);
	hData->Draw("P");
	stack->Draw("HISTSAME");
	hData->Draw("PSAME");
	leg->Draw();
	c.Print(dirOut+"/"+name+".pdf","pdf");
      }
    }
  }

  return 0;
}

Int_t writeHistos(TString tag, M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D)
{

  TString dirOut="/user/ndaci/Results/Monojet/QCD/";
  TFile* outfile = new TFile(dirOut+"/"+tag+"/plots_"+tag+".root","recreate");
  outfile->cd();

  M_SEL_PROC_CUT_VAR_H1D::iterator itSelProcCutVarH1D;
  M_PROC_CUT_VAR_H1D::iterator     itProcCutVarH1D;
  M_CUT_VAR_H1D::iterator          itCutVarH1D;
  M_VAR_H1D::iterator              itVarH1D;

  M_SEL_PROC_CUT_VAR_H2D::iterator itSelProcCutVarH2D;
  M_PROC_CUT_VAR_H2D::iterator     itProcCutVarH2D;
  M_CUT_VAR_H2D::iterator          itCutVarH2D;
  M_VAR_H2D::iterator              itVarH2D;

  for((itSelProcCutVarH1D=histos1D.begin()) ; 
      (itSelProcCutVarH1D!=histos1D.end()) ; 
      itSelProcCutVarH1D++) {
    for((itProcCutVarH1D=(itSelProcCutVarH1D->second).begin()) ; 
	(itProcCutVarH1D!=(itSelProcCutVarH1D->second).end()) ; 
	itProcCutVarH1D++) {
      for((itCutVarH1D=(itProcCutVarH1D->second).begin()) ; 
	  (itCutVarH1D!=(itProcCutVarH1D->second).end()) ; 
	  itCutVarH1D++) {
	for((itVarH1D=(itCutVarH1D->second).begin()) ; 
	    (itVarH1D!=(itCutVarH1D->second).end()) ; 
	    itVarH1D++) {
	  (itVarH1D->second)->Write();
	}
      }
    }
  }

  for((itSelProcCutVarH2D=histos2D.begin()) ; 
      (itSelProcCutVarH2D!=histos2D.end()) ; 
      itSelProcCutVarH2D++) {
    for((itProcCutVarH2D=(itSelProcCutVarH2D->second).begin()) ; 
	(itProcCutVarH2D!=(itSelProcCutVarH2D->second).end()) ; 
	itProcCutVarH2D++) {
      for((itCutVarH2D=(itProcCutVarH2D->second).begin()) ; 
	  (itCutVarH2D!=(itProcCutVarH2D->second).end()) ; 
	  itCutVarH2D++) {
	for((itVarH2D=(itCutVarH2D->second).begin()) ; 
	    (itVarH2D!=(itCutVarH2D->second).end()) ; 
	    itVarH2D++) {
	  (itVarH2D->second)->Write();
	}
      }
    }
  }

  outfile->Write();
  return 0;
}

Int_t analyze(Float_t lumi1D, Float_t lumi2D, TString tag, TString region1D, TString region2D,
	      M_PROC_CH chains, M_SEL_PROC_CUT_VAR_H1D &histos1D, M_SEL_PROC_CUT_VAR_H2D &histos2D,
	      const UInt_t nS,    const UInt_t nC, TString* selection, TString* extra)
{

  cout << "- START analyze()" << endl
       << "- Getting variable names: GetVarName(histos1D, histos2D, 'both');" << endl;

  TString mode;
  if(tag.Contains("1D") && tag.Contains("2D")) mode="both";
  else if(tag.Contains("1D"))                  mode="1D";
  else if(tag.Contains("2D"))                  mode="2D";
  else mode="1D";

  vector<TString> var = GetVarName(histos1D, histos2D, mode);
  const UInt_t nV=var.size();

  // Loop over chains
  TChain* ch;
  TEntryList* skim;
  TH1F* hTemp1;
  TH2F* hTemp2;
  TString process, tskim, select, hname, locVar, region;
  M_PROC_CH::iterator itCh;
  TCut theCut, theWeight;
  Float_t theScale;
  Int_t   nGen;
  //
  cout << "- Start loop: chains" << endl;
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {

    process = itCh->first;
    ch      = itCh->second;

    cout << "- " << process << endl;

    for(UInt_t iS=0 ; iS<nS ; iS++) {

      cout << "-- " << selection[iS] << endl;

      for(UInt_t iC=0 ; iC<nC ; iC++) {

	cout << "--- " << extra[iC] << endl;
	
	// Skim the chain
	select = selection[iS]+"_"+extra[iC];
	theCut = defineCut(process, select);
	theWeight = defineWeight(process);
	cout << "--- cut: "    << theCut    << endl
	     << "--- weight: " << theWeight << endl;
	//
	ch->SetEntryList(0);
	tskim  = "skim_"+process+"_"+select;
	if(FAST) ch->Draw(">>+"+tskim, theCut, "entrylist", NFAST);
	else     ch->Draw(">>+"+tskim, theCut, "entrylist");
	skim = (TEntryList*)gDirectory->Get(tskim);
	if(skim) {
	  ch->SetEntryList(skim);
	  cout << "--- skim size: " << skim->GetN() << " entries" << endl;
	}

	// Draw the variable
	for(UInt_t iV=0 ; iV<nV ; iV++) {

	  cout << "----- " << process << " " << selection[iS] << " " << extra[iC] << " " << var[iV] << endl;

	  // Translate variable names
	  locVar=var[iV];
	  if(     var[iV]=="jetmetdphimin_vs_t1mumet")    locVar="abs(jetmetdphimin):t1mumet";
	  else if(var[iV]=="incjetmetdphimin_vs_t1mumet") locVar="abs(incjetmetdphimin):t1mumet";
	  else if(var[iV].Contains("phi")) locVar = "abs("+var[iV]+")";

	  // Establish the rescaling factor to the final MC histogram
	  theScale = 1; // lumi is applied later depending on the requested region
	  nGen     = GetNGen(process);
	  if(!process.Contains("data")) {
	    theScale = nGen!=0 ? 1/nGen : 1;
	  }
	  if(process.Contains("zll") || process.Contains("znn")) {
	    theScale *= 1.23;
	  }
	  if(process.Contains("wln")) {
	    theScale *= 1.21;
	  }

	  // Launch calls to TTree::Draw()
	  if(var[iV].Contains("_vs_")) { // 2D
	    hTemp2 = histos2D[selection[iS]][process][extra[iC]][var[iV]];
	    hname  = hTemp2->GetName();
	    region = region2D;
	    theCut = defineCut(process, select+"_"+region); // "HT" => use HT triggers + offline HT cut
	    cout << "------ 2D: " << theCut << endl;
	    ch->Draw(locVar+">>"+hname, theCut*theWeight);
	    if(!process.Contains("data")) theScale *= lumi2D;
	    hTemp2->Scale(theScale);
	  }
	  else { // 1D
	    hTemp1 = histos1D[selection[iS]][process][extra[iC]][var[iV]];
	    hname  = hTemp1->GetName();
	    region = region1D;
	    if(!var[iV].Contains("phimin") && region1D.Contains("MET")) region+="_JetMet0p5";
	    theCut = defineCut(process, select+"_"+region); 
	    cout << "------ 1D: " << theCut << endl;
	    ch->Draw(locVar+">>"+hname, theCut*theWeight);
	    if(!process.Contains("data")) theScale *= lumi1D;
	    hTemp1->Scale(theScale);
	  }

	} // end loop: var
      } // end loop: cuts
    } // end loop: selection
  } // end loop: processes

  return 0;
}

TCut defineWeight(TString process)
{

  TCut weight="1";

  if(!process.Contains("data")) {
    //weight = "xsec*puweight*(wgt/wgtsum)";
    weight = "xsec*puweight"; // fixme
  }

  return weight;
}

Int_t defineHistos1D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H1D &histos1D,
		     const UInt_t nS,    const UInt_t nC, 
		     TString* selection, TString* extra)
{

  /// Stack 1D plots
  const UInt_t nV=6;
  TString var[nV]     ={"t1mumet", "jetmetdphimin", "incjetmetdphimin", "nvtx", "ht", "njets"};
  TString nameAxisX[nV]={"Type1 PFMETNoMu [GeV]", "Min #Delta#phi(M,J_{i}^{C})", "Min #Delta#phi(M,J_{i})",
			 "Number of Vertices", "Reco PFHT [GeV]" , "Number of Jets" };
  TString nameAxisY[nV]={"Events","Events","Events","Events","Events","Events"};

  UInt_t   nBins[nV]={ 100,   64,   64, 40,   50, 10};
  Float_t xFirst[nV]={ 200,    0,    0,  0,    0,  0};
  Float_t xLast[ nV]={1000,  3.2,  3.2, 40, 1000, 10};
  Bool_t regular[nV]={true,true,true,true,true,true};

  // Tuned binning
  Float_t* binsTuned[nV];  
  /// MET
  regular[0] = false;
  nBins[0] = 8;
  Float_t bins_met[8] = {200, 250, 300, 350, 400, 500, 600, 1000};
  binsTuned[0]=bins_met;

  vector<vector<Float_t>> x_bins;
  vector<Float_t> theBins;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {
    theBins.clear();
    for(UInt_t iB=0 ; iB<nBins[iV] ; iB++) {
      if(regular[iV]) {
	if(nBins[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBins[iV]);
	else             binval = 0;
	theBins.push_back(binval);
      }
      else {
	theBins.push_back(binsTuned[iV][iB]);
      }
    }
    x_bins.push_back(theBins);
  }

  TString name,title;

  // Get process names from chains map
  vector<TString> process;
  M_PROC_CH::iterator itCh;
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {
    process.push_back( itCh->first );
  }
  const UInt_t nP = process.size();

  // Create histos and fill the map
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      for(UInt_t iC=0 ; iC<nC ; iC++) {

	for(UInt_t iV=0 ; iV<nV ; iV++) { 
	  //for(UInt_t iV=1 ; iV<nV ; iV++) {

	  name ="h1D_"+selection[iS]+"_"+process[iP]+"_"+extra[iC]+"_"+var[iV];
	  title=selection[iS]+" "+process[iP]+" "+extra[iC]+" "+var[iV];

	  const UInt_t nBins = x_bins[iV].size();
	  Float_t theBins[nBins];
	  for(UInt_t iB=0 ; iB<nBins ; iB++) {
	    theBins[iB] = x_bins[iV][iB];
	  }

	  histos1D[selection[iS]][process[iP]][extra[iC]][var[iV]] = 
	    new TH1F(name, title, nBins-1, theBins);

	  histos1D[selection[iS]][process[iP]][extra[iC]][var[iV]]->Sumw2();	  
	}
      }
    }
  }

  return 0;
}

Int_t defineHistos2D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H2D &histos2D,
		     const UInt_t nS,    const UInt_t nC, 
		     TString* selection, TString* extra)
{

  /// Stack 2D plots
  const UInt_t nV=2;
  TString var[nV]      ={"jetmetdphimin_vs_t1mumet", "incjetmetdphimin_vs_t1mumet"};
  TString nameAxisX[nV]={"Type1 PFMETNoMu [GeV]", "Type1 PFMETNoMu [GeV]"};
  TString nameAxisY[nV]={"Min #Delta#phi(M,J_{i}^{C})", "Min #Delta#phi(M,J_{i})"};

  UInt_t   nBinsX[nV]={ 100,  100};  
  Float_t  xFirst[nV]={ 200,  200};
  Float_t  xLast[ nV]={1000, 1000};
  Bool_t regularX[nV]={true,true};

  UInt_t   nBinsY[nV]={64,   64};
  Float_t yFirst[ nV]={ 0,    0};
  Float_t yLast[  nV]={ 3.2,  3.2};
  Bool_t regularY[nV]={true,true};

  // Tuned binning
  Float_t* binsTunedX[nV];
  Float_t* binsTunedY[nV];    

  /// MET
  Float_t bins_met[17] = {50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 1000};
  regularX[0] = false;
  nBinsX[0] = 17;
  binsTunedX[0]=bins_met;
  // FIXME
  regularX[1] = false;
  nBinsX[1] = 17;
  binsTunedX[1]=bins_met;

  vector<vector<Float_t>> x_bins, y_bins;
  vector<Float_t> theBinsX, theBinsY;
  Float_t  binval=0;
  //
  for(UInt_t iV=0 ; iV<nV ; iV++) {

    // x bins
    theBinsX.clear();
    for(UInt_t iB=0 ; iB<nBinsX[iV] ; iB++) {
      if(regularX[iV]) {
	if(nBinsX[iV]!=0) binval = xFirst[iV] + (iB*(xLast[iV]-xFirst[iV])/nBinsX[iV]);
	else             binval = 0;
	theBinsX.push_back(binval);
      }
      else {
	theBinsX.push_back(binsTunedX[iV][iB]);
      }
    }
    x_bins.push_back(theBinsX);

    // y bins
    theBinsY.clear();
    for(UInt_t iB=0 ; iB<nBinsY[iV] ; iB++) {
      if(regularY[iV]) {
	if(nBinsY[iV]!=0) binval = yFirst[iV] + (iB*(yLast[iV]-yFirst[iV])/nBinsY[iV]);
	else              binval = 0;
	theBinsY.push_back(binval);
      }
      else {
	theBinsY.push_back(binsTunedY[iV][iB]);
      }
    }
    y_bins.push_back(theBinsY);
  }

  TString name,title;

  // Get process names from chains map
  vector<TString> process;
  M_PROC_CH::iterator itCh;
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {
    process.push_back( itCh->first );
  }
  const UInt_t nP = process.size();

  // Create histos and fill the map
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      for(UInt_t iC=0 ; iC<nC ; iC++) {

	for(UInt_t iV=0 ; iV<nV ; iV++) {

	  name ="h2D_"+selection[iS]+"_"+process[iP]+"_"+extra[iC]+"_"+var[iV];
	  title=selection[iS]+" "+process[iP]+" "+extra[iC]+" "+var[iV];

	  const UInt_t nBins = x_bins[iV].size();
	  Float_t theBins[nBins];
	  for(UInt_t iB=0 ; iB<nBins ; iB++) {
	    theBins[iB] = x_bins[iV][iB];
	  }

	  const UInt_t nBinsY = y_bins[iV].size(); 
	  Float_t theBinsY[nBinsY];
	  for(UInt_t iB=0 ; iB<nBinsY ; iB++) {
	    theBinsY[iB] = y_bins[iV][iB];
	  }

	  histos2D[selection[iS]][process[iP]][extra[iC]][var[iV]] = 
	    new TH2F(name, title, nBins-1, theBins, nBinsY-1, theBinsY);
	}
      }
    }
  }

  return 0;
}

Int_t defineChains(M_PROC_CH &chains, TString region1D, TString region2D)
{

  cout << "- START defineChains(chains, " << region1D << ", " << region2D << ");" << endl;

  TString pathMC="/user/ndaci/Data/XMET/AdishTrees_29Oct2015/Spring15MC_25ns/";
  TString nameFileMC="skim.root";
  TString nameTreeMC= "tree/tree";

  TString pathData= "/user/ndaci/Data/XMET/AdishTrees_29Oct2015/Run2015D/";
  TString pathDataHT= "/user/ndaci/Data/XMET/NadirTrees_26Oct2015/Run2015D/";
  TString nameFileData="tree.root";
  TString nameTreeData= "tree/tree";

  const UInt_t nMC=28;
  TString dirMC[nMC]={"znn100to200","znn200to400","znn400to600","znn600toinf",
		      "wln100to200","wln200to400","wln400to600","wln600toinf",
		      "zll100to200","zll200to400","zll400to600","zll600toinf",
		      "ttbar","ww","wz","zz",
		      "singletbart","singletbarw","singlett","singletw",
		      "qcdht1000to1500","qcdht100to200","qcdht1500to2000","qcdht2000toinf",
		      "qcdht200to300","qcdht300to500","qcdht500to700","qcdht700to1000"};

  // const UInt_t nData=2;
  // TString nameData[nData]={"data_jetht","data_met"};
  // TString dirData[nData] ={"jetht_05Oct2015","met"};

  // Define chains
  cout << "- loop over MC: ";
  for(UInt_t i=0 ; i<nMC ; i++) {
    cout << dirMC[i] << ":";
    chains[dirMC[i]] = new TChain(nameTreeMC);
    chains[dirMC[i]]->Add(pathMC+"/"+dirMC[i]+"/"+nameFileMC);
    cout << "done! ";
  }
  cout << endl;
  //
  if(region1D=="HT" || region2D=="HT") {
    cout << "- data_jetht:";
    chains["data_jetht"] = new TChain(nameTreeData);
    chains["data_jetht"]->Add(pathDataHT+"/jetht_05Oct2015/skimJSON_t1mumet50.root");
    cout << "done! ";
  }
  cout << endl;

  if(region1D=="MET" || region2D=="MET") {
    cout << "- datamet:";
    chains["data_met"] = new TChain(nameTreeData);
    chains["data_met"]->Add(pathData+"/met/tree.root");
    cout << "done! ";
  }
  cout << endl;

  cout << "- END defineChains" << endl;
  return 0;
}

TCut defineCut(TString sample, TString region)
{

  //TCut trig ="hltmet90>0 || hltmet120>0";
  TCut trig="";
  TCut noise="flaghbheloose>0 && flagcsctight>0 && flageebadsc>0";
  TCut maxrun=""; 
  TCut metCut="";
  TCut leptons="nelectrons == 0 && ntaus == 0 && nmuons == 0";
  TCut photons="nphotons==0";
  TCut bveto="nbjets==0";
  TCut jetID="signaljetpt>100 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && signaljetNHfrac < 0.8";
  TCut jetBin="njets>=1";
  TCut noqcd="";
  TCut fwdveto="";

  // QCD Killer
  if(region.Contains("JetMet0p5")) noqcd="abs(jetmetdphimin)>0.5";

  // FwdVeto
  if(region.Contains("FwdVeto")) {
    fwdveto="abs(leadingjeteta)<2.5"; // Adish trees
    if(sample.Contains("jetht")) fwdveto="(signaljetpt>Max$(jet_pt * (abs(jet_eta)>2.5)))"; // Nadir trees
  }

  // HT or MET region
  if(region.Contains("HT")) {
    if(sample.Contains("jetht")) trig = "hltpfht200 && ht>300"; // Nadir trees
    else                         trig = "ht>300"; // no hltpfht branch in Adish trees
  }  

  else if(region.Contains("MET")) {
    //trig = "hltmet90>0 || hltmet120>0";
    trig = "hltmet90>0"; // to be sure of the lumi measurement
    if(region.Contains("Met200")) metCut="t1mumet>200";
    maxrun="run<257599"; // stay blind only in MET region
    noise *= "flaghbheiso>0"; // fixme: add to everybody once I get updated JetHT trees
  }

  // Jet ID //
  TCut jetID1 = "(signaljetpt>100 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && signaljetNHfrac < 0.8)";
  TCut jetID2 = "(secondjetpt>30  && abs(secondjeteta)<2.5 && secondjetCHfrac > 0.1 && secondjetNHfrac < 0.8)";
  TCut jetID3 = "( thirdjetpt>30  && abs(thirdjeteta)<2.5  &&  thirdjetCHfrac > 0.1 &&  thirdjetNHfrac < 0.8)";

  if(     region.Contains("1jet"))   {jetBin="njets==1"; jetID = jetID1;}
  else if(region.Contains("2jet"))   {jetBin="njets==2"; jetID = jetID1*jetID2;}
  else if(region.Contains("3jet"))   {jetBin="njets>=3"; jetID = jetID1*jetID2*jetID3;}
  //else if(region.Contains("4jet"))   {jetBin="njets>=4"; jetID = jetID1*jetID2*jetID3;}

  return trig*noise*maxrun*metCut*leptons*photons*bveto*jetID*jetBin*noqcd*fwdveto;
  //return "t1mumet>200"; // FIXME: debug test

}

vector<TString> GetVarName(M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D, TString mode)
{

  vector<TString> var;

  M_SEL_PROC_CUT_VAR_H1D::iterator itSelProcCutVarH1D;
  M_PROC_CUT_VAR_H1D::iterator     itProcCutVarH1D;
  M_CUT_VAR_H1D::iterator          itCutVarH1D;
  M_VAR_H1D::iterator              itVarH1D;

  M_SEL_PROC_CUT_VAR_H2D::iterator itSelProcCutVarH2D;
  M_PROC_CUT_VAR_H2D::iterator     itProcCutVarH2D;
  M_CUT_VAR_H2D::iterator          itCutVarH2D;
  M_VAR_H2D::iterator              itVarH2D;

  if(mode=="1D" || mode=="both") {
    itSelProcCutVarH1D=histos1D.begin();
    itProcCutVarH1D=(itSelProcCutVarH1D->second).begin();
    itCutVarH1D=(itProcCutVarH1D->second).begin();
  
    for((itVarH1D=(itCutVarH1D->second).begin()) ; 
	(itVarH1D!=(itCutVarH1D->second).end()) ; 
	itVarH1D++) {
      var.push_back(itVarH1D->first);
    }
  }

  if(mode=="2D" || mode=="both") {

    itSelProcCutVarH2D=histos2D.begin();
    itProcCutVarH2D=(itSelProcCutVarH2D->second).begin();
    itCutVarH2D=(itProcCutVarH2D->second).begin();
  
    for((itVarH2D=(itCutVarH2D->second).begin()) ; 
	(itVarH2D!=(itCutVarH2D->second).end()) ; 
	itVarH2D++) {
      var.push_back(itVarH2D->first);
    }
  }
  
  return var;

}

Int_t GetNGen(TString sample)
{
  
  if(     sample=="qcdht1000to1500") return 4963895;
  else if(sample=="qcdht100to200")   return 81637494;
  else if(sample=="qcdht1500to2000") return 3868886;
  else if(sample=="qcdht2000toinf")  return 1912529;
  else if(sample=="qcdht200to300")   return 18718905;
  else if(sample=="qcdht300to500")   return 19826197;
  else if(sample=="qcdht500to700")   return 19664159;
  else if(sample=="qcdht700to1000")  return 15356448;
  else if(sample=="singletbarw")     return 988500;
  else if(sample=="singletw")        return 995600;
  else if(sample=="singletbart")     return 1680200;
  else if(sample=="singlett")        return 3299800;
  else if(sample=="ttbar")           return 11344206;
  else if(sample=="wln100to200")     return 10152718;
  else if(sample=="wln200to400")     return 5221599;
  else if(sample=="wln400to600")     return 1745914;
  else if(sample=="wln600toinf")     return 1039152;
  else if(sample=="ww")              return 993640;
  else if(sample=="wz")              return 978512;
  else if(sample=="zll100to200")     return 2725655;
  else if(sample=="zll200to400")     return 973937;
  else if(sample=="zll400to600")     return 1067758;
  else if(sample=="zll600toinf")     return 998912;
  else if(sample=="znn100to200")     return 5154824;
  else if(sample=="znn200to400")     return 4998316;
  else if(sample=="znn400to600")     return 1018882;
  else if(sample=="znn600toinf")     return 1008333;
  else if(sample=="zz")              return 996944;

  else return 1;
}


TGraphErrors TransferFactor(TH2F *hTemp2, Float_t cut)
{
  
  UInt_t nBinsX = hTemp2->GetNbinsX();
  UInt_t nBinsY = hTemp2->GetNbinsY();

  Float_t low1=0;
  Float_t low2=0;
  UInt_t  idxB=0;
  
  for(UInt_t iB=1 ; iB<nBinsY ; iB++) {
    low1 = hTemp2->GetYaxis()->GetBinLowEdge(iB);
    low2 = hTemp2->GetYaxis()->GetBinLowEdge(iB+1);
    if(low1<cut && low2>=cut) idxB=iB;
  }

  Double_t integral1, integral2, error1, error2, ratio, error;
  integral1 = integral2 = error1 = error2 = ratio = error = 0;

  Double_t posX[nBinsX];
  Double_t posY[nBinsX];
  Double_t errX[nBinsX];
  Double_t errY[nBinsX];

  for(UInt_t iB=1 ; iB<nBinsX+1 ; iB++) {

    integral1 = hTemp2->IntegralAndError(iB, iB, 0, idxB, error1);
    integral2 = hTemp2->IntegralAndError(iB, iB, idxB+1, nBinsY+1, error2);
    ratio = integral1!=0 ? integral2/integral1 : -0.1;
    error = ErrorRatio(integral2, integral1, error2, error1);
    //error = QuadSum(error1, error2);
    //hSF->SetBinContent(iB, ratio);
    //hSF->SetBinError(  iB, error);

    posX[iB-1] = hTemp2->GetXaxis()->GetBinCenter(iB);
    errX[iB-1] = abs(posX[iB-1] - hTemp2->GetXaxis()->GetBinLowEdge(iB));
    posY[iB-1] = ratio;
    errY[iB-1] = error;

    cout << "MET=" << posX[iB-1] << "("   << integral2 << "+/-" << error2 << "/" << integral1 << "+/-" << error1 
	 << "="    << ratio      << "+/-" << error     << endl;
  }

  //gSF = new TGraphErrors( nBinsX, posX, posY, errX, errY );
  TGraphErrors gSF( nBinsX, posX, posY, errX, errY );

  return gSF;
}
