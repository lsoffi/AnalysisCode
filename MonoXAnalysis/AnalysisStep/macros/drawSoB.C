#include "myIncludes.h"

vector<Int_t> findIdx(TString name);

Int_t drawSoB(TString prefix="v3_AN15_JetMet_", TString postfix="", TString var="t1mumet")
{

  ofstream outlog("plots/soverb/"+prefix+"_"+postfix+"/eff_"+prefix+"_"+postfix+".txt");

  outlog << setw(9) << "MET"
	 << setw(7) << "#Jets"
	 << setw(5) << "Cut";

  /*
  for(UInt_t iSrc=0 ; iSrc<nSrc ; iSrc++) {
    for(UInt_t iSel=0 ; iSel<nSel ; iSel++) {
      for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
	for(UInt_t iS=0 ; iS<nSig ; iS++) {
	  for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {

  */


  // GET INPUTS //
  const UInt_t nSrc=3;
  TString metRange[nSrc] = {"MetFrom200to250","MetFrom250to350","Met350"};
  TString shortRange[nSrc]={"200-250","250-350",">350"};

  TFile* file[nSrc];
  TString tag[nSrc];
  TString pathFile;
  for(UInt_t iSrc=0 ; iSrc<nSrc; iSrc++) {
    tag[iSrc]  = prefix+metRange[iSrc]+postfix;
    pathFile = "plots/"+tag[iSrc]+"/plots_"+tag[iSrc]+".root";
    cout << "- Open file: " << pathFile << " ... " ;
    file[iSrc] = new TFile(pathFile,"read");
    cout << "file opened!" << endl;
  }

  // Processes to use
  cout << "- Define processes" << endl;
  vector<TString> signal, background;

  // Backgrounds: MC + W DD
  //background.push_back("zll"); 
  //background.push_back("znn"); 
  //background.push_back("wjets"); 
  //background.push_back("ttbar"); 
  //background.push_back("top"); 
  //background.push_back("vv"); 
  background.push_back("qcd"); 
  //
  // W DD
  //background.push_back("wj_sr_DD");      // FIXME
  //background.push_back("wj_sr_corr_DD"); // FIXME

  // Signal + Irreducible Bkg
  //signal.push_back("qcd"); 
  signal.push_back("znn"); 
  //signal.push_back("wjets");
  //signal.push_back("zn_sr_DD"); // FIXME
  //signal.push_back("zn_sr_corr_DD"); // FIXME
  //
  // Signal
  /* // FIXME
  const UInt_t nSpin=2;
  const UInt_t nMass=9;
  TString nameSpin[nSpin] = {"V","AV"};
  TString nameMass[nMass] = {"0p1","1","10","100","200","300","400","700","1000"};
  TString nameProcess;
  for(UInt_t iS=0 ; iS<nSpin ; iS++) {
    for(UInt_t iM=0 ; iM<nMass ; iM++) {
      nameProcess = "dm_"+nameSpin[iS]+"_"+nameMass[iM];
      signal.push_back(nameProcess);
    }
  }
  */

  const UInt_t nSig = signal.size();
  const UInt_t nBkgIn = background.size();

  // Selections and variables
  const UInt_t nSel=4;
  TString sel[nSel] = {"1jet","2jet","3jet","4jet"};
  TString shortsel[nSel] = {"1","2","3","4"};

  const UInt_t nCut=5;
  TString cut[  nCut] = {"NoJmCut","JetMet0p2","JetMet0p4","JetMet0p6","JetMet0p8"};
  TString shortcut[nCut]={"No","0.2","0.4","0.6","0.8"};

  // Loop over input files to get inputs
  TString nameH;
  TH1F*   hTemp;
  Double_t integral;

  const UInt_t nSB=2; //0:Signal 1:Background
  vector<Double_t> vIntegral[nCut][nSrc][nSel][nSB];

  // Define S/B treatment
  //const UInt_t nSig=3+(nSpin*nMass); // ZNN:MC,DD,DDcorr;DM
  //const UInt_t nBkgOut=4; // QCD, SumBkg, SumBkg(W=DD), SumBkg(W=DDcorr) // FIXME
  //TString bkgout[nBkgOut] = {"QCD","SumBkg","SumBkgWDD","SumBkgWDDc"}; // FIXME

  //const UInt_t nBkgOut=2; // QCD, SumBkg
  //TString bkgout[nBkgOut] = {"QCD","SumBkg"};
  const UInt_t nBkgOut=1; // QCD
  //TString bkgout[nBkgOut] = {"QCD"};
  TString bkgout[nBkgOut] = {"SumBkg"};
  vector<Int_t> idxBkg;
  Double_t theS, theB, theSoB, totS, totB;
  theS = theB = theSoB = totS = totB = 0;

  Double_t effS[nSig], effB[nBkgOut];

  for(UInt_t iS=0 ; iS<nSig ; iS++) {\
    outlog << setw(12) << signal[iS];
  }
 
  for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
    outlog << setw(12) << bkgout[iB];
  }

  outlog << endl;
 
  // Define S/B histograms
  cout << "- Define S/B histograms" << endl;
  TString nameSoB, titleSoB;
  TH2F* hSoB[nCut][nSig][nBkgOut];
  //
  // bins
  const UInt_t nBinsMet = 4;
  const UInt_t nBinsJet = 5;
  Float_t metBins[nBinsMet] = {200, 250, 350, 450};
  Float_t jetBins[nBinsJet] = {0.5, 1.5, 2.5, 3.5, 4.5};
  //
  Float_t jetBinX[nSel]  = {1, 2, 3, 4};
  Float_t metBinY[nSrc]  = {225, 300, 400};
  // 
  // Loop over histograms to declare/define them
  for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
    for(UInt_t iS=0 ; iS<nSig ; iS++) {
      for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	nameSoB  = "h2_"+cut[iCut]+"_"+signal[iS]+"_"+bkgout[iB] ;
	titleSoB = nameSoB ;
	cout << "---- name : " << nameSoB << endl;
	hSoB[iCut][iS][iB] = 
	  new TH2F(nameSoB, titleSoB, nBinsJet-1, jetBins, nBinsMet-1, metBins);
	hSoB[iCut][iS][iB]->SetXTitle("Jet multiplicity");
	hSoB[iCut][iS][iB]->SetYTitle("PFMETNoMu [GeV]");
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;

  // Loop over input histograms to compute S/B
  cout << "- Loop over input histograms" << endl;
  for(UInt_t iSrc=0 ; iSrc<nSrc ; iSrc++) {
    for(UInt_t iSel=0 ; iSel<nSel ; iSel++) {
      for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {

	cout << "---- Case: " 
	     << metRange[iSrc] << " " 
	     << sel[iSel] << " " 
	     << cut[iCut] << endl;

	outlog << setw(9) << shortRange[iSrc]
	       << setw(7) << shortsel[iSel]
	       << setw(5) << shortcut[iCut];


	// Signal
	cout << "---- Look at signals: " ;
	for(UInt_t iS=0 ; iS<nSig ; iS++) {
	  nameH = "h_"+var+"_"+signal[iS]+"_"+sel[iSel]+"_"+cut[iCut];
	  hTemp = (TH1F*) file[iSrc]->Get(nameH);
	  integral = hTemp->Integral(0, hTemp->GetNbinsX() + 1);
	  vIntegral[iCut][iSrc][iSel][0].push_back(integral);
	  cout << "#Int(" << nameH << ")=" << integral << " | ";
	  effS[iS]=0;
	}// end loop:signal
	cout << endl << "---- ...done!" << endl;

	// Background: initialize output sums
	for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	  vIntegral[iCut][iSrc][iSel][1].push_back(0);
	  effB[iB]=0;
	}

	// Background: Loop over input bkg
	cout << "---- Look at backgrounds: " ;
	for(UInt_t iB=0 ; iB<nBkgIn ; iB++) {
	  nameH = "h_"+var+"_"+background[iB]+"_"+sel[iSel]+"_"+cut[iCut];
	  hTemp = (TH1F*) file[iSrc]->Get(nameH);
	  integral = hTemp->Integral(0, hTemp->GetNbinsX() + 1);
	  cout << "#Int(" << nameH << ")=" << integral << " | ";

	  idxBkg = findIdx(background[iB]);
	  for(UInt_t i=0 ; i<idxBkg.size() ; i++) {
	    if(idxBkg[i]<0 || idxBkg[i]>=(Int_t)nBkgOut) continue;
	    vIntegral[iCut][iSrc][iSel][1][idxBkg[i]] += integral;
	  }
	}// end loop:backgrounds
	cout << endl << "---- ...done!" << endl;
	
	// Fill histogram
	cout << "---- Fill histogram: " << endl;

	for(UInt_t iS=0 ; iS<nSig ; iS++) {

	  // Compute signal informations
	  theS = vIntegral[iCut][iSrc][iSel][0][iS];
	  if(iCut==0) totS = theS;
	  effS[iS] = totS!=0 ? theS / totS : -1;

	  for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {

	    // compute background informations
	    theB = vIntegral[iCut][iSrc][iSel][1][iB];
	    if(iCut==0) totB = theB;
	    effB[iB] = totB!=0 ? theB / totB : -1;

	    // compute S/B
	    theSoB = theB!=0 ? theS / theB : -999;
	    hSoB[iCut][iS][iB]->Fill( jetBinX[iSel] , metBinY[iSrc] , theSoB );
	    cout << hSoB[iCut][iS][iB]->GetName() 
		 << "->Fill(" << jetBinX[iSel] 
		 << ","       << metBinY[iSrc] 
		 << ","       << theSoB 
		 << ") | " ;
	  }
	}

	// put efficiencies in the outlog
	for(UInt_t iS=0 ; iS<nSig ; iS++) {
	  outlog << setw(12) << effS[iS];
	}
	for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	  outlog << setw(12) << effB[iB];
	}
	outlog << endl;

	cout << endl << "---- ...done!" << endl;

      }// end loop:scanned cuts
      cout << endl;
    }// end loop:selections
    cout << endl;
  }// end loop:files
  cout << endl;


  // DRAW HISTOGRAMS //
  TFile *outfile = new TFile("plots/soverb/"+prefix+"_"+postfix+"/sob_"+prefix+"_"+postfix+".root","recreate");
  outfile->cd();
  //
  TCanvas *c[nCut][nSig][nBkgOut];
  TString nameC;
  gStyle->SetOptStat(0);
  //
  for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
    for(UInt_t iS=0 ; iS<nSig ; iS++) {
      for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	outfile->cd();
	hSoB[iCut][iS][iB]->Write();

	nameC = "c_"+signal[iS]+"_"+bkgout[iB]+"_"+cut[iCut];
	c[iCut][iS][iB] = new TCanvas(nameC,nameC,20,20,600,600);
	setStyle(c[iCut][iS][iB]);
	c[iCut][iS][iB]->cd();

	//hSoB[iCut][iS][iB]->Draw("colztexte");
	hSoB[iCut][iS][iB]->Draw("colztext");
      }
    }
  }

  // PRINT PLOTS IN PDF FILES //

  // All plots in a single pdf file
  for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
    for(UInt_t iS=0 ; iS<nSig ; iS++) {
      for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	TString pdftitle = "Title:"+cut[iCut]+" "+signal[iS]+" "+bkgout[iB];
	TString pdfp = "";
	if(iCut==0 && iS==0 && iB==0) pdfp = "(";
	else if(iCut==nCut-1 && iS==nSig-1 && iB==nBkgOut-1) pdfp = ")";
	c[iCut][iS][iB]->Print("plots/soverb/"+prefix+"_"+postfix+"/sob_"+prefix+"_"+postfix+".pdf"+pdfp, pdftitle);
      }
    }
  }

  // 1 pdf file per cut
  for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
    for(UInt_t iS=0 ; iS<nSig ; iS++) {
      for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	TString pdftitle = "Title:"+cut[iCut]+" "+signal[iS]+" "+bkgout[iB];
	TString pdfp = "";
	if(iS==0 && iB==0) pdfp = "(";
	else if(iS==nSig-1 && iB==nBkgOut-1) pdfp = ")";
	c[iCut][iS][iB]->Print("plots/soverb/"+prefix+"_"+postfix+"/sob_"+prefix+"_"+cut[iCut]+"_"+postfix+".pdf"+pdfp, pdftitle);
      }
    }
  }

  // 1 pdf file per background
  for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
    for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
      for(UInt_t iS=0 ; iS<nSig ; iS++) {
	TString pdftitle = "Title:"+cut[iCut]+" "+signal[iS]+" "+bkgout[iB];
	TString pdfp = "";
	if(iCut==0 && iS==0) pdfp = "(";
	else if(iCut==nCut-1 && iS==nSig-1) pdfp = ")";
	c[iCut][iS][iB]->Print("plots/soverb/"+prefix+"_"+postfix+"/sob_"+prefix+"_"+bkgout[iB]+"_"+postfix+".pdf"+pdfp, pdftitle);
      }
    }
  }

  // 1 pdf file per signal
  for(UInt_t iS=0 ; iS<nSig ; iS++) {
    for(UInt_t iCut=0 ; iCut<nCut ; iCut++) {
      for(UInt_t iB=0 ; iB<nBkgOut ; iB++) {
	TString pdftitle = "Title:"+cut[iCut]+" "+signal[iS]+" "+bkgout[iB];
	TString pdfp = "";
	if(iCut==0 && iB==0) pdfp = "(";
	else if(iCut==nCut-1 && iB==nBkgOut-1) pdfp = ")";
	c[iCut][iS][iB]->Print("plots/soverb/"+prefix+"_"+postfix+"/sob_"+prefix+"_"+signal[iS]+"_"+postfix+".pdf"+pdfp, pdftitle);
      }
    }
  }


  // END //
  cout << "THE END" << endl;

  return 0;
}

vector<Int_t> findIdx(TString name)
{

  vector<Int_t> res;

  if(name.Contains("DD")) {
    if(name.Contains("corr")) {
      res.push_back(3);
    }
    else {
      res.push_back(2);
    }
  }
  else if(name=="qcd") {
    res.push_back(0);
    res.push_back(1);
    res.push_back(2);
    res.push_back(3);
  }
  else if(name=="SumBkg") {
    res.push_back(0);
  }
  else {
    res.push_back(1);
    res.push_back(2);
    res.push_back(3);
  }

  return res;
}
