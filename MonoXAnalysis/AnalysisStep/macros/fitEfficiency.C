#include "myIncludes.h"

bool DEBUG = false;

//           run        lumi      events counter
typedef map< int , map< int , map< int , int> > > MAP_RLE;

struct STEP{
  double  T; // threshold
  TString c; // collection
  TString n; // name
  TString t; // title

  Int_t   C; // color
  Int_t   S; // style
  
  double pt;
  double phi;
  bool   pass;
  bool   serial;
};

struct PATH{
  TString nameP;
  TString namePath;
  UInt_t nSteps;
  vector<STEP> steps;
};

pair<Int_t, Int_t> getStyle(TString name);
Double_t evaluate(double *x, double *par);
Double_t evaluate2(double *x, double *par);
Double_t ApproxErf(Double_t arg);
Double_t dichotomy(double eff, double a0, double b0, double relErr,
		   TF1 f, bool verbose);

Int_t fitEfficiency(TString resultName="v0_test",
		    TString seed="ETM50",
		    TString HBHECleaning="NoHBHE")
{

  TFile* infile = new TFile("results/"+resultName+"/f_"+resultName+".root","read");
  if(infile->IsZombie()) {
    cout << "ERROR: input file is a zombie. Exit ==> []" << endl;
    return -1;
  }

  ///////////////////////////////////
  // Define trigger interpretation //
  ///////////////////////////////////

  const UInt_t nV=6; // mumet, t1mumet, pfmet, t1pfmet, signaljetpt, signaljetNHfrac
  const UInt_t nF=2; // denom, num
  const UInt_t nP=11; // PFMNoMu90,120; Jet80PFMNoMu90,120; PFM90,120; PFM170; CM200; PFMNoMu90_Up, Jet80PFMNoMu90_Up; OR90GeV
  UInt_t nS;

  cout << "Define steps" << endl;

  TString theStep;

  // Level-1 Seeds
  STEP s_L1ETM60 ={.T=60,.c="hltL1extraParticles:MET:HLT",.n="L1",.t="L1",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};               
  STEP s_L1ETM50 ={.T=50,.c="hltL1extraParticles:MET:HLT",.n="L1",.t="L1",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};               

  // MonoPFJet
  STEP s_Jet65   ={.T=65,.c="hltAK4CaloJetsCorrectedIDPassed::HLT",.n="Jet",  .t="Jet",  .C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_PFJet80 ={.T=80,.c="hltAK4PFJetsTightIDCorrected::HLT",   .n="PFJet",.t="PFJet",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_bJ80MuPFM90 ={.T=0,.c="hltjetmet90",                   .n="Full",.t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_bJ80MuPFM120={.T=0,.c="hltjetmet120",                  .n="Full",.t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // OR
  STEP s_bOR90GeV = {.T=0,.c="OR90GeV",                  .n="Full",.t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMNoMu90
  STEP s_MET65   ={.T=65,.c="hltMet::HLT",               .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  STEP s_METC55  ={.T=55,.c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};      
  STEP s_METJ55  ={.T=55,.c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_MHT65   ={.T=65,.c="hltMht::HLT",               .n="MHT",  .t="MHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_MuMHT90 ={.T=90,.c="hltPFMHTNoMuTightID::HLT",  .n="PFMHT",.t="PFMHTNoMu",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_MuMET90 ={.T=90,.c="hltPFMETNoMuProducer::HLT", .n="PFMET",.t="PFMETNoMu",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_bMuPFM90={.T=0, .c="hltmet90",                  .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMNoMu120
  STEP s_MET80    ={.T=80, .c="hltMet::HLT",               .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METC70   ={.T=70, .c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METJ70   ={.T=70, .c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_MHT80    ={.T=80, .c="hltMht::HLT",               .n="MHT",  .t="MHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_MuMHT120 ={.T=120,.c="hltPFMHTNoMuTightID::HLT",  .n="PFMHT",.t="PFMHTNoMu",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_MuMET120 ={.T=120,.c="hltPFMETNoMuProducer::HLT", .n="PFMET",.t="PFMETNoMu",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_bMuPFM120={.T=0,  .c="hltmet120",                 .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFM90
  STEP s_MET70   ={.T=70,.c="hltMet::HLT",          .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  STEP s_MHT70   ={.T=70,.c="hltMht::HLT",          .n="MHT",  .t="MHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_PFMHT90 ={.T=90,.c="hltPFMHTTightID::HLT", .n="PFMHT",.t="PFMHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_PFMET90 ={.T=90,.c="hltPFMETProducer::HLT",.n="PFMET",.t="PFMET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_bPFM90  ={.T=0, .c="hltmetwithmu90",       .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFM120
  STEP s_MET90   ={.T=90,.c="hltMet::HLT",            .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};              
  STEP s_MHT90   ={.T=90,.c="hltMht::HLT",            .n="MHT",  .t="MHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_PFMHT120={.T=120,.c="hltPFMHTTightID::HLT", .n="PFMHT",.t="PFMHT",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_PFMET120={.T=120,.c="hltPFMETProducer::HLT",.n="PFMET",.t="PFMET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0}; 
  STEP s_bPFM120 ={.T=0, .c="hltmetwithmu120",       .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PFMET170
  //STEP s_MET90   ={.T=90, .c="hltMet::HLT",               .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METC80  ={.T=80, .c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METJ80  ={.T=80, .c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_PFMET170={.T=170,.c="hltPFMETProducer::HLT",     .n="PFMET",.t="PFMET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_bMET170 ={.T=0,  .c="hltmetwithmu170",           .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // CaloMET200
  STEP s_MET210  ={.T=210, .c="hltMet::HLT",               .n="MET",  .t="MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METC200 ={.T=200, .c="hltMetClean::HLT",          .n="METC", .t="Cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  STEP s_METJ200 ={.T=200, .c="hltMetCleanUsingJetID::HLT",.n="METJ", .t="JetID-cleaned MET",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};
  //STEP s_bMET170 ={.T=0,   .c="hltmetwithmu170",           .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  // PATHS //
  Bool_t useHBHECleaning = (HBHECleaning!="NoHBHE");
  cout << "Define paths" << endl;
  PATH myPaths[nP];
  vector<STEP> vStepEmpty;
  
  myPaths[0]={.nameP="PFMNoMu90",.namePath="HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[0].steps.clear();
  //
  if(     seed=="ETM50") myPaths[0].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[0].steps.push_back(s_L1ETM60);
  //
  myPaths[0].steps.push_back(s_MET65);
  if(useHBHECleaning) myPaths[0].steps.push_back(s_METC55);
  myPaths[0].steps.push_back(s_METJ55);
  myPaths[0].steps.push_back(s_MHT65);
  myPaths[0].steps.push_back(s_MuMHT90);
  myPaths[0].steps.push_back(s_MuMET90);
  myPaths[0].steps.push_back(s_bMuPFM90);
  myPaths[0].nSteps = myPaths[0].steps.size();

  myPaths[1]={.nameP="PFMNoMu120",.namePath="HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[1].steps.clear();
  //
  if(     seed=="ETM50") myPaths[1].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[1].steps.push_back(s_L1ETM60);
  //
  myPaths[1].steps.push_back(s_MET80);
  if(useHBHECleaning) myPaths[1].steps.push_back(s_METC70);
  myPaths[1].steps.push_back(s_METJ70);
  myPaths[1].steps.push_back(s_MHT80);
  myPaths[1].steps.push_back(s_MuMHT120);
  myPaths[1].steps.push_back(s_MuMET120);
  myPaths[1].steps.push_back(s_bMuPFM120);
  myPaths[1].nSteps = myPaths[1].steps.size();

  myPaths[2]={.nameP="Jet80PFMNoMu90",.namePath="HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[2].steps.clear();
  //
  if(     seed=="ETM50") myPaths[2].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[2].steps.push_back(s_L1ETM60);
  //
  myPaths[2].steps.push_back(s_MET65);
  myPaths[2].steps.push_back(s_Jet65);
  if(useHBHECleaning) myPaths[2].steps.push_back(s_METC55);
  myPaths[2].steps.push_back(s_METJ55);
  myPaths[2].steps.push_back(s_MHT65);
  myPaths[2].steps.push_back(s_PFJet80);
  myPaths[2].steps.push_back(s_MuMHT90);
  myPaths[2].steps.push_back(s_MuMET90);
  myPaths[2].steps.push_back(s_bJ80MuPFM90);
  myPaths[2].nSteps = myPaths[2].steps.size();

  myPaths[3]={.nameP="Jet80PFMNoMu120",.namePath="HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[3].steps.clear();
  //
  if(     seed=="ETM50") myPaths[3].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[3].steps.push_back(s_L1ETM60);
  //
  myPaths[3].steps.push_back(s_MET80);
  myPaths[3].steps.push_back(s_Jet65);
  if(useHBHECleaning) myPaths[3].steps.push_back(s_METC70);
  myPaths[3].steps.push_back(s_METJ70);
  myPaths[3].steps.push_back(s_MHT80);
  myPaths[3].steps.push_back(s_PFJet80);
  myPaths[3].steps.push_back(s_MuMHT120);
  myPaths[3].steps.push_back(s_MuMET120);
  myPaths[3].steps.push_back(s_bJ80MuPFM120);
  myPaths[3].nSteps = myPaths[3].steps.size();

  myPaths[4]={.nameP="PFM90",.namePath="HLT_PFMET90_PFMHT90_IDTight",.nSteps=6,.steps=vStepEmpty};
  myPaths[4].steps.clear();
  //
  if(     seed=="ETM50") myPaths[4].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[4].steps.push_back(s_L1ETM60);
  //
  myPaths[4].steps.push_back(s_MET70);
  myPaths[4].steps.push_back(s_MHT70);
  myPaths[4].steps.push_back(s_PFMHT90);
  myPaths[4].steps.push_back(s_PFMET90);
  myPaths[4].steps.push_back(s_bPFM90);
  myPaths[4].nSteps = myPaths[4].steps.size();

  myPaths[5]={.nameP="PFM120",.namePath="HLT_PFMET120_PFMHT120_IDTight",.nSteps=6,.steps=vStepEmpty};
  myPaths[5].steps.clear();
  //
  if(     seed=="ETM50") myPaths[5].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[5].steps.push_back(s_L1ETM60);
  //
  myPaths[5].steps.push_back(s_MET90);
  myPaths[5].steps.push_back(s_MHT90);
  myPaths[5].steps.push_back(s_PFMHT120);
  myPaths[5].steps.push_back(s_PFMET120);
  myPaths[5].steps.push_back(s_bPFM120);
  myPaths[5].nSteps = myPaths[5].steps.size();

  myPaths[6]={.nameP="PFMET170",.namePath="HLT_PFMET170_JetIdCleaned",.nSteps=6,.steps=vStepEmpty};
  myPaths[6].steps.clear();
  //
  if(     seed=="ETM50") myPaths[6].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[6].steps.push_back(s_L1ETM60);
  //
  myPaths[6].steps.push_back(s_MET90);
  if(useHBHECleaning) myPaths[6].steps.push_back(s_METC80);
  myPaths[6].steps.push_back(s_METJ80);
  myPaths[6].steps.push_back(s_PFMET170);
  myPaths[6].steps.push_back(s_bMET170);
  myPaths[6].nSteps = myPaths[6].steps.size();

  myPaths[7]={.nameP="CaloMET200",.namePath="HLT_CaloMET200_JetIdCleaned",.nSteps=3,.steps=vStepEmpty};
  myPaths[7].steps.clear();
  //
  if(     seed=="ETM50") myPaths[7].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[7].steps.push_back(s_L1ETM60);
  //
  myPaths[7].steps.push_back(s_MET210);
  if(useHBHECleaning) myPaths[7].steps.push_back(s_METC200);
  myPaths[7].steps.push_back(s_METJ200);
  myPaths[7].nSteps = myPaths[7].steps.size();

  myPaths[8]={.nameP="PFMNoMu90_Up",.namePath="HLT_CaloUp_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[8].steps.clear();
  //
  if(     seed=="ETM50") myPaths[8].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[8].steps.push_back(s_L1ETM60);
  //
  myPaths[8].steps.push_back(s_MET80);
  if(useHBHECleaning) myPaths[8].steps.push_back(s_METC70);
  myPaths[8].steps.push_back(s_METJ70);
  myPaths[8].steps.push_back(s_MHT80);
  myPaths[8].steps.push_back(s_MuMHT90);
  myPaths[8].steps.push_back(s_MuMET90);
  myPaths[8].steps.push_back(s_bMuPFM90);
  myPaths[8].nSteps = myPaths[8].steps.size();

  myPaths[9]={.nameP="Jet80PFMNoMu90_Up",.namePath="HLT_CaloUp_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[9].steps.clear();
  //
  if(     seed=="ETM50") myPaths[9].steps.push_back(s_L1ETM50);
  else if(seed=="ETM60") myPaths[9].steps.push_back(s_L1ETM60);
  //
  myPaths[9].steps.push_back(s_MET80);
  myPaths[9].steps.push_back(s_Jet65);
  if(useHBHECleaning) myPaths[9].steps.push_back(s_METC70);
  myPaths[9].steps.push_back(s_METJ70);
  myPaths[9].steps.push_back(s_MHT80);
  myPaths[9].steps.push_back(s_PFJet80);
  myPaths[9].steps.push_back(s_MuMHT90);
  myPaths[9].steps.push_back(s_MuMET90);
  myPaths[9].steps.push_back(s_bMuPFM90);
  myPaths[9].nSteps = myPaths[9].steps.size();

  myPaths[10]={.nameP="OR90GeV",.namePath="HLT_PFMNoMu90_PFMET170",.nSteps=1,.steps=vStepEmpty};
  myPaths[10].steps.clear();
  myPaths[10].steps.push_back(s_bOR90GeV);
  myPaths[10].nSteps = myPaths[10].steps.size();

  cout << "Set style" << endl;

  // Set style per step in all paths
  pair<Int_t, Int_t> theStyle = make_pair(0,0);

  for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
    nS = myPaths[iP].nSteps;
    for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
      theStyle = getStyle(myPaths[iP].steps[iS].n);
      myPaths[iP].steps[iS].C = theStyle.first;
      myPaths[iP].steps[iS].S = theStyle.second;
    }
  }
  // end set style

  ////////////////
  // HISTOGRAMS //
  ////////////////

  cout << "Define histograms." << endl;

  vector<TH1F*> h[nV][nF][nP];

  TString hname, title;
  TString nameV[nV]={"mumet","t1mumet","pfmet","t1pfmet","signaljetpt","signaljetNHfrac"};

  TString nameAxis[nV]={"Reco PFMETNoMu [GeV]",
			"Type1 PFMETNoMu [GeV]",
			"Reco PFMET [GeV]",
			"Type1 PFMET",
			"Leading PFJet p_{T}",
			"Leading PFJet NHEF"};
  /*
  int   xbins_reg[nV] = {40,  40,  40,  40,  40,  50};
  float xlow_reg[ nV] = {100, 100, 100, 100, 100, 0};
  float xup_reg[  nV] = {900, 900, 900, 900, 900, 1};

  const UInt_t xbins[nV]    = {21, 21, 21, 21, 21, 27};

  float bins_met[]  = {50,  75,  100, 110, 120, 
		       130, 140, 150, 160, 170, 
		       180, 190, 200, 220, 250, 
		       300, 350, 400, 500, 650, 
		       1000};

  float bins_nhef[] = {0.00, 0.02, 0.04, 0.06, 0.08, 
		       0.10, 0.12, 0.14, 0.16, 0.18,
		       0.20, 0.22, 0.24, 0.26, 0.28,
		       0.30, 0.35, 0.40, 0.45, 0.50,
		       0.55, 0.60, 0.65, 0.70, 0.80,
		       0.90, 1.00};
  
  float* v_xlow[nV] = {bins_met , bins_met , bins_met , bins_met , bins_met , bins_nhef};
  */

  TString nameF[nF]={"denom","num"};

  TH1F* hTemp;

  cout << "- Get histograms" << endl;

  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths

	nS = myPaths[iP].nSteps;
	cout << "---- path: " << myPaths[iP].nameP << " ; nS=" << nS << endl;

	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths

	  cout << "----- declare #" << iS << " : ";

	  hname  = "h_"+nameV[iV]+"_"+nameF[iF]+"_"+myPaths[iP].nameP+"_"+myPaths[iP].steps[iS].n;
	  cout << hname << endl;
	  hTemp = (TH1F*) infile->Get(hname);
	  h[iV][iF][iP].push_back(hTemp);
	}
      }
    }
  }

  /////////////////////////
  // BUILD TEFFICIENCIES //
  /////////////////////////

  TFile* outfile = new TFile("results/"+resultName+"/fits_"+resultName+".root","recreate");
  outfile->cd();

  TH1F *hNum, *hDen;

  const UInt_t nFunc=2;
  TF1 *f[nFunc];
  TString nameFunc[nFunc] = {"cb","sigmoid"};
  vector<TEfficiency*> pEff[nFunc][nV][nP];
  vector<TF1*>         fitEff[nFunc][nV][nP];

  TEfficiency* pEffTemp;
  TF1* fitEffTemp;

  f[0] = new TF1(nameFunc[0],evaluate,0,1000,5);
  f[0]->SetParName(0, "m0");
  f[0]->SetParName(1, "sigma");
  f[0]->SetParName(2, "alpha");
  f[0]->SetParName(3, "n");
  f[0]->SetParName(4, "norm");
  f[0]->SetParameter(0, 200);
  f[0]->SetParameter(1, 1);
  f[0]->SetParameter(2, 1);
  f[0]->SetParameter(3, 5);
  f[0]->SetParameter(4, 1);
  f[0]->SetParLimits(1, 0.01, 50);
  f[0]->SetParLimits(2, 0.01, 8);
  f[0]->SetParLimits(3, 1.1, 35);
  f[0]->SetParLimits(4, 0.6, 1);
  f[0]->SetLineWidth(2);
  
  f[1] = new TF1(nameFunc[1],evaluate2,0,1000,3);
  f[1]->SetParName(0, "midpoint");
  f[1]->SetParName(1, "steepness");
  f[1]->SetParName(2, "max");
  f[1]->SetParameter(0, 200);
  f[1]->SetParameter(1, 0.06);
  f[1]->SetParameter(2, 1);
  f[1]->SetParLimits(2, 0.995, 1);
  f[1]->SetLineWidth(2);

  // Set style //
  gROOT->Reset();
  setTDRStyle();
  gROOT->ForceStyle();

  Double_t eff95=0; 
	  
  // Loop over histograms
  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	
	// Get numerator and denominator histos
	hDen = h[iV][0][iP][iS];
	hNum = h[iV][1][iP][iS];

	for(UInt_t iFunc=0 ; iFunc<nFunc ; iFunc++) {

	  // Produce TEfficiency and fit it
	  if(hNum && hDen && TEfficiency::CheckConsistency(*hNum, *hDen) ) {
	    pEffTemp = new TEfficiency(*hNum,*hDen);
	    pEffTemp->
	      SetNameTitle( "t_"+TString(hNum->GetName())+nameFunc[iFunc], 
			    myPaths[iP].namePath+";"+nameAxis[iV]+";Efficiency" );

	    if( nameV[iV].Contains("met") || nameV[iV].Contains("pt") ) {

	      if(iFunc==0) {
		f[0]->SetParameter(0, 120);
		f[0]->SetParameter(1, 1);
		f[0]->SetParameter(2, 1);
		f[0]->SetParameter(3, 5);
		f[0]->SetParameter(4, 1);
		f[0]->SetParLimits(1, 0.01, 50);
		f[0]->SetParLimits(2, 0.01, 8);
		f[0]->SetParLimits(3, 1.1, 35);
		f[0]->SetParLimits(4, 0.6, 1);
	      }
	      else if(iFunc==1) {
		f[1]->SetParameter(0, 120);
		f[1]->SetParameter(1, 0.06);
		f[1]->SetParameter(2, 1);
		f[1]->SetParLimits(2, 0.995, 1);
	      }

	      pEffTemp->Fit(f[iFunc],"R");

	      eff95 = dichotomy(0.95, 0, 1000, 0.0000001, *f[iFunc], true);

	      fitEffTemp = 
		(TF1*)(pEffTemp->GetListOfFunctions()->FindObject(nameFunc[iFunc]));

	      fitEffTemp->SetLineColor(  myPaths[iP].steps[iS].C);
	      fitEffTemp->SetMarkerColor(myPaths[iP].steps[iS].C);
	      fitEffTemp->SetMarkerStyle(myPaths[iP].steps[iS].S);
	      
	    }
	  }

	  pEffTemp->SetLineColor(  myPaths[iP].steps[iS].C);
	  pEffTemp->SetMarkerColor(myPaths[iP].steps[iS].C);
	  pEffTemp->SetMarkerStyle(myPaths[iP].steps[iS].S);

	  pEffTemp->Write();
	  TCanvas c("c","c",0,0,600,600);
	  pEffTemp->Draw("AP");

	  pEff[iFunc][iV][iP].push_back(pEffTemp);
	  fitEff[iFunc][iV][iP].push_back(fitEffTemp);
	
	  gStyle->SetStatX(0.85);
	  gStyle->SetStatY(0.4);
	  gStyle->SetStatW(0.2);
	  gStyle->SetStatH(0.1);
	
	  Float_t nPass = (Float_t)hNum->Integral();
	  Float_t nTot  = (Float_t)hDen ->Integral();
	  Float_t globalEff = nTot!=0 ? nPass/nTot : -1.0;
	  TString s_globalEff = "#epsilon = "+TString(Form("%.1f",100*globalEff))+" %";
	
	  TString s_eff95 = "#epsilon = 95% @ "+TString(Form("%.0f", eff95))+" GeV";
	
	  TPaveText *pt2 = new TPaveText(0.58,0.15,0.85,0.22,"brNDC"); 
	  pt2->SetLineColor(1);
	  pt2->SetTextColor(1);
	  pt2->SetTextFont(42);
	  pt2->SetTextSize(0.03);
	  pt2->SetFillColor(kWhite);
	  pt2->SetShadowColor(kWhite);
	  //pt2->AddText(s_globalEff);
	  pt2->AddText(s_eff95);
	  if( nameV[iV].Contains("met") || 
	      nameV[iV].Contains("pt") ) {
	    pt2->Draw();
	  }
	
	  //c.Print("results/"+resultName+"/"+TString(hNum->GetName())+".png","png");
	  c.Print("results/"+resultName+"/"+TString(hNum->GetName())+"_"+nameFunc[iFunc]+".pdf","pdf");

	} // end loop: fit functions nFunc
      } // end loop: steps nS
    } // end loop: paths nP
  } // end loop: variables nV
  // end loop over histograms


  // Write histograms //
  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	nS = myPaths[iP].nSteps;
	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	  h[iV][iF][iP][iS]->Write();
	}
      }
    }
  }



  // PRODUCE PLOTS //
  //const UInt_t nSPlot = 7;
  //UInt_t idxSPlot[  nSPlot] = {0,1,2,3,4,5,7};
  //UInt_t idxFitStep[nSPlot] = {1,1,1,0,0,0,0};

  const UInt_t nVPlot = 3;
  UInt_t idxVPlot[nVPlot] = {0,1,2};

  gStyle->SetOptStat(0);

  TCanvas *cMerge;

  TString namePlot;

  for(UInt_t iP=0 ; iP<nP ; iP++) {

    for(UInt_t iVPlot=0 ; iVPlot<nVPlot ; iVPlot++) {
      
      UInt_t iV=idxVPlot[iVPlot];
      
      namePlot = "plot_"+nameV[iV]+"_"+myPaths[iP].nameP;
      cMerge = new TCanvas("c_"+namePlot,"c_"+namePlot,0,0,600,600);
      gStyle->SetOptStat(0);
      gPad->SetLogx();
      gPad->RangeAxis(0,0,300,1.05); // xmin, ymin, xmax, ymax
      gPad->Update();

      TLegend *leg = new TLegend(0.50,0.50,0.70,0.70);
      leg->SetFillColor(kWhite);
      leg->SetBorderSize(1);

      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) {

	theStep=myPaths[iP].steps[iS].n;
	UInt_t iF=0;
	if(theStep=="L1" || 
	   theStep=="MET" || 
	   theStep=="METC") {
	  iF=1;
	}

	if(iS==0) pEff[iF][iV][iP][iS]->Draw();
	else      pEff[iF][iV][iP][iS]->Draw("SAME");
	leg->AddEntry(pEff[iF][iV][iP][iS], myPaths[iP].steps[iS].t,"P");
      }

      leg->Draw();
      cMerge->Update();
      cMerge->Print("results/"+resultName+"/"+namePlot+".pdf","pdf");
    }
  }


  // Fit only 3 stages: L1, CaloMHT, entire path
  //const UInt_t nSPlotLess = 3;
  //UInt_t idxSPlotLess[  nSPlotLess] = {0,4,7};
  //UInt_t idxFitStepLess[nSPlotLess] = {1,1,0};

  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iVPlot=0 ; iVPlot<nVPlot ; iVPlot++) {
      
      UInt_t iV=idxVPlot[iVPlot];
      
      namePlot = "plot3stages_"+nameV[iV]+"_"+myPaths[iP].nameP;
      cMerge = new TCanvas("c_"+namePlot,"c_"+namePlot,0,0,600,600);
      gStyle->SetOptStat(0);
      gPad->SetLogx();
      gPad->RangeAxis(0,0,300,1.05); // xmin, ymin, xmax, ymax

      TLegend *leg = new TLegend(0.50,0.50,0.70,0.70);
      leg->SetFillColor(kWhite);
      leg->SetBorderSize(1);

      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) {
      
	theStep=myPaths[iP].steps[iS].n;
	UInt_t iF=0;
	if(theStep=="L1" || theStep=="MET" || theStep=="METC") {
	  iF=1;
	}

	if(theStep!="L1" && theStep!="MHT" && theStep!="Full") {
	  continue;
	}

	if(iS==0)     pEff[iF][iV][iP][iS]->Draw();
	else          pEff[iF][iV][iP][iS]->Draw("SAME");
	leg->AddEntry(pEff[iF][iV][iP][iS], myPaths[iP].steps[iS].t,"P");
      }

      leg->Draw();
      cMerge->Update();
      cMerge->Print("results/"+resultName+"/"+namePlot+".pdf","pdf");
    }
  }

  // END //
  outfile->Write();

  cout << "=== THE END ===" << endl
       << "==> [] o/" << endl
       << "++++++++++++++" << endl;

  return 0;
}


Double_t evaluate(double *x, double *par)
{ 
  double m = x[0];
  double m0 = par[0];
  double sigma = par[1];
  double alpha = par[2];
  double n = par[3];
  double norm = par[4];
  
  const double sqrtPiOver2 = 1.2533141373; // sqrt(pi/2)
  const double sqrt2 = 1.4142135624;

  Double_t sig = fabs((Double_t) sigma);
  Double_t t = (m - m0)/sig ;
  
  if (alpha < 0)
    t = -t;

  Double_t absAlpha = fabs(alpha / sig);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  Double_t b = absAlpha - n/absAlpha;

  Double_t aireGauche = (1 + ApproxErf( absAlpha / sqrt2 )) * sqrtPiOver2 ;
  Double_t aireDroite = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  Double_t aire = aireGauche + aireDroite;

  if ( t <= absAlpha ){
    return norm * (1 + ApproxErf( t / sqrt2 )) * sqrtPiOver2 / aire ;
  }
  else{
    return norm * (aireGauche +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / aire ;
  }
  
} 

Double_t ApproxErf(Double_t arg)
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return TMath::Erf(arg);
}

Double_t evaluate2(double *x, double *par)
{ 
  return par[2] / (1 + TMath::Exp(-par[1]*(x[0] - par[0])));
} 

Double_t dichotomy(double eff, double a0, double b0, double relErr,
		   TF1 f, bool verbose) 
{
  
  double dicho, effApprox, a, b;
  
  if(a0<b0) {
    a = a0;
    b = b0;
  } else if(a0>b0) {
    a = b0;
    b = a0;
  }
  else {
    cout << "PLEASE CHOOSE DIFFERENT VALUES FOR a AND b" << endl;
    return -999;
  }

  // Test bounds
  if( (f.Eval(a) > eff) || (f.Eval(b) < eff) ) {
    cout << "Bounds not large enough : eff(a)=" << f.Eval(a) 
	 << " ; eff(b)=" << f.Eval(b) << " ; tested eff=" << eff
	 << endl;
    return -999;
  }

  do {
    dicho = (a+b)/2 ;
    effApprox = f.Eval(dicho);

    if( effApprox < eff ) {
      a = dicho;
    } else {
      b = dicho;
    }
  }
  while( (fabs(effApprox-eff) / eff) > relErr );

  if(verbose) {
    cout << "relative precision asked (" << relErr*100 << " %) reached !"
	 << endl
	 << "found value of eT : " << dicho << " GeV" 
	 << endl
      //<< "efficiency value : " << 100*efficiency(dicho,mean,sigma,alpha,n,norm) << " %"
	 << endl;
  }

  return dicho;

}

pair<Int_t, Int_t> getStyle(TString name)
{

  if(name=="L1")    return make_pair(kBlack , kOpenSquare);
  if(name=="MET")   return make_pair(kBlue+2, kOpenTriangleUp);
  if(name=="METC")  return make_pair(kBlue  , kOpenTriangleDown);
  if(name=="METJ")  return make_pair(kCyan+2, kFullTriangleUp);
  if(name=="MHT")   return make_pair(kGreen+2,kFullTriangleDown);
  if(name=="PFMHT") return make_pair(kRed+2,  kOpenCircle);
  if(name=="PFMET") return make_pair(kRed,    kFullCircle);
  if(name=="Full")  return make_pair(kRed,    kFullCircle);

  return make_pair(kBlack,kFullCircle);
}
