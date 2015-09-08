#include "myIncludes.h"

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

Int_t myTrigger(TString resultName="v1_test", 
		TString json="DCS",
		TString binning="regular",
		TString offlineSel="TightMuon")
		//TString fitfunc="sigmoid")
{

  // Output log //
  ofstream outIneff("results/"+resultName+"/outIneff.txt");

  // External configuration
  //
  // Define JSON file
  bool applyJson=false;
  map<int, vector<pair<int, int> > > jsonMap;
  MAP_RLE mapRunLumiEvents;

  //Bool_t useSigmoid = fitfunc.Contains("sigmoid");
  //Bool_t useCB = fitfunc.Contains("CB");

  if(json=="DCS") {
    applyJson = true;
    jsonMap = readJSONFile("/user/ndaci/Data/json/json_DCSONLY_Run2015B.txt");
  }
  else if(json=="Prompt") {
    applyJson = true;
    jsonMap = readJSONFile("/user/ndaci/Data/json/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt");    
  }

  //////////////////
  // Define Chain //
  //////////////////

  TChain* ch = new TChain("tree");
  //ch->Add("/user/ndaci/Data/XMET/Run2015B/DoubleMuon/skim.root");
  //ch->Add("/user/ndaci/Data/XMET/Run2015B/SingleMuon/skim.root");

  if(resultName.Contains("0T"))       
    ch->Add("/user/ndaci/Data/XMET/Run2015C/SingleMuon/V2/skim_met30.root");
  else if(resultName.Contains("38T")) 
    ch->Add("/user/ndaci/Data/XMET/Run2015C/SingleMuon/38T_V5/skim_met30.root");
  else 
    ch->Add("/user/ndaci/Data/XMET/Run2015C/SingleMuon/38T_V5/skim_met30.root");

  ///////////////////////////////////
  // Define trigger interpretation //
  ///////////////////////////////////

  const UInt_t nV=6; // mumet, t1mumet, pfmet, t1pfmet, signaljetpt, signaljetNHfrac
  const UInt_t nF=2; // denom, num
  const UInt_t nP=6; // PFMNoMu90,120; PFM90,120; PFM170; CM200;
  UInt_t nS;

  // PFMNoMu: L1, MET, METClean, METJetID, MHT, PFMHT, PFMET, HLT bit
  // PFM:     L1, MET, METClean, METJetID, MHT, PFMHT, PFMET, HLT bit
  // PFM170:  L1, MET, METClean, METJetID, PFMET
  // CM200 :  L1, MET, METClean, METJetID

  //TString nameS[nS]={"L1","MET","METClean","METJetID","MHT","PFMHT","PFMET","bit"};
  //TString nameStep[nS]={"L1","MET","Cleaned MET", "JetID-cleaned MET", "MHT", "PFMHTNoMu", "PFMETNoMu", "Full path"};
  /*
  const UInt_t nS[nP] = {8,8,8,8,5,4};

  vector<TString> _myCol[nP];
  vector<double>  _thresh[nP];
  */

  STEP s_L1ETM60 ={.T=60,.c="hltL1extraParticles:MET:HLT",.n="L1",.t="L1",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};               

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
  //STEP s_bMET170 ={.T=0,   .c="hltmetwithmu170",           .n="Full", .t="Full path",.C=1,.S=1,.pt=0,.phi=0,.pass=0,.serial=0};

  PATH myPaths[nP];
  vector<STEP> vStepEmpty;
  myPaths[0]={.nameP="PFMNoMu90",.namePath="HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[0].steps.clear();
  myPaths[0].steps.push_back(s_L1ETM60);
  myPaths[0].steps.push_back(s_MET65);
  myPaths[0].steps.push_back(s_METC55);
  myPaths[0].steps.push_back(s_METJ55);
  myPaths[0].steps.push_back(s_MHT65);
  myPaths[0].steps.push_back(s_MuMHT90);
  myPaths[0].steps.push_back(s_MuMET90);
  myPaths[0].steps.push_back(s_bMuPFM90);

  myPaths[1]={.nameP="PFMNoMu120",.namePath="HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight",.nSteps=8,.steps=vStepEmpty};
  myPaths[1].steps.clear();
  myPaths[1].steps.push_back(s_L1ETM60);
  myPaths[1].steps.push_back(s_MET80);
  myPaths[1].steps.push_back(s_METC70);
  myPaths[1].steps.push_back(s_METJ70);
  myPaths[1].steps.push_back(s_MHT80);
  myPaths[1].steps.push_back(s_MuMHT120);
  myPaths[1].steps.push_back(s_MuMET120);
  myPaths[1].steps.push_back(s_bMuPFM120);

  myPaths[2]={.nameP="PFM90",.namePath="HLT_PFMET90_PFMHT90_IDTight",.nSteps=6,.steps=vStepEmpty};
  myPaths[2].steps.clear();
  myPaths[2].steps.push_back(s_L1ETM60);
  myPaths[2].steps.push_back(s_MET70);
  myPaths[2].steps.push_back(s_MHT70);
  myPaths[2].steps.push_back(s_PFMHT90);
  myPaths[2].steps.push_back(s_PFMET90);
  myPaths[2].steps.push_back(s_bPFM90);

  myPaths[3]={.nameP="PFM120",.namePath="HLT_PFMET120_PFMHT120_IDTight",.nSteps=6,.steps=vStepEmpty};
  myPaths[3].steps.clear();
  myPaths[3].steps.push_back(s_L1ETM60);
  myPaths[3].steps.push_back(s_MET90);
  myPaths[3].steps.push_back(s_MHT90);
  myPaths[3].steps.push_back(s_PFMHT120);
  myPaths[3].steps.push_back(s_PFMET120);
  myPaths[3].steps.push_back(s_bPFM120);

  myPaths[4]={.nameP="PFMET170",.namePath="HLT_PFMET170_JetIdCleaned",.nSteps=6,.steps=vStepEmpty};
  myPaths[4].steps.clear();
  myPaths[4].steps.push_back(s_L1ETM60);
  myPaths[4].steps.push_back(s_MET90);
  myPaths[4].steps.push_back(s_METC80);
  myPaths[4].steps.push_back(s_METJ80);
  myPaths[4].steps.push_back(s_PFMET170);
  myPaths[4].steps.push_back(s_bMET170);

  myPaths[5]={.nameP="CaloMET200",.namePath="HLT_CaloMET200_JetIdCleaned",.nSteps=3,.steps=vStepEmpty};
  myPaths[5].steps.clear();
  myPaths[5].steps.push_back(s_L1ETM60);
  myPaths[5].steps.push_back(s_MET210);
  myPaths[5].steps.push_back(s_METC200);

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

  // Trigger inputs //
  double _toPt, _toEta, _toPhi;
  TString _toCol, _toLab, _toPathFF, _toPathFT, _toPathTF, _toPathTT;

  // Trigger outputs
  //vector<double> _hlt_pt[nP], _hlt_phi[nP];
  //vector<bool>   _pass[nP];
  //vector<bool>   _serial[nP][nS];
  TString theColl, theStep;
  bool fired;
  
  ////////////////
  // HISTOGRAMS //
  ////////////////

  vector<TH1F*> h[nV][nF][nP];

  TString hname, title;
  TString nameV[nV]={"mumet","t1mumet","pfmet","t1pfmet","signaljetpt","signaljetNHfrac"};
  TString nameAxis[nV]={"Reco PFMETNoMu [GeV]",
			"Type1 PFMETNoMu [GeV]",
			"Reco PFMET [GeV]",
			"Type1 PFMET",
			"Leading PFJet p_{T}",
			"Leading PFJet NHEF"};

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

  TString nameF[nF]={"denom","num"};

  TH1F* hTemp;

  cout << "- Declare histograms" << endl;

  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths

	nS = myPaths[iP].nSteps;
	cout << "---- path: " << myPaths[iP].nameP << " ; nS=" << nS << endl;

	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths

	  cout << "----- declare #" << iS << " : ";

	  hname  = "h_"+nameV[iV]+"_"+nameF[iF]+"_"+myPaths[iP].nameP+"_"+myPaths[iP].steps[iS].n;
	  title = nameV[iV]+" "+nameF[iF]+" "+myPaths[iP].nameP+" "+myPaths[iP].steps[iS].n;

	  cout << hname << endl;

	  if(binning=="regular") {
	    hTemp = new TH1F(hname, title, xbins_reg[iV], xlow_reg[iV], xup_reg[iV]);				    
	  }
	  else {
	    hTemp = new TH1F(hname, title, xbins[iV]-1, v_xlow[iV]);
	  }

	  setStyle(hTemp, kBlack);
	  hTemp->SetXTitle(nameAxis[iV]);

	  h[iV][iF][iP].push_back(hTemp);
	}
      }
    }
  }

  ///////////////
  // SET CHAIN //
  ///////////////

  cout << "-Set Chain" << endl;

  // branch variables //
  Int_t           trig_obj_n;
  vector<double> *trig_obj_pt  = new vector<double> ();
  vector<double> *trig_obj_eta = new vector<double> ();
  vector<double> *trig_obj_phi = new vector<double> ();

  vector<string>  *trig_obj_col= new vector<string> ();
  //vector<vector<int> > *trig_obj_ids= new vector<> ();
  vector<string>  *trig_obj_lab= new vector<string> ();
  vector<string>  *trig_obj_path_FF= new vector<string> ();
  vector<string>  *trig_obj_path_FT= new vector<string> ();
  vector<string>  *trig_obj_path_TF= new vector<string> ();
  vector<string>  *trig_obj_path_TT= new vector<string> ();

  // other
  int32_t  puobs, putrue; 
  int32_t  wzid, l1id, l2id, i1id, i2id, i3id, mu1pid, mu2pid, mu1id, mu2id, el1pid, el2pid, el1id, el2id; 
  //uint32_t event, run, lumi;
  int event, run, lumi;
  uint32_t nvtx, nmuons, nelectrons, ntaus, ntightmuons, ntightelectrons, nphotons, njets, nbjets, nfatjets;
  uint8_t  hltmet90, hltmet120, hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel;
  uint8_t  flagcsctight, flaghbhenoise, flaghcallaser, flagecaltrig, flageebadsc, flagecallaser, flagtrkfail, flagtrkpog, flaghnoiseloose, flaghnoisetight, flaghnoisehilvl;
  double   pfmet, pfmetphi, t1pfmet, t1pfmetphi, pfmupt, pfmuphi, mumet, mumetphi, phmet, phmetphi, t1mumet, t1mumetphi, t1phmet, t1phmetphi;
  double   fatjetpt, fatjeteta, fatjetphi, fatjettau2, fatjettau1, fatjetCHfrac, fatjetNHfrac, fatjetEMfrac, fatjetCEMfrac, fatjetmetdphi, fatjetprmass, fatjetsdmass, fatjettrmass, fatjetftmass;
  double   signaljetpt, signaljeteta, signaljetphi, signaljetbtag, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi;
  double   secondjetpt, secondjeteta, secondjetphi, secondjetbtag, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetCEMfrac, secondjetmetdphi;
  double   thirdjetpt , thirdjeteta , thirdjetphi , thirdjetbtag , thirdjetCHfrac , thirdjetNHfrac , thirdjetEMfrac , thirdjetCEMfrac , thirdjetmetdphi ;
  double   jetjetdphi, jetmetdphimin, incjetmetdphimin;
  double   ht, dht, mht, alphat, apcjetmetmax, apcjetmetmin; 
  double   wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi, i1pt, i1eta, i1phi, i2pt, i2eta, i2phi, i3pt, i3eta, i3phi;
  double   mu1pt, mu1eta, mu1phi, mu2pt, mu2eta, mu2phi, el1pt, el1eta, el1phi, el2pt, el2eta, el2phi, phpt, pheta, phphi;
  double   zmass, zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt;
  double   loosephpt, loosepheta, loosephphi, loosephsieie, loosephrndiso;
  double   xsec, wgt, kfact, puwgt;

  // branch status //
  ch->SetBranchStatus("*", 1);  
  /*
  ch->SetBranchStatus("*", 0);  
  ch->SetBranchStatus("trig_obj_n",   1);
  ch->SetBranchStatus("trig_obj_pt",  1);
  ch->SetBranchStatus("trig_obj_eta", 1);
  ch->SetBranchStatus("trig_obj_phi", 1);

  ch->SetBranchStatus("trig_obj_eta", 1); 
  ch->SetBranchStatus("trig_obj_phi", 1); 
  ch->SetBranchStatus("trig_obj_col", 1); 
  //ch->SetBranchStatus("trig_obj_ids", 1); 
  ch->SetBranchStatus("trig_obj_lab", 1); 
  ch->SetBranchStatus("trig_obj_path_FF", 1); 
  ch->SetBranchStatus("trig_obj_path_FT", 1); 
  ch->SetBranchStatus("trig_obj_path_TF", 1); 
  ch->SetBranchStatus("trig_obj_path_TT", 1); 

  ch->SetBranchStatus("hltmet90", 1); // ", &hltmet90); // , &b_hltmet90);
  ch->SetBranchStatus("hltmet120", 1); // ", &hltmet120); // , &b_hltmet120);
  ch->SetBranchStatus("hltmetwithmu90", 1); // ", &hltmetwithmu90); // , &b_hltmetwithmu90);
  ch->SetBranchStatus("hltmetwithmu120", 1); // ", &hltmetwithmu120); // , &b_hltmetwithmu120);
  ch->SetBranchStatus("hltmetwithmu170", 1); // ", &hltmetwithmu170); // , &b_hltmetwithmu170);
  ch->SetBranchStatus("hltmetwithmu300", 1); // ", &hltmetwithmu300); // , &b_hltmetwithmu300);
  ch->SetBranchStatus("hltjetmet90", 1); // ", &hltjetmet90); // , &b_hltjetmet90);
  ch->SetBranchStatus("hltjetmet120", 1); // ", &hltjetmet120); // , &b_hltjetmet120);
  ch->SetBranchStatus("hltphoton165", 1); // ", &hltphoton165); // , &b_hltphoton165);
  ch->SetBranchStatus("hltphoton175", 1); // ", &hltphoton175); // , &b_hltphoton175);
  ch->SetBranchStatus("hltdoublemu", 1); // ", &hltdoublemu); // , &b_hltdoublemu);
  ch->SetBranchStatus("hltsinglemu", 1); // ", &hltsinglemu); // , &b_hltsinglemu);
  ch->SetBranchStatus("hltdoubleel", 1); // ", &hltdoubleel); // , &b_hltdoubleel);
  ch->SetBranchStatus("hltsingleel", 1); // ", &hltsingleel); // , &b_hltsingleel);

  ch->SetBranchStatus("event", 1); 
  ch->SetBranchStatus("run", 1); 
  ch->SetBranchStatus("lumi", 1); 
  ch->SetBranchStatus("puwgt", 1); 
  ch->SetBranchStatus("puobs", 1); 
  ch->SetBranchStatus("putrue", 1); 
  ch->SetBranchStatus("nvtx", 1); 

  ch->SetBranchStatus("mumet", 1);
  ch->SetBranchStatus("t1mumet", 1);
  ch->SetBranchStatus("pfmet", 1);
  ch->SetBranchStatus("t1pfmet", 1);
  ch->SetBranchStatus("signaljetpt", 1);
  ch->SetBranchStatus("signaljetNHfrac", 1);

  ch->SetBranchStatus("ntightmuons", 1);
  */

  // branch address //

  // trigger objects
  ch->SetBranchAddress("trig_obj_n", &trig_obj_n); // , &b_trig_obj_n);
  ch->SetBranchAddress("trig_obj_pt", &trig_obj_pt); // , &b_trig_obj_pt);
  ch->SetBranchAddress("trig_obj_eta", &trig_obj_eta); // , &b_trig_obj_eta);
  ch->SetBranchAddress("trig_obj_phi", &trig_obj_phi); // , &b_trig_obj_phi);
  //
  ch->SetBranchAddress("trig_obj_col", &trig_obj_col); // , &b_trig_obj_col);
  //ch->SetBranchAddress("trig_obj_ids", &trig_obj_ids); // , &b_trig_obj_ids);
  ch->SetBranchAddress("trig_obj_lab", &trig_obj_lab); // , &b_trig_obj_lab);
  ch->SetBranchAddress("trig_obj_path_FF", &trig_obj_path_FF); // , &b_trig_obj_path_FF);
  ch->SetBranchAddress("trig_obj_path_FT", &trig_obj_path_FT); // , &b_trig_obj_path_FT);
  ch->SetBranchAddress("trig_obj_path_TF", &trig_obj_path_TF); // , &b_trig_obj_path_TF);
  ch->SetBranchAddress("trig_obj_path_TT", &trig_obj_path_TT); // , &b_trig_obj_path_TT);

  // global informations
  ch->SetBranchAddress("event", &event); // , &b_event);
  ch->SetBranchAddress("run", &run); // , &b_run);
  ch->SetBranchAddress("lumi", &lumi); // , &b_lumi);

  // pile-up
  ch->SetBranchAddress("puwgt", &puwgt); // , &b_puwgt);
  ch->SetBranchAddress("puobs", &puobs); // , &b_puobs);
  ch->SetBranchAddress("putrue", &putrue); // , &b_putrue);
  ch->SetBranchAddress("nvtx", &nvtx); // , &b_nvtx);

  // trigger bits
  ch->SetBranchAddress("hltmet90", &hltmet90); // , &b_hltmet90);
  ch->SetBranchAddress("hltmet120", &hltmet120); // , &b_hltmet120);
  ch->SetBranchAddress("hltmetwithmu90", &hltmetwithmu90); // , &b_hltmetwithmu90);
  ch->SetBranchAddress("hltmetwithmu120", &hltmetwithmu120); // , &b_hltmetwithmu120);
  ch->SetBranchAddress("hltmetwithmu170", &hltmetwithmu170); // , &b_hltmetwithmu170);
  ch->SetBranchAddress("hltmetwithmu300", &hltmetwithmu300); // , &b_hltmetwithmu300);
  ch->SetBranchAddress("hltjetmet90", &hltjetmet90); // , &b_hltjetmet90);
  ch->SetBranchAddress("hltjetmet120", &hltjetmet120); // , &b_hltjetmet120);
  ch->SetBranchAddress("hltphoton165", &hltphoton165); // , &b_hltphoton165);
  ch->SetBranchAddress("hltphoton175", &hltphoton175); // , &b_hltphoton175);
  ch->SetBranchAddress("hltdoublemu", &hltdoublemu); // , &b_hltdoublemu);
  ch->SetBranchAddress("hltsinglemu", &hltsinglemu); // , &b_hltsinglemu);
  ch->SetBranchAddress("hltdoubleel", &hltdoubleel); // , &b_hltdoubleel);
  ch->SetBranchAddress("hltsingleel", &hltsingleel); // , &b_hltsingleel);

  ch->SetBranchAddress("pfmet", &pfmet); // , &b_pfmet);
  ch->SetBranchAddress("pfmetphi", &pfmetphi); // , &b_pfmetphi);
  ch->SetBranchAddress("t1pfmet", &t1pfmet); // , &b_t1pfmet);
  ch->SetBranchAddress("t1pfmetphi", &t1pfmetphi); // , &b_t1pfmetphi);
  ch->SetBranchAddress("pfmupt", &pfmupt); // , &b_pfmupt);
  ch->SetBranchAddress("pfmuphi", &pfmuphi); // , &b_pfmuphi);
  ch->SetBranchAddress("mumet", &mumet); // , &b_mumet);
  ch->SetBranchAddress("mumetphi", &mumetphi); // , &b_mumetphi);
  ch->SetBranchAddress("phmet", &phmet); // , &b_phmet);
  ch->SetBranchAddress("phmetphi", &phmetphi); // , &b_phmetphi);
  ch->SetBranchAddress("t1mumet", &t1mumet); // , &b_t1mumet);
  ch->SetBranchAddress("t1mumetphi", &t1mumetphi); // , &b_t1mumetphi);
  ch->SetBranchAddress("t1phmet", &t1phmet); // , &b_t1phmet);
  ch->SetBranchAddress("t1phmetphi", &t1phmetphi); // , &b_t1phmetphi);

  ch->SetBranchAddress("signaljetpt", &signaljetpt); // , &b_signaljetpt);
  ch->SetBranchAddress("signaljeteta", &signaljeteta); // , &b_signaljeteta);
  ch->SetBranchAddress("signaljetphi", &signaljetphi); // , &b_signaljetphi);
  ch->SetBranchAddress("signaljetbtag", &signaljetbtag); // , &b_signaljetbtag);
  ch->SetBranchAddress("signaljetCHfrac", &signaljetCHfrac); // , &b_signaljetCHfrac);
  ch->SetBranchAddress("signaljetNHfrac", &signaljetNHfrac); // , &b_signaljetNHfrac);
  ch->SetBranchAddress("signaljetEMfrac", &signaljetEMfrac); // , &b_signaljetEMfrac);
  ch->SetBranchAddress("signaljetCEMfrac", &signaljetCEMfrac); // , &b_signaljetCEMfrac);
  ch->SetBranchAddress("signaljetmetdphi", &signaljetmetdphi); // , &b_signaljetmetdphi);
  ch->SetBranchAddress("secondjetpt", &secondjetpt); // , &b_secondjetpt);
  ch->SetBranchAddress("secondjeteta", &secondjeteta); // , &b_secondjeteta);
  ch->SetBranchAddress("secondjetphi", &secondjetphi); // , &b_secondjetphi);
  ch->SetBranchAddress("secondjetbtag", &secondjetbtag); // , &b_secondjetbtag);
  ch->SetBranchAddress("secondjetCHfrac", &secondjetCHfrac); // , &b_secondjetCHfrac);
  ch->SetBranchAddress("secondjetNHfrac", &secondjetNHfrac); // , &b_secondjetNHfrac);
  ch->SetBranchAddress("secondjetEMfrac", &secondjetEMfrac); // , &b_secondjetEMfrac);
  ch->SetBranchAddress("secondjetCEMfrac", &secondjetCEMfrac); // , &b_secondjetCEMfrac);
  ch->SetBranchAddress("secondjetmetdphi", &secondjetmetdphi); // , &b_secondjetmetdphi);
  ch->SetBranchAddress("thirdjetpt", &thirdjetpt); // , &b_thirdjetpt);
  ch->SetBranchAddress("thirdjeteta", &thirdjeteta); // , &b_thirdjeteta);
  ch->SetBranchAddress("thirdjetphi", &thirdjetphi); // , &b_thirdjetphi);
  ch->SetBranchAddress("thirdjetbtag", &thirdjetbtag); // , &b_thirdjetbtag);
  ch->SetBranchAddress("thirdjetCHfrac", &thirdjetCHfrac); // , &b_thirdjetCHfrac);
  ch->SetBranchAddress("thirdjetNHfrac", &thirdjetNHfrac); // , &b_thirdjetNHfrac);
  ch->SetBranchAddress("thirdjetEMfrac", &thirdjetEMfrac); // , &b_thirdjetEMfrac);
  ch->SetBranchAddress("thirdjetCEMfrac", &thirdjetCEMfrac); // , &b_thirdjetCEMfrac);
  ch->SetBranchAddress("thirdjetmetdphi", &thirdjetmetdphi); // , &b_thirdjetmetdphi);

  ch->SetBranchAddress("jetjetdphi", &jetjetdphi); // , &b_jetjetdphi);
  ch->SetBranchAddress("jetmetdphimin", &jetmetdphimin); // , &b_jetmetdphimin);
  ch->SetBranchAddress("incjetmetdphimin", &incjetmetdphimin); // , &b_incjetmetdphimin);
  ch->SetBranchAddress("ht", &ht); // , &b_ht);
  ch->SetBranchAddress("dht", &dht); // , &b_dht);
  ch->SetBranchAddress("mht", &mht); // , &b_mht);
  ch->SetBranchAddress("alphat", &alphat); // , &b_alphat);
  ch->SetBranchAddress("apcjetmetmax", &apcjetmetmax); // , &b_apcjetmetmax);
  ch->SetBranchAddress("apcjetmetmin", &apcjetmetmin); // , &b_apcjetmetmin);
  ch->SetBranchAddress("mu1pid", &mu1pid); // , &b_mu1pid);
  ch->SetBranchAddress("mu1pt", &mu1pt); // , &b_mu1pt);
  ch->SetBranchAddress("mu1eta", &mu1eta); // , &b_mu1eta);
  ch->SetBranchAddress("mu1phi", &mu1phi); // , &b_mu1phi);
  ch->SetBranchAddress("mu1id", &mu1id); // , &b_mu1id);
  ch->SetBranchAddress("mu2pid", &mu2pid); // , &b_mu2pid);
  ch->SetBranchAddress("mu2pt", &mu2pt); // , &b_mu2pt);
  ch->SetBranchAddress("mu2eta", &mu2eta); // , &b_mu2eta);
  ch->SetBranchAddress("mu2phi", &mu2phi); // , &b_mu2phi);
  ch->SetBranchAddress("mu2id", &mu2id); // , &b_mu2id);

  ch->SetBranchAddress("nmuons", &nmuons); // , &b_nmuons);
  ch->SetBranchAddress("nelectrons", &nelectrons); // , &b_nelectrons);
  ch->SetBranchAddress("ntightmuons", &ntightmuons); // , &b_ntightmuons);
  ch->SetBranchAddress("ntightelectrons", &ntightelectrons); // , &b_ntightelectrons);
  ch->SetBranchAddress("ntaus", &ntaus); // , &b_ntaus);
  ch->SetBranchAddress("njets", &njets); // , &b_njets);
  ch->SetBranchAddress("nbjets", &nbjets); // , &b_nbjets);
  ch->SetBranchAddress("nfatjets", &nfatjets); // , &b_nfatjets);
  ch->SetBranchAddress("nphotons", &nphotons); // , &b_nphotons);

  /*
  ch->SetBranchAddress("flagcsctight", &flagcsctight); // , &b_flagcsctight);
  ch->SetBranchAddress("flaghbhenoise", &flaghbhenoise); // , &b_flaghbhenoise);
  ch->SetBranchAddress("flaghcallaser", &flaghcallaser); // , &b_flaghcallaser);
  ch->SetBranchAddress("flagecaltrig", &flagecaltrig); // , &b_flagecaltrig);
  ch->SetBranchAddress("flageebadsc", &flageebadsc); // , &b_flageebadsc);
  ch->SetBranchAddress("flagecallaser", &flagecallaser); // , &b_flagecallaser);
  ch->SetBranchAddress("flagtrkfail", &flagtrkfail); // , &b_flagtrkfail);
  ch->SetBranchAddress("flagtrkpog", &flagtrkpog); // , &b_flagtrkpog);
  ch->SetBranchAddress("flaghnoiseloose", &flaghnoiseloose); // , &b_flaghnoiseloose);
  ch->SetBranchAddress("flaghnoisetight", &flaghnoisetight); // , &b_flaghnoisetight);
  ch->SetBranchAddress("flaghnoisehilvl", &flaghnoisehilvl); // , &b_flaghnoisehilvl);
  */

  /////////////////////
  // LOOP OVER CHAIN //
  /////////////////////

  UInt_t entries=ch->GetEntries();
  bool printOut=false;
  cout << "Start processing: " << entries << " entries." << endl;

  // initialize input variables
  double _var[nV] = {0, 0, 0, 0, 0, 0};
  UInt_t nIneff=0;
  UInt_t nEff=0;

  cout << "- Start looping over the chain" << endl;

  // START LOOP //
  for(UInt_t iE=0 ; iE<entries ; iE++) {

    // PRINT OUT //
    if(iE%1000==0) printOut=true;
    else printOut=false;
    if(printOut) {
      cout << "- processing entry #" 
	   << iE << "/" << entries
	   << endl;
    }

    // INITIALIZE //
    trig_obj_n = 0;
    trig_obj_pt->clear();
    trig_obj_eta->clear();
    trig_obj_phi->clear();
    //
    trig_obj_col->clear();
    trig_obj_lab->clear();
    trig_obj_path_FF->clear();
    trig_obj_path_FT->clear();
    trig_obj_path_TF->clear();
    trig_obj_path_TT->clear();
    //
    /*
    for(UInt_t i=0;i<trig_obj_ids.size();i++) {
      (trig_obj_ids)[i].clear();
    }
    trig_obj_ids.clear();
    */
    //
    hltmet90 = hltmet120 = hltmetwithmu90 = hltmetwithmu120 = hltmetwithmu170 = hltmetwithmu300 = hltjetmet90 = hltjetmet120 = hltphoton165 = hltphoton175 = hltdoublemu = hltsinglemu = hltdoubleel = hltsingleel = 0;
    ntightmuons=0;
    //
    // end initialization //
    

    // GET ENTRY //
    ch->GetEntry(iE);

    // json selection
    if(applyJson) {
      if( ! AcceptEventByRunAndLumiSection(run, lumi, jsonMap) ) {
	continue;
      }
    }

    // output json
    if(mapRunLumiEvents[run][lumi][event]==1)
      continue;
    else mapRunLumiEvents[run][lumi][event]=1;

    // event selection
    if(resultName.Contains("25ns")) {
      if(run==254833) continue;
    }
    
    if(offlineSel.Contains("TightMuon")) {
      if(ntightmuons==0) continue;
    }

    if(offlineSel.Contains("HLTMu")) {
      if(!hltsinglemu) continue;      
    }

    if(offlineSel.Contains("HLTDiMu")) {
      if(!hltdoublemu) continue;      
    }

        // print out every 1000 events
    if(printOut) {
      cout << "- trig_obj_n=" << trig_obj_n 
	//<< " trig_obj_pt->size()=" << trig_obj_pt->size()
	//<< " trig_obj_eta->size()=" << trig_obj_eta->size()
	//<< " trig_obj_phi->size()=" << trig_obj_phi->size()
	   << endl;
    }

    // get x-axis variables //
    _var[0] = mumet;
    _var[1] = t1mumet;
    _var[2] = pfmet;
    _var[3] = t1pfmet;
    _var[4] = signaljetpt;
    _var[5] = signaljetNHfrac;


    // PROCESS TRIGGER OBJECTS //

    // initialize trigger output
    for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	myPaths[iP].steps[iS].pt = myPaths[iP].steps[iS].phi = 0;
	myPaths[iP].steps[iS].pass = false;
      }
    }
    
    // loop: trigger objects
    for(UInt_t iObj=0 ; iObj<(UInt_t)trig_obj_n ; iObj++) {
      
      // Get object's properties
      _toPt  = (*trig_obj_pt )[iObj];
      _toEta = (*trig_obj_eta)[iObj];
      _toPhi = (*trig_obj_phi)[iObj];
      //
      _toCol    = (TString)(*trig_obj_col)[iObj];

      // loop: paths
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	nS = myPaths[iP].nSteps;
	// loop: steps
	for(UInt_t iS=0 ; iS<nS ; iS++) {
	  if(_toCol==myPaths[iP].steps[iS].c) {
	    myPaths[iP].steps[iS].pt  = _toPt;
	    myPaths[iP].steps[iS].phi = _toPhi;
	  }
	} // end loop: steps
      } // end loop: paths

    } // end loop: trigger objects

    // trigger outputs
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      //
      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) {
	//
	theColl=myPaths[iP].steps[iS].c;
	theStep=myPaths[iP].steps[iS].n;
	fired=false;
	//
	if(theStep=="Full") { // check trigger bit
	  if(     theColl=="hltmet90")        fired=hltmet90;
	  else if(theColl=="hltmet120")       fired=hltmet120;
	  else if(theColl=="hltmetwithmu90")  fired=hltmetwithmu90;
	  else if(theColl=="hltmetwithmu120") fired=hltmetwithmu120;
	  else if(theColl=="hltmetwithmu170") fired=hltmetwithmu170;
	  else fired=false;
	}
	//
	else { // check trigger objects
	  fired = (myPaths[iP].steps[iS].pt>myPaths[iP].steps[iS].T);
	}
	//
	myPaths[iP].steps[iS].pass = fired;
      }
    }

    // serial trigger
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      nS = myPaths[iP].nSteps;
      for(UInt_t iS=0 ; iS<nS ; iS++) {
	if(iS==0 || iS==nS-1) {
	  myPaths[iP].steps[iS].serial = myPaths[iP].steps[iS].pass;
	}
	else {
	  myPaths[iP].steps[iS].serial = myPaths[iP].steps[iS-1].serial && myPaths[iP].steps[iS].pass;
	}
      }
    }

    // FILL HISTOGRAMS //
    for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables

      // forget about energy fractions if MET<=200
      if(nameV[iV].Contains("frac") && mumet<=200) continue;

      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	nS = myPaths[iP].nSteps;
	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths

	  /*
	  // conditional denominator (eff1)
	  if(iS==0 || iS==nS-1) { // L1 or entire path
	    h[iV][0][iP][iS]->Fill(_var[iV]);
	  }
	  else { // intermediate HLT filters
	    if(_serial[iP][iS-1]) // events that fire up to step iS-1
	      h[iV][0][iP][iS]->Fill(_var[iV]);
	  }
	  */

	  // denominator (eff2)
	  h[iV][0][iP][iS]->Fill(_var[iV]);
	  
	  // numerator
	  if(myPaths[iP].steps[iS].serial) { // event fired step iS of path iP
	    h[iV][1][iP][iS]->Fill(_var[iV]);
	  }

	} // loop:nS
      } // loop:nP
    } // loop:nV
    // end fill histograms //

    // INVESTIGATE INEFFICIENCIES //
    if(pfmet>400) {
      if(!hltmet120 || !hltmetwithmu120 || !hltmetwithmu170) {
	nIneff++ ;
	outIneff << "INEFF: run " << run 
		 << " lumi "      << lumi
		 << " event "     << event
		 << endl
		 << "hltmet120:"  << (hltmet120 ? 1 : 0)
		 << " hltmetwithmu120:" << (hltmetwithmu120 ? 1 : 0)
		 << " hltmetwithmu170:" << (hltmetwithmu170 ? 1 : 0)
		 << endl << endl
		 << "1) Trigger objects" << endl
		 << "L1ETM="     << myPaths[0].steps[0].pt
		 << "   phi=" << myPaths[0].steps[0].phi 
		 << endl
		 << "MET="       << myPaths[0].steps[1].pt
 		 << "   phi=" << myPaths[0].steps[1].phi 
		 << endl
		 << "HBHE-Cleaned MET="       << myPaths[0].steps[2].pt
 		 << "   phi=" << myPaths[0].steps[2].phi 
		 << endl
		 << "JetID-Cleaned MET="       << myPaths[0].steps[3].pt
 		 << "   phi=" << myPaths[0].steps[3].phi 
		 << endl
		 << "MHT="       << myPaths[0].steps[4].pt
 		 << "   phi=" << myPaths[0].steps[4].phi 
		 << endl
		 << "PFMHTNoMu="       << myPaths[0].steps[5].pt
 		 << "   phi=" << myPaths[0].steps[5].phi 
		 << endl
		 << "PFMETNoMu="       << myPaths[0].steps[6].pt
 		 << "   phi=" << myPaths[0].steps[6].phi 
		 << endl
		 << "PFMHT="       << myPaths[2].steps[3].pt
 		 << "   phi=" << myPaths[2].steps[3].phi 
		 << endl
		 << "PFMET="       << myPaths[2].steps[4].pt
 		 << "   phi=" << myPaths[2].steps[4].phi 
		 << endl << endl
		 << "2) Offline MET objects"
		 << endl
		 << "mumet="     << mumet
		 << "   phi=" << mumetphi   
		 << endl
		 << "t1mumet="  << t1mumet
		 << "   phi=" << t1mumetphi 
		 << endl
		 << "pfmet="     << pfmet
		 << "   phi=" << pfmetphi   
		 << endl
		 << "t1pfmet="  << t1pfmet
		 << "   phi=" << t1pfmetphi 
		 << endl << endl
		 << "3) Offline objects"
		 << endl
		 << "Jet1(pt,eta,phi)="
		 << "(" << signaljetpt 
		 << "," << signaljeteta 
		 << "," << signaljetphi 
		 << ")" << endl
		 << "Jet2(pt,eta,phi)="
		 << "(" << secondjetpt 
		 << "," << secondjeteta 
		 << "," << secondjetphi 
		 << ")" << endl
		 <<"Jet3(pt,eta,phi)="
		 << "(" << thirdjetpt 
		 << "," << thirdjeteta 
		 << "," << thirdjetphi 
		 << ")" << endl
		 << "Mu1(pt,eta,phi)="
		 << "(" << mu1pt 
		 << "," << mu1eta 
		 << "," << mu1phi 
		 << ")" << endl
		 << "Mu2(pt,eta,phi)="
		 << "(" << mu2pt 
		 << "," << mu2eta 
		 << "," << mu2phi 
		 << ")" << endl
		 << "El1(pt,eta,phi)="
		 << "(" << el1pt 
		 << "," << el1eta 
		 << "," << el1phi 
		 << ")" << endl
		 << "El2(pt,eta,phi)="
		 << "(" << el2pt 
		 << "," << el2eta 
		 << "," << el2phi 
		 << ")" << endl
		 << "Ph1(pt,eta,phi)="
		 << "(" << loosephpt 
		 << "," << loosepheta 
		 << "," << loosephphi 
		 << ")" << endl
		 << endl;
      }
      else nEff++ ;
    }


  } // end loop:entries


  /////////////////////////
  // BUILD TEFFICIENCIES //
  /////////////////////////
  TFile* outfile = new TFile("results/"+resultName+"/f_"+resultName+".root","recreate");
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
  f[0]->SetParameter(0, 120);
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
  f[1]->SetParameter(0, 120);
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


  // Final printouts //
  cout << "INEFFICIENCY SUMMARY :"
       << " nIneff=" << nIneff
       << " nEff="   << nEff
       << endl;

  outIneff << "INEFFICIENCY SUMMARY :"
	   << " nIneff=" << nIneff
	   << " nEff="   << nEff
	   << endl;

  // END //
  outfile->Write();
  delete ch;
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
