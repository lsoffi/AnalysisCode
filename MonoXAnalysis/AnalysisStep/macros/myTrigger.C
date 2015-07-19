#include "myIncludes.h"

Double_t evaluate(double *x, double *par);
Double_t evaluate2(double *x, double *par);
Double_t ApproxErf(Double_t arg);
Double_t dichotomy(double eff, double a0, double b0, double relErr,
		   TF1 f, bool verbose);

Int_t myTrigger(TString resultName="v1_test", 
		TString json="DCS",
		TString binning="regular",
		TString offlineSel="TightMuon",
		TString fitfunc="sigmoid")
{

  // External configuration
  //
  // Define JSON file
  bool applyJson=false;
  map<int, vector<pair<int, int> > > jsonMap;

  Bool_t useSigmoid = fitfunc.Contains("sigmoid");
  Bool_t useCB = fitfunc.Contains("CB");

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
  ch->Add("/user/ndaci/Data/XMET/Run2015B/SingleMuon/skim.root");


  ///////////////////////////////////
  // Define trigger interpretation //
  ///////////////////////////////////

  const UInt_t nV=4; // mumet, t1mumet, signaljetpt, signaljetNHfrac
  const UInt_t nF=2; // denom, num
  const UInt_t nP=2; // 90GeV, 120GeV
  const UInt_t nS=8; // L1, MET, METClean, METJetID, MHT, PFMHT, PFMET, HLT bit

  // Trigger inputs //
  double _toPt, _toEta, _toPhi;
  TString _toCol, _toLab, _toPathFF, _toPathFT, _toPathTF, _toPathTT;

  // Trigger outputs
  double _hlt_pt[nS], _hlt_phi[nS];
  bool _pass[nP][nS];
  bool _serial[nP][nS];
  double _thresh[nP][nS] = { {60, 65, 55, 55, 65,  90,  90, 0},
			     {60, 80, 70, 70, 80, 120, 120, 0} };
  TString _myCol[nS] = {"hltL1extraParticles:MET:HLT", 
			"hltMet::HLT", 
			"hltMetClean::HLT",
			"hltMetCleanUsingJetID::HLT",
			"hltMht::HLT",
			"hltPFMHTNoMuTightID::HLT",
			"hltPFMETNoMuProducer::HLT",
			"bit"};
  double _var[nV] = {0, 0, 0, 0};


  ////////////////
  // HISTOGRAMS //
  ////////////////

  TH1F* h[nV][nF][nP][nS];

  TString name, title;
  TString nameV[nV]={"mumet","t1mumet","signaljetpt","signaljetNHfrac"};
  TString nameAxis[nV]={"Reco PFMETNoMu [GeV]",
			"Type1 PFMETNoMu [GeV]",
			"Leading PFJet p_{T}",
			"Leading PFJet NHEF"};

  int   xbins_reg[nV] = {40,  40,  40,  50};
  float xlow_reg[ nV] = {100, 100, 100, 0};
  float xup_reg[  nV] = {900, 900, 900, 1};

  const UInt_t xbins[nV]    = {21, 21, 21, 27};

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
  
  float* v_xlow[nV] = {bins_met , bins_met , bins_met , bins_nhef};

  TString nameF[nF]={"denom","num"};
  TString nameP[nP]={"MET90","MET120"};
  TString namePath[nP]={"HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
			"HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight"};
  TString nameS[nS]={"L1","MET","METClean","METJetID","MHT","PFMHT","PFMET","bit"};

  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths

	  name  = "h_"+nameV[iV]+"_"+nameF[iF]+"_"+nameP[iP]+"_"+nameS[iS];
	  title = nameV[iV]+" "+nameF[iF]+" "+nameP[iP]+" "+nameS[iS];

	  if(binning=="regular") {
	    h[iV][iF][iP][iS] = 
	      new TH1F(name, title, xbins_reg[iV], xlow_reg[iV], xup_reg[iV]);
	  }
	  else {
	    h[iV][iF][iP][iS] = 
	      new TH1F(name, title, xbins[iV]-1, v_xlow[iV]);
	  }

	  setStyle(h[iV][iF][iP][iS], kBlack);
	  h[iV][iF][iP][iS]->SetXTitle(nameAxis[iV]);
	}
      }
    }
  }


  ///////////////
  // SET CHAIN //
  ///////////////

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
  ch->SetBranchStatus("signaljetpt", 1);
  ch->SetBranchStatus("signaljetNHfrac", 1);

  ch->SetBranchStatus("ntightmuons", 1);

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

  // x-axis variables
  ch->SetBranchAddress("mumet", &mumet); // , &b_mumet);
  ch->SetBranchAddress("t1mumet", &t1mumet); // , &b_mumet);
  ch->SetBranchAddress("signaljetpt", &signaljetpt); // , &b_signaljetpt);
  ch->SetBranchAddress("signaljetNHfrac", &signaljetNHfrac); // , &b_signaljetNHfrac);

  // selection variables
  ch->SetBranchAddress("ntightmuons", &ntightmuons); // , &b_ntightmuons);

  /*
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

    if(offlineSel=="TightMuon") {
      if(ntightmuons==0) continue;
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
    _var[2] = signaljetpt;
    _var[3] = signaljetNHfrac;


    // PROCESS TRIGGER OBJECTS //

    // initialize trigger outputs
    for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
      _hlt_pt[iS] = _hlt_phi[iS] = 0;
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	_pass[iP][iS] = false;
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

      // loop: steps in paths
      for(UInt_t iS=0 ; iS<nS ; iS++) {
	if(_toCol==_myCol[iS]) {
	  _hlt_pt[iS]  = _toPt;
	  _hlt_phi[iS] = _toPhi;
	}
      } //loop:nS

    } // end loop: trigger objects

    // trigger outputs
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      for(UInt_t iP=0 ; iP<nP ; iP++) {
	// check trigger objects
	if(_hlt_pt[iS]>_thresh[iP][iS]) {
	  _pass[iP][iS] = true;
	}
	// check trigger bit
	if(iS==nS-1) {
	  if(iP==0 && hltmet90)  _pass[iP][iS]=true;
	  if(iP==1 && hltmet120) _pass[iP][iS]=true;
	}
      }
    }

    // serial trigger
    for(UInt_t iP=0 ; iP<nP ; iP++) {
      for(UInt_t iS=0 ; iS<nS ; iS++) {
	if(iS==0 || iS==nS-1) _serial[iP][iS] = _pass[iP][iS];
	else _serial[iP][iS] = _serial[iP][iS-1] && _pass[iP][iS];
      }
    }

    // FILL HISTOGRAMS //
    for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
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
	  if(_serial[iP][iS]) { // event fired step iS of path iP
	    h[iV][1][iP][iS]->Fill(_var[iV]);
	  }

	} // loop:nS
      } // loop:nP
    } // loop:nV
    // end fill histograms //

  } // end loop:entries


  /////////////////////////
  // BUILD TEFFICIENCIES //
  /////////////////////////

  TFile* outfile = new TFile("results/"+resultName+"/f_"+resultName+".root","recreate");
  outfile->cd();

  TH1F *hNum, *hDen;
  TEfficiency *pEff;

  TF1 *f1 = new TF1("fit",evaluate,0,1000,5);
  f1->SetParName(0, "m0");
  f1->SetParName(1, "sigma");
  f1->SetParName(2, "alpha");
  f1->SetParName(3, "n");
  f1->SetParName(4, "norm");
  f1->SetParameter(0, 120);
  f1->SetParameter(1, 1);
  f1->SetParameter(2, 1);
  f1->SetParameter(3, 5);
  f1->SetParameter(4, 1);
  f1->SetParLimits(1, 0.01, 50);
  f1->SetParLimits(2, 0.01, 8);
  f1->SetParLimits(3, 1.1, 35);
  f1->SetParLimits(4, 0.6, 1);
  f1->SetLineWidth(2);
  
  TF1 *f2 = new TF1("fit2",evaluate2,0,1000,3);
  f2->SetParName(0, "midpoint");
  f2->SetParName(1, "steepness");
  f2->SetParName(2, "max");
  f2->SetParameter(0, 120);
  f2->SetParameter(1, 0.06);
  f2->SetParameter(2, 1);
  f2->SetParLimits(2, 0.995, 1);
  f2->SetLineWidth(2);

  // Set style //
  gROOT->Reset();
  setTDRStyle();
  gROOT->ForceStyle();

  // Loop over histograms
  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
      for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	
	// Get numerator and denominator histos
	hDen = h[iV][0][iP][iS];
	hNum = h[iV][1][iP][iS];

	// Produce TEfficiency and fit it
	if(hNum && hDen && TEfficiency::CheckConsistency(*hNum, *hDen) ) {
	  pEff = new TEfficiency(*hNum,*hDen);
	  pEff->SetNameTitle( "t_"+TString(hNum->GetName()) , 
			      namePath[iP] );

	  if( nameV[iV].Contains("met") || nameV[iV].Contains("pt") ) {
	    if(useSigmoid) pEff->Fit(f2,"R"); // use function's definition Range
	    else if(useCB) pEff->Fit(f1,"R"); // use function's definition Range
	  }
	  
	  Double_t eff95=0; 
	  if( nameV[iV].Contains("met") || nameV[iV].Contains("pt") ) {
	    if(useSigmoid) eff95 = dichotomy(0.95, 0, 1000, 0.0000001, *f2, true);
	    else if(useCB) eff95 = dichotomy(0.95, 0, 1000, 0.0000001, *f1, true);
	  }

	  pEff->Write();
	  TCanvas c("c","c",0,0,600,600);
	  pEff->Draw("AP");

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
	  c.Print("results/"+resultName+"/"+TString(hNum->GetName())+".pdf","pdf");
	  
	}
	// end TEfficiency

      }
    }
  }
  // end loop over histograms


  // Write histograms //
  for(UInt_t iV=0 ; iV<nV ; iV++) { // x-axis variables
    for(UInt_t iF=0 ; iF<nF ; iF++) { // num/den
      for(UInt_t iP=0 ; iP<nP ; iP++) { // paths
	for(UInt_t iS=0 ; iS<nS ; iS++) { // steps in the paths
	  h[iV][iF][iP][iS]->Write();
	}
      }
    }
  }


  outfile->Write();


  // clean memory
  delete ch;

  // END //
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
