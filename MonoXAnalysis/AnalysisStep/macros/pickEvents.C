#include "myIncludes.h"

Int_t pickEvents(TString source)
{

  ofstream outlog("check_1jet_lowdphi_"+source+".txt",ios::out);
  TString dir="/user/ndaci/Data/XMET/AdishTrees_15Oct2015/";

  TChain* tree = new TChain("tree/tree");
  if(source=="data_met") {
    tree->Add(dir+"/Run2015D/met/skimMumet100WgtSum.root");
  }
  else if(source=="qcd") {
    /*
    tree->Add(dir+"/Spring15MC_25ns/qcdht1000to1500/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht100to200/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht1500to2000/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht2000toinf/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht200to300/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht300to500/skimMumet100WgtSum.root");
    tree->Add(dir+"/Spring15MC_25ns/qcdht500to700/skimMumet100WgtSum.root");
    */
    tree->Add(dir+"/Spring15MC_25ns/qcdht700to1000/skimMumet100WgtSum.root");
  }
    

  double   jetmetdphimin, incjetmetdphimin, t1mumet, signaljetpt, signaljeteta, signaljetphi, signaljetbtag, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi, signaljetqgl, signaljetqgs2, signaljetqgptd;
  uint32_t event, run, lumi;
  uint32_t nvtx, nmuons, nelectrons, ntaus, ntightmuons, ntightelectrons, nphotons, njets, njets80, nbjets, nfatjets, nfwdjets, nfwdjets80, nsoftjets, nsoftbjets, nsoftfwdjets;
  uint8_t  hltmet90, hltmet120, hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel, hltpfht200, hltpfht250, hltpfht300, hltpfht350, hltpfht400, hltpfht475, hltpfht600, hltpfht650, hltpfht800;

  tree->SetBranchAddress("event"                , &event                ); // , "event/i");
  tree->SetBranchAddress("run"                  , &run                  ); // , "run/i");
  tree->SetBranchAddress("lumi"                 , &lumi                 ); // , "lumi/i");
  tree->SetBranchAddress("hltmet90"             , &hltmet90             ); // , "hltmet90/b");
  tree->SetBranchAddress("hltmet120"            , &hltmet120            ); // , "hltmet120/b");
  tree->SetBranchAddress("hltmetwithmu90"       , &hltmetwithmu90       ); // , "hltmetwithmu90/b");
  tree->SetBranchAddress("hltmetwithmu120"      , &hltmetwithmu120      ); // , "hltmetwithmu120/b");
  tree->SetBranchAddress("hltmetwithmu170"      , &hltmetwithmu170      ); // , "hltmetwithmu170/b");
  tree->SetBranchAddress("hltmetwithmu300"      , &hltmetwithmu300      ); // , "hltmetwithmu300/b");
  tree->SetBranchAddress("hltjetmet90"          , &hltjetmet90          ); // , "hltjetmet90/b");
  tree->SetBranchAddress("hltjetmet120"         , &hltjetmet120         ); // , "hltjetmet120/b");
  tree->SetBranchAddress("hltphoton165"         , &hltphoton165         ); // , "hltphoton165/b");
  tree->SetBranchAddress("hltphoton175"         , &hltphoton175         ); // , "hltphoton175/b");
  tree->SetBranchAddress("hltdoublemu"          , &hltdoublemu          ); // , "hltdoublemu/b");
  tree->SetBranchAddress("hltsinglemu"          , &hltsinglemu          ); // , "hltsinglemu/b");
  tree->SetBranchAddress("hltdoubleel"          , &hltdoubleel          ); // , "hltdoubleel/b");
  tree->SetBranchAddress("hltsingleel"          , &hltsingleel          ); // , "hltsingleel/b");
  tree->SetBranchAddress("t1mumet"              , &t1mumet              ); // , "t1mumet/D");
  tree->SetBranchAddress("nmuons"               , &nmuons               ); // , "nmuons/i");
  tree->SetBranchAddress("nelectrons"           , &nelectrons           ); // , "nelectrons/i");
  tree->SetBranchAddress("ntightmuons"          , &ntightmuons          ); // , "ntightmuons/i");
  tree->SetBranchAddress("ntightelectrons"      , &ntightelectrons      ); // , "ntightelectrons/i");
  tree->SetBranchAddress("ntaus"                , &ntaus                ); // , "ntaus/i");
  tree->SetBranchAddress("njets"                , &njets                ); // , "njets/i");
  tree->SetBranchAddress("nbjets"               , &nbjets               ); // , "nbjets/i");
  tree->SetBranchAddress("nfatjets"             , &nfatjets             ); // , "nfatjets/i");
  tree->SetBranchAddress("nphotons"             , &nphotons             ); // , "nphotons/i");
  tree->SetBranchAddress("signaljetpt"          , &signaljetpt          ); // , "signaljetpt/D");
  tree->SetBranchAddress("signaljeteta"         , &signaljeteta         ); // , "signaljeteta/D");
  tree->SetBranchAddress("signaljetCHfrac"      , &signaljetCHfrac      ); // , "signaljetCHfrac/D");
  tree->SetBranchAddress("signaljetNHfrac"      , &signaljetNHfrac      ); // , "signaljetNHfrac/D");
  tree->SetBranchAddress("signaljetEMfrac"      , &signaljetEMfrac      ); // , "signaljetEMfrac/D");
  tree->SetBranchAddress("signaljetCEMfrac"     , &signaljetCEMfrac     ); // , "signaljetCEMfrac/D");
  tree->SetBranchAddress("jetmetdphimin"        , &jetmetdphimin        ); // , "jetmetdphimin/D");
  tree->SetBranchAddress("incjetmetdphimin"     , &incjetmetdphimin     ); // , "incjetmetdphimin/D");

  TCut SR_1jet="(hltmet90>0 || hltmet120>0) && t1mumet>250 && nelectrons==0 && ntaus==0 && nmuons==0 && nphotons==0 && signaljetpt>30 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && njets==1";

  //UInt_t nEntries = tree->GetEntries();

  //tree->SetEntryList(0);
  tree->Draw(">>+skim", SR_1jet, "entrylist");
  TEntryList* skim = (TEntryList*)gDirectory->Get("skim");
  tree->SetEntryList(skim);

  Int_t nEntries = skim->GetN();
  Int_t iEntry, chainEntry, treenum;

  for(Int_t iE=0 ; iE<nEntries ; iE++) {

    //tree->GetEntry(iE);

    iEntry     = skim->GetEntryAndTree(iE, treenum);
    chainEntry = iEntry + (tree->GetTreeOffset())[treenum];
    tree->GetEntry(chainEntry);

    if(iE%1000==0) cout << "Processed entry #" << iE << "/" << nEntries << endl;

    if(TMath::Abs(jetmetdphimin)<0.1) {
      outlog << run << ":" << lumi << ":" << event << endl;

      cout << tree->GetFile()->GetName() << endl;
    }

  }

  return 0;
}
