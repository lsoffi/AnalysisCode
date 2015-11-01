#include "/user/ndaci/WorkArea/DarkMatter/CMSSW_7413_Update/src/MonoXAnalysis/AnalysisStep/macros/myIncludes.h"

double sumwgt(TTree* tree);

void createWgtSumAndPU() {

    TFile* pufile = new TFile("/user/ndaci/Data/XMET/purwt.root");
    TH1*   puhist = (TH1*)pufile->Get("puhist");

    TFile* infile = new TFile("tree.root");
    TTree* intree = (TTree*)infile->Get("gentree/gentree");
    TTree* frtree = (TTree*)infile->Get("tree/tree");

    double wgtsum = sumwgt(intree);
    double puweight=1;
    
    //const char* cut = "nmuons == 2 && zmass > 60 && zmass < 120 && mu1pt > 20 && mu1id == 1";
    const char* cut = "t1mumet>50";

    TFile* outfile = new TFile("skim.root", "RECREATE");
    outfile->cd();
    TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
    treedir->cd();
    TTree* outtree = frtree->CopyTree(cut);

    TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    TBranch* bpuwgt  = outtree->Branch("puweight", &puweight, "puweight/D");

    uint32_t nvtx;
    outtree->SetBranchAddress("nvtx", &nvtx);
    
    Long64_t nEntries=outtree->GetEntries();
    cout << "nEntries=" << nEntries << endl;

    for (Long64_t i = 0; i < nEntries; i++) {
      bwgtsum->Fill();

      outtree->GetEntry(i);
      puweight = 0;
      if(nvtx<=35) puweight = puhist->GetBinContent(nvtx);
      bpuwgt->Fill();
      if(i%1000) cout << "i=" << i << " nvtx=" << nvtx << " puweight=" << puweight << endl;
    }

    outfile->Write();

    gSystem->Exit(0);
}

double sumwgt(TTree* tree) {
    TBranch *bweight = tree->GetBranch("wgtsign");

    double vweight  = 0.0;

    bweight->SetAddress(&vweight);

    double weightsum = 0.;
    for (int i = 0; i < tree->GetEntries(); i++) {
        bweight->GetEvent(i);
        weightsum += vweight;
    }

    return weightsum;
}
