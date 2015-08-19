void maketemplate(TTree* tree, int nbins, double ptmin, double ptmax, double etamin, double etamax, const char* outfilename) {

    TFile* pufile = new TFile("purwt.root");
    TH1* puhist = (TH1*)pufile->Get("puhist");

    double xmin = 65.0;
    double xmax = 115.0;

    double wgtsum = 3.129e+11;

    TH1F hpass("hpass", "", nbins, xmin, xmax);
    TH1F hfail("hfail", "", nbins, xmin, xmax);
    TH1F hp("hp" , "", 1, xmin, xmax);
    TH1F ha("ha" , "", 1, xmin, xmax);
    TH1F hr("hr" , "", 1, xmin, xmax);

    hp.Sumw2();
    ha.Sumw2();

    TBranch  *bmass       = tree->GetBranch("mass");
    TBranch  *bpt         = tree->GetBranch("pt");
    TBranch  *beta        = tree->GetBranch("abseta");
    TBranch  *bmc         = tree->GetBranch("mcTrue");
    TBranch  *bid         = tree->GetBranch("tightid");
    TBranch  *bnvtx       = tree->GetBranch("nvtx");
    TBranch  *bwgt        = tree->GetBranch("wgt");

    Float_t  mass         = 0.0;
    Float_t  pt           = 0.0;
    Float_t  eta          = 0.0;
    Float_t  nvtx         = 0.0;
    Float_t  wgt          = 0.0;
    Int_t    mc           = 0;
    Int_t    id           = 0;

    bmass->SetAddress(&mass);
    bpt  ->SetAddress(&pt);
    beta ->SetAddress(&eta);
    bmc  ->SetAddress(&mc);
    bid  ->SetAddress(&id);
    bnvtx->SetAddress(&nvtx);
    bwgt ->SetAddress(&wgt);

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        bmass->GetEvent(i);
        bpt  ->GetEvent(i);
        beta ->GetEvent(i);
        bmc  ->GetEvent(i);
        bid  ->GetEvent(i);
        bnvtx->GetEvent(i);
        bwgt ->GetEvent(i);

        Float_t puwgt = 0.;
        if (nvtx <= 49) puwgt = puhist->GetBinContent(nvtx);

        if (pt > ptmin && pt <= ptmax && eta > etamin && eta <= etamax && mc > 0) {
            if (id == 0) hfail.Fill(mass, puwgt*wgt/wgtsum);
            if (id >  0) hpass.Fill(mass, puwgt*wgt/wgtsum);
            if (id >  0) hp.Fill(mass, puwgt*wgt/wgtsum);
            ha.Fill(mass, puwgt*wgt/wgtsum);
        }
    }

    for (int j = 1; j <= nbins; j++) {
        if (hpass.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
        if (hfail.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
    }

    hr.Divide(&hp, &ha, 1.0, 1.0, "B");

    hpass.Smooth();
    hfail.Smooth();

    hpass.Smooth();
    hfail.Smooth();

    RooRealVar m("mass", "", xmin, xmax);
    RooDataHist dhpass("dhpass", "", RooArgList(m), RooFit::Import(hpass), 0);
    RooDataHist dhfail("dhfail", "", RooArgList(m), RooFit::Import(hfail), 0);
    RooHistPdf signalPassMC("signalPassMC", "", RooArgSet(m), dhpass);
    RooHistPdf signalFailMC("signalFailMC", "", RooArgSet(m), dhfail);

    RooWorkspace w("w", "");
    w.import(signalPassMC);
    w.import(signalFailMC);
        
    TFile outfile(outfilename, "RECREATE");
    w.Write();
    hpass.Write();
    hfail.Write();

    cout << "Efficiency -- pT [" << ptmin << ", " << ptmax << "], eta [" << etamin << ", " << etamax << "]   :   " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << endl;
}

void tnpsignaltemplates() {

    TFile* file = new TFile("tnptree.root");
    TTree* tree = (TTree*)file->Get("tnptree/fitter_tree");

    maketemplate(tree, 100, 10., 20., 1.2, 2.4, "templates0.root");    
    maketemplate(tree, 100, 20., 30., 1.2, 2.4, "templates1.root");    
    maketemplate(tree, 100, 30., 40., 1.2, 2.4, "templates2.root");    
    maketemplate(tree, 100, 40., 50., 1.2, 2.4, "templates3.root");    
    maketemplate(tree, 100, 50., 70., 1.2, 2.4, "templates4.root");    
    maketemplate(tree, 100, 70.,100., 1.2, 2.4, "templates5.root");    

}
