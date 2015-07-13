#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <TTree.h>

class LHEWeightsTreeMaker : public edm::EDAnalyzer {
    public:
        explicit LHEWeightsTreeMaker(const edm::ParameterSet&);
        ~LHEWeightsTreeMaker();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        // InputTags
        edm::InputTag lheInfoTag;
        edm::InputTag genInfoTag;

        // Tokens
        edm::EDGetTokenT<LHEEventProduct> lheInfoToken;
        edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;

        bool uselheweights, addqcdpdfweights;

        uint32_t event, run, lumi;
        double   wgtsign, wgtxsec, wgtpdf1, wgtpdf2, wgtpdf3, wgtpdf4, wgtpdf5;
        double*  wgtpdf;
        double*  wgtqcd;
        TTree* tree;
};

LHEWeightsTreeMaker::LHEWeightsTreeMaker(const edm::ParameterSet& iConfig): 
    lheInfoTag(iConfig.getParameter<edm::InputTag>("lheinfo")),
    genInfoTag(iConfig.getParameter<edm::InputTag>("geninfo")),
    uselheweights(iConfig.getParameter<bool>("uselheweights")),
    addqcdpdfweights(iConfig.getParameter<bool>("addqcdpdfweights"))
{
    // Token consumes instructions
    lheInfoToken = consumes<LHEEventProduct>(lheInfoTag);
    genInfoToken = consumes<GenEventInfoProduct>(genInfoTag);

    wgtqcd = new double[8];
    wgtpdf = new double[100];
}


LHEWeightsTreeMaker::~LHEWeightsTreeMaker() {
    delete wgtqcd;
    delete wgtpdf;
}

void LHEWeightsTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    // Get handles to all the requisite collections
    Handle<LHEEventProduct> lheInfoH;
    if (uselheweights) iEvent.getByToken(lheInfoToken, lheInfoH);

    Handle<GenEventInfoProduct> genInfoH;
    if (uselheweights) iEvent.getByToken(genInfoToken, genInfoH);

    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();

    // Weights info
    wgtsign = 1.0;
    wgtxsec = 1.0;
    if (uselheweights) {
        wgtsign = genInfoH->weight();
        wgtxsec = lheInfoH->originalXWGTUP();
    }
    wgtpdf1 = 0.0;
    wgtpdf2 = 0.0;
    wgtpdf3 = 0.0;
    wgtpdf4 = 0.0;
    wgtpdf5 = 0.0;
    
    for (size_t i = 0; i < 8  ; i++) wgtqcd[i] = 0.;
    for (size_t i = 0; i < 100; i++) wgtpdf[i] = 0.;

    if (addqcdpdfweights) {
        vector<gen::WeightsInfo> weights = lheInfoH->weights();
        for (size_t i = 0; i < weights.size(); i++) {
            if (weights[i].id == "315") wgtpdf1 = weights[i].wgt; // cteq6l1
            if (weights[i].id == "316") wgtpdf2 = weights[i].wgt; // MMHT2014lo68cl
            if (weights[i].id == "370") wgtpdf3 = weights[i].wgt; // HERAPDF15LO
            if (weights[i].id == "393") wgtpdf4 = weights[i].wgt; // CT10nlo
            if (weights[i].id == "446") wgtpdf5 = weights[i].wgt; // MMHT2014nlo68cl
        
            for (size_t j = 2; j <= 9; j++) {
                stringstream ss;
                ss << j;
                if (weights[i].id == ss.str()) wgtqcd[j-2]  = weights[i].wgt;
            }
            for (size_t j = 11; j <= 110; j++) {
                stringstream ss;
                ss << j;
                if (weights[i].id == ss.str()) wgtpdf[j-11] = weights[i].wgt;
            }
        }
    }

    tree->Fill();

}


void LHEWeightsTreeMaker::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("gentree"    , "gentree");
    // Run, Lumi, Event info
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
    // Event weights
    tree->Branch("wgtsign"              , &wgtsign              , "wgtsign/D");
    tree->Branch("wgtxsec"              , &wgtxsec              , "wgtxsec/D");
    if (addqcdpdfweights) {
    tree->Branch("wgtpdf1"              , &wgtpdf1              , "wgtpdf1/D");
    tree->Branch("wgtpdf2"              , &wgtpdf2              , "wgtpdf2/D");
    tree->Branch("wgtpdf3"              , &wgtpdf3              , "wgtpdf3/D");
    tree->Branch("wgtpdf4"              , &wgtpdf4              , "wgtpdf4/D");
    tree->Branch("wgtpdf5"              , &wgtpdf5              , "wgtpdf5/D");
    tree->Branch("wgtpdf"               ,  wgtpdf               , "wgtpdf[100]/D");
    tree->Branch("wgtqcd"               ,  wgtqcd               , "wgtqcd[8]/D");
    }
}

void LHEWeightsTreeMaker::endJob() {
}

void LHEWeightsTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void LHEWeightsTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void LHEWeightsTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void LHEWeightsTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void LHEWeightsTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LHEWeightsTreeMaker);

