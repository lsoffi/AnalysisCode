#include <memory>
#include <vector>
#include <iostream>

#include <TRandom3.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

class PFCleaner : public edm::EDProducer {
    public:
        explicit PFCleaner(const edm::ParameterSet&);
        ~PFCleaner();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        bool randomConeOverlaps(double, double, double, std::vector<pat::Jet>);

        edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;
        edm::EDGetTokenT<edm::View<reco::Candidate> > pfcandsToken;
        edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
        edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken;
        edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken;
        edm::EDGetTokenT<std::vector<pat::Photon> > photonsToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > photonLooseIdMapToken;

        bool userandomphi;
};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 
    verticesToken            (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
    pfcandsToken             (consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("pfcands"))),
    jetsToken                (consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
    muonsToken               (consumes<std::vector<pat::Muon> > (iConfig.getParameter<edm::InputTag>("muons"))), 
    electronsToken           (consumes<std::vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons"))),
    photonsToken             (consumes<std::vector<pat::Photon> > (iConfig.getParameter<edm::InputTag>("photons"))),
    electronVetoIdMapToken   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronidveto"))),
    electronMediumIdMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronidmedium"))),
    photonLooseIdMapToken    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonidloose"))),
    userandomphi             (iConfig.existsAs<bool>("userandomphiforRC") ? iConfig.getParameter<bool>("userandomphiforRC") : false)
{
    produces<pat::MuonRefVector>("muons");
    produces<pat::ElectronRefVector>("electrons");
    produces<pat::PhotonRefVector>("photons");
    produces<pat::MuonRefVector>("tightmuons");
    produces<pat::ElectronRefVector>("tightelectrons");
    produces<pat::PhotonRefVector>("tightphotons");
    produces<pat::PhotonRefVector>("loosephotons");
    produces<edm::ValueMap<float> >("rndgammaiso");
    produces<edm::ValueMap<float> >("rndchhadiso");
}


PFCleaner::~PFCleaner() {
}

void PFCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<std::vector<reco::Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);

    Handle<edm::View<reco::Candidate> > pfcandsH;
    iEvent.getByToken(pfcandsToken, pfcandsH);

    Handle<std::vector<pat::Jet> > jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<std::vector<pat::Muon> > muonsH;
    iEvent.getByToken(muonsToken, muonsH);

    Handle<std::vector<pat::Electron> > electronsH;
    iEvent.getByToken(electronsToken, electronsH);

    Handle<std::vector<pat::Photon> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    Handle<edm::ValueMap<bool> > electronVetoIdH;
    iEvent.getByToken(electronVetoIdMapToken, electronVetoIdH);

    Handle<edm::ValueMap<bool> > electronMediumIdH;
    iEvent.getByToken(electronMediumIdMapToken, electronMediumIdH);

    Handle<edm::ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdMapToken, photonLooseIdH);

    std::auto_ptr<pat::MuonRefVector> outputmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputphotons(new pat::PhotonRefVector);
    std::auto_ptr<pat::MuonRefVector> outputtightmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputtightelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputtightphotons(new pat::PhotonRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputloosephotons(new pat::PhotonRefVector);
    std::auto_ptr<edm::ValueMap<float> > outputgammaisomap(new ValueMap<float>());
    std::auto_ptr<edm::ValueMap<float> > outputchhadisomap(new ValueMap<float>());

    TRandom3 rand;

    for (vector<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        if (verticesH->size() == 0) continue;
        bool passeskincuts = (muons_iter->pt() > 10 && fabs(muons_iter->eta()) < 2.4);
        float isoval = muons_iter->pfIsolationR04().sumNeutralHadronEt;
        isoval += muons_iter->pfIsolationR04().sumPhotonEt;
        isoval -= 0.5*muons_iter->pfIsolationR04().sumPUPt;
        if (isoval < 0.) isoval = 0.;
        isoval += muons_iter->pfIsolationR04().sumChargedHadronPt;
        isoval /= muons_iter->pt();

        if (passeskincuts) {
            if (muon::isLooseMuon(*muons_iter) && isoval < 0.2) outputmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
            if (muon::isTightMuon(*muons_iter, *(verticesH->begin())) && isoval < 0.12) outputtightmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
        }
    }

    for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        const Ptr<pat::Electron> electronPtr(electronsH, electrons_iter - electronsH->begin());
        bool passeskincuts = (electrons_iter->pt() > 10 && fabs(electrons_iter->superCluster()->eta()) < 2.5);
        bool passesvetoid = (*electronVetoIdH)[electronPtr];
        bool passesmediumid = (*electronMediumIdH)[electronPtr];

        if (passeskincuts && passesvetoid) outputelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
        if (passeskincuts && passesmediumid) outputtightelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }

    std::vector<float> rndgammaiso;
    std::vector<float> rndchhadiso;


    for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
        float gaisoval = 0.;
        float chisoval = 0.;
        double rndphi = photons_iter->phi() + M_PI/2.0;
        if (userandomphi) {
            rndphi = rand.Uniform(-M_PI, M_PI);
            while (randomConeOverlaps(rndphi, photons_iter->eta(), photons_iter->phi(), *jetsH)) rndphi = rand.Uniform(-M_PI, M_PI);
        }

        for(size_t i = 0; i < pfcandsH->size(); i++) {
            const auto& pfcand = pfcandsH->ptrAt(i);
            if (    pfcand->pdgId()  ==  22 && deltaR(photons_iter->eta(), rndphi, pfcand->eta(), pfcand->phi()) <= 0.3) gaisoval += pfcand->pt();
            if (abs(pfcand->pdgId()) == 211 && deltaR(photons_iter->eta(), rndphi, pfcand->eta(), pfcand->phi()) <= 0.3) chisoval += pfcand->pt();
        }
        rndgammaiso.push_back(gaisoval);
        rndchhadiso.push_back(chisoval);

        if (fabs(photons_iter->superCluster()->eta()) > 2.5 || photons_iter->pt() < 15) continue;
        if (photons_iter->r9() > 0.8 || photons_iter->chargedHadronIso() < 20. || photons_iter->chargedHadronIso() < photons_iter->pt()*0.3) outputloosephotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin())); 

        const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
        bool passeslooseid = (*photonLooseIdH)[photonPtr];
        if (passeslooseid) {
            outputphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
            if (photons_iter->pt() > 175) outputtightphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
        }
    }

    edm::ValueMap<float>::Filler gafiller(*outputgammaisomap);
    gafiller.insert(photonsH, rndgammaiso.begin(), rndgammaiso.end());
    gafiller.fill();

    edm::ValueMap<float>::Filler chfiller(*outputchhadisomap);
    chfiller.insert(photonsH, rndchhadiso.begin(), rndchhadiso.end());
    chfiller.fill();

    iEvent.put(outputmuons, "muons");
    iEvent.put(outputelectrons, "electrons");
    iEvent.put(outputphotons, "photons");
    iEvent.put(outputtightmuons, "tightmuons");
    iEvent.put(outputtightelectrons, "tightelectrons");
    iEvent.put(outputtightphotons, "tightphotons");
    iEvent.put(outputloosephotons, "loosephotons");
    iEvent.put(outputgammaisomap, "rndgammaiso");
    iEvent.put(outputchhadisomap, "rndchhadiso");
}

void PFCleaner::beginJob() {
}

void PFCleaner::endJob() {
}

void PFCleaner::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void PFCleaner::endRun(edm::Run const&, edm::EventSetup const&) {
}

void PFCleaner::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void PFCleaner::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void PFCleaner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

bool PFCleaner::randomConeOverlaps(double randomphi, double photoneta, double photonphi, std::vector<pat::Jet> jets) {
    if (reco::deltaR(photoneta, randomphi, photoneta, photonphi) < 0.4) return true;
    for (std::size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > 30. && reco::deltaR(photoneta, randomphi, jets[i].eta(), jets[i].phi()) < 0.4) return true;
    }
    return false;
}

DEFINE_FWK_MODULE(PFCleaner);
