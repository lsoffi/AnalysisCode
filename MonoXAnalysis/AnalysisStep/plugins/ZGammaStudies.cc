#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <TH1F.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

class ZGammaStudies : public edm::EDAnalyzer {
    public:
        explicit ZGammaStudies(const edm::ParameterSet&);
        ~ZGammaStudies();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        edm::InputTag gensTag;
        edm::InputTag genjetsTag;
        edm::InputTag triggerResultsTag;
        edm::InputTag muonsTag;
        edm::InputTag electronsTag;
        edm::InputTag photonsTag;
        edm::InputTag tausTag;
        edm::InputTag jetsTag;
        edm::InputTag t1pfmetTag;
        bool isPhotonSample;
        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        TTree* tree;

        int32_t  puobs, putrue; 
        int32_t  vid, l1id, l2id; 
        uint32_t event, run, lumi;
        uint32_t hltmet90, hltmet120, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175;
        uint32_t nmuons, nelectrons, ntaus, nphotons, njets, ngenjets;
        double   t1pfmet, t1pfmetphi, t1phmet, t1phmetphi;
        double   phpt, pheta, phphi;
        double   signaljetpt, signaljeteta, signaljetphi, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi;
        double   secondjetpt, secondjeteta, secondjetphi, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetCEMfrac, secondjetmetdphi;
        double   jetjetdphi;
        double   signalgenjetpt, signalgenjeteta, signalgenjetphi;
        double   secondgenjetpt, secondgenjeteta, secondgenjetphi;
        double   genjetjetdphi;
        double   vmass, vmt, vpt, veta, vphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
        double   wgt;

        struct PatJetPtSorter {
            bool operator() (const pat::Jet& i, const pat::Jet& j) {
                return (i.pt() > j.pt());
            }
        } patjetsorter;

        struct GenJetPtSorter {
            bool operator() (const reco::GenJet& i, const reco::GenJet& j) {
                return (i.pt() > j.pt());
            }
        } genjetsorter;


};

ZGammaStudies::ZGammaStudies(const edm::ParameterSet& iConfig): 
    gensTag(iConfig.getParameter<edm::InputTag>("gens")),
    genjetsTag(iConfig.getParameter<edm::InputTag>("genjets")),
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
    electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
    photonsTag(iConfig.getParameter<edm::InputTag>("photons")),
    tausTag(iConfig.getParameter<edm::InputTag>("taus")),
    jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
    t1pfmetTag(iConfig.getParameter<edm::InputTag>("t1pfmet")),
    isPhotonSample(iConfig.getParameter<bool>("isPhotonSample")),
    wgt(iConfig.getParameter<double>("weight"))
{
}


ZGammaStudies::~ZGammaStudies() {
}

void ZGammaStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get handles to all the requisite collections
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByLabel(triggerResultsTag, triggerResultsH);

    Handle<View<GenParticle> > gensH;
    iEvent.getByLabel(gensTag, gensH);

    Handle<View<GenJet> > genjetsH;
    iEvent.getByLabel(genjetsTag, genjetsH);

    Handle<pat::MuonRefVector> muonsH;
    iEvent.getByLabel(muonsTag, muonsH);
    pat::MuonRefVector muons = *muonsH;

    Handle<pat::ElectronRefVector> electronsH;
    iEvent.getByLabel(electronsTag, electronsH);
    pat::ElectronRefVector electrons = *electronsH;

    Handle<pat::PhotonRefVector> photonsH;
    iEvent.getByLabel(photonsTag, photonsH);
    pat::PhotonRefVector photons = *photonsH;

    Handle<View<pat::Tau> > tausH;
    iEvent.getByLabel(tausTag, tausH);

    Handle<View<pat::Jet> > jetsH;
    iEvent.getByLabel(jetsTag, jetsH);

    Handle<View<MET> > t1pfmetH;
    iEvent.getByLabel(t1pfmetTag, t1pfmetH);

    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();

    // Trigger info
    hltmet90     = 0;
    hltmet120    = 0;
    hltjetmet90  = 0;
    hltjetmet120 = 0;
    hltphoton165 = 0;
    hltphoton175 = 0;

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
        if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90     = 1; // MET trigger
        if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120    = 1; // MET trigger
        if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90  = 1; // Jet-MET trigger
        if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120 = 1; // Jet-MET trigger
        if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165 = 1; // Photon trigger
        if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175 = 1; // Photon trigger
    }

    // Generator-level information
    vid    = 0;
    vmass  = 0.0;
    vmt    = 0.0;
    vpt    = 0.0;
    veta   = 0.0;
    vphi   = 0.0;
    l1id   = 0;
    l1pt   = 0.0;
    l1eta  = 0.0;
    l1phi  = 0.0;
    l2id   = 0;
    l2pt   = 0.0;
    l2eta  = 0.0;
    l2phi  = 0.0;

    if (!isPhotonSample) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() > 1 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) {
                vid    = gens_iter->pdgId();
                vmass  = gens_iter->mass();
                vpt    = gens_iter->pt();
                veta   = gens_iter->eta();
                vphi   = gens_iter->phi();
                l1id   = gens_iter->daughter(0)->pdgId();
                l1pt   = gens_iter->daughter(0)->pt();
                l1eta  = gens_iter->daughter(0)->eta();
                l1phi  = gens_iter->daughter(0)->phi();
                l2id   = gens_iter->daughter(1)->pdgId();
                l2pt   = gens_iter->daughter(1)->pt();
                l2eta  = gens_iter->daughter(1)->eta();
                l2phi  = gens_iter->daughter(1)->phi();
                vmt    = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi)))); 
            }
        }
    }

    else {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if (gens_iter->pdgId() == 22 && gens_iter->status() == 1 && gens_iter->pt() > vpt) {
                vid    = gens_iter->pdgId();
                vpt    = gens_iter->pt();
                veta   = gens_iter->eta();
                vphi   = gens_iter->phi();
            }
        }
    }

    vector<GenJet> genjets;
    for(View<GenJet>::const_iterator genjets_iter = genjetsH->begin(); genjets_iter != genjetsH->end(); ++genjets_iter) {
        bool skipjet = false;
        if (fabs(genjets_iter->eta()) > 2.5) skipjet = true;
        if (isPhotonSample  &&  vpt > 0. && deltaR( veta,  vphi, genjets_iter->eta(), genjets_iter->phi()) < 0.4) skipjet = true;
        //if (!isPhotonSample && l1pt > 0. && deltaR(l1eta, l1phi, genjets_iter->eta(), genjets_iter->phi()) < 0.4) skipjet = true;
        //if (!isPhotonSample && l2pt > 0. && deltaR(l2eta, l2phi, genjets_iter->eta(), genjets_iter->phi()) < 0.4) skipjet = true;
        if (!isPhotonSample) {
            for (size_t i = 0; i < genjets_iter->numberOfDaughters(); i++) {
                if (!genjets_iter->daughterPtr(i).isAvailable()) continue;
                if (genjets_iter->daughter(i)->pdgId() == l1id && genjets_iter->daughter(i)->pt() == l1pt && genjets_iter->daughter(i)->eta() == l1eta && genjets_iter->daughter(i)->phi() == l1phi) skipjet = true;
                if (genjets_iter->daughter(i)->pdgId() == l2id && genjets_iter->daughter(i)->pt() == l2pt && genjets_iter->daughter(i)->eta() == l2eta && genjets_iter->daughter(i)->phi() == l2phi) skipjet = true;
            }
        }            
        GenJet genjet = *genjets_iter;
        if (!skipjet) genjets.push_back(genjet);
    }
    sort(genjets.begin(), genjets.end(), genjetsorter);

    ngenjets = 0;
    for (size_t i = 0; i < genjets.size(); i++) {
        if (genjets[i].pt() > 30) ngenjets++;
    }

    if (genjets.size() > 0) {
        signalgenjetpt      = genjets[0].pt();
        signalgenjeteta     = genjets[0].eta();
        signalgenjetphi     = genjets[0].phi();
    }

    if (genjets.size() > 1) {
        secondgenjetpt      = genjets[1].pt();
        secondgenjeteta     = genjets[1].eta();
        secondgenjetphi     = genjets[1].phi();
    }

    if ( genjets.size() > 1) genjetjetdphi = deltaPhi(signalgenjetphi, secondgenjetphi);

    // MET information 
    t1pfmet      = t1pfmetH->front().et();
    t1pfmetphi   = t1pfmetH->front().phi();
  
    t1phmet      = t1pfmetH->front().et();
    t1phmetphi   = t1pfmetH->front().phi();

    // Jet information
    int hardestPhotonIndex = -1;
    double hardestPhotonPt = 0.0;
    for (size_t i = 0; i < photons.size(); i++) {
        if (photons[i]->pt() > hardestPhotonPt) {
            hardestPhotonIndex = i;
            hardestPhotonPt = photons[i]->pt();
        }
    }

    if (hardestPhotonIndex >= 0) {
        double t1phmetx = 0.;
        double t1phmety = 0.;

        phpt       = photons[hardestPhotonIndex]->pt();
        phphi      = photons[hardestPhotonIndex]->phi();

        t1phmetx   = t1phmet * cos(t1phmetphi);
        t1phmety   = t1phmet * sin(t1phmetphi);

        t1phmetx  += phpt * cos(phphi);
        t1phmety  += phpt * sin(phphi);

        t1phmet    = t1phmetx*t1phmetx + t1phmety*t1phmety;
        t1phmetphi = atan2(t1phmety, t1phmetx);
    } 

    pat::JetCollection jets;
    for (View<pat::Jet>::const_iterator jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
        if (fabs(jets_iter->eta()) > 2.5) continue;
        bool skipjet = false;
        if (isPhotonSample && vpt > 0. && deltaR(veta, vphi, jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < photons.size(); j++) {
            if (deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        if (skipjet) continue;
        bool passjetid = false;
        if (jets_iter->neutralHadronEnergyFraction() < 0.99 && jets_iter->neutralEmEnergyFraction() < 0.99 && (jets_iter->chargedMultiplicity() + jets_iter->neutralMultiplicity()) > 1 && jets_iter->muonEnergyFraction() < 0.8) {
            if (fabs(jets_iter->eta()) > 2.4) passjetid = true;
            else if (fabs(jets_iter->eta()) <= 2.4 && jets_iter->chargedHadronEnergyFraction() > 0. && jets_iter->chargedEmEnergyFraction() < 0.99 && jets_iter->chargedMultiplicity() > 0) passjetid = true;
        }
        if (!passjetid) continue;
        bool passpuid = false;
        double puidval = jets_iter->userFloat("pileupJetId:fullDiscriminant");
        double jetabseta = fabs(jets_iter->eta());
        if (jetabseta >= 0.00 && jetabseta < 2.50 && puidval > -0.63) passpuid = true;
        if (jetabseta >= 2.50 && jetabseta < 2.75 && puidval > -0.60) passpuid = true;
        if (jetabseta >= 2.75 && jetabseta < 3.00 && puidval > -0.55) passpuid = true;
        if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.45) passpuid = true;
        if (!passpuid) continue;
        pat::Jet jet = *jets_iter;
        jets.push_back(jet);
    }
    sort(jets.begin(), jets.end(), patjetsorter);

    njets     = 0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > 30) njets++;
    }

    signaljetpt       = 0.0;
    signaljeteta      = 0.0;
    signaljetphi      = 0.0;
    signaljetCHfrac   = 0.0;
    signaljetNHfrac   = 0.0;
    signaljetEMfrac   = 0.0;
    signaljetCEMfrac  = 0.0;
    signaljetmetdphi  = 0.0;
    secondjetpt       = 0.0;
    secondjeteta      = 0.0;
    secondjetphi      = 0.0;
    secondjetCHfrac   = 0.0;
    secondjetNHfrac   = 0.0;
    secondjetEMfrac   = 0.0;
    secondjetCEMfrac  = 0.0;
    jetjetdphi        = 0.0;

    if (jets.size() > 0) {
        signaljetpt      = jets[0].pt();
        signaljeteta     = jets[0].eta();
        signaljetphi     = jets[0].phi();
        signaljetCHfrac  = jets[0].chargedHadronEnergyFraction();
        signaljetNHfrac  = jets[0].neutralHadronEnergyFraction();
        signaljetEMfrac  = jets[0].neutralEmEnergyFraction();
        signaljetCEMfrac = jets[0].chargedEmEnergyFraction();
    }

    if (jets.size() > 1) {
        secondjetpt      = jets[1].pt();
        secondjeteta     = jets[1].eta();
        secondjetphi     = jets[1].phi();
        secondjetCHfrac  = jets[1].chargedHadronEnergyFraction();
        secondjetNHfrac  = jets[1].neutralHadronEnergyFraction();
        secondjetEMfrac  = jets[1].neutralEmEnergyFraction();
        secondjetCEMfrac = jets[1].chargedEmEnergyFraction();
    }

    if ( jets.size() > 1) jetjetdphi        = deltaPhi(signaljetphi, secondjetphi);

    // Lepton counts
    nmuons          = muonsH->size();
    nelectrons      = electronsH->size();
    nphotons        = photonsH->size();
    ntaus           = 0;
    for (View<pat::Tau>::const_iterator taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
        if (taus_iter->pt() > 18 && fabs(taus_iter->eta()) < 2.3 && taus_iter->tauID("decayModeFinding") > 0.5 && taus_iter->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 5) ntaus++;
    }

    tree->Fill();
}


void ZGammaStudies::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"       , "tree");
    // Run, Lumi, Event info
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
    // Event weight
    tree->Branch("wgt"                  , &wgt                  , "wgt/D");
    // Triggers
    tree->Branch("hltmet90"             , &hltmet90             , "hltmet90/i");
    tree->Branch("hltmet120"            , &hltmet120            , "hltmet120/i");
    tree->Branch("hltjetmet90"          , &hltjetmet90          , "hltjetmet90/i");
    tree->Branch("hltjetmet120"         , &hltjetmet120         , "hltjetmet120/i");
    tree->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/i");
    tree->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/i");
    // Object counts
    tree->Branch("nmuons"               , &nmuons               , "nmuons/i");
    tree->Branch("nelectrons"           , &nelectrons           , "nelectrons/i");
    tree->Branch("nphotons"             , &nphotons             , "nphotons/i");
    tree->Branch("ntaus"                , &ntaus                , "ntaus/i");
    tree->Branch("njets"                , &njets                , "njets/i");
    tree->Branch("ngenjets"             , &ngenjets             , "ngenjets/i");
    // MET info
    tree->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/D");
    tree->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/D");
    tree->Branch("t1phmet"              , &t1phmet              , "t1phmet/D");
    tree->Branch("t1phmetphi"           , &t1phmetphi           , "t1phmetphi/D");
    // Jet info
    tree->Branch("signalgenjetpt"       , &signalgenjetpt       , "signalgenjetpt/D");
    tree->Branch("signalgenjeteta"      , &signalgenjeteta      , "signalgenjeteta/D");
    tree->Branch("signalgenjetphi"      , &signalgenjetphi      , "signalgenjetphi/D");
    tree->Branch("secondgenjetpt"       , &secondgenjetpt       , "secondgenjetpt/D");
    tree->Branch("secondgenjeteta"      , &secondgenjeteta      , "secondgenjeteta/D");
    tree->Branch("secondgenjetphi"      , &secondgenjetphi      , "secondgenjetphi/D");
    tree->Branch("signaljetpt"          , &signaljetpt          , "signaljetpt/D");
    tree->Branch("signaljeteta"         , &signaljeteta         , "signaljeteta/D");
    tree->Branch("signaljetphi"         , &signaljetphi         , "signaljetphi/D");
    tree->Branch("signaljetCHfrac"      , &signaljetCHfrac      , "signaljetCHfrac/D");
    tree->Branch("signaljetNHfrac"      , &signaljetNHfrac      , "signaljetNHfrac/D");
    tree->Branch("signaljetEMfrac"      , &signaljetEMfrac      , "signaljetEMfrac/D");
    tree->Branch("signaljetCEMfrac"     , &signaljetCEMfrac     , "signaljetCEMfrac/D");
    tree->Branch("secondjetpt"          , &secondjetpt          , "secondjetpt/D");
    tree->Branch("secondjeteta"         , &secondjeteta         , "secondjeteta/D");
    tree->Branch("secondjetphi"         , &secondjetphi         , "secondjetphi/D");
    tree->Branch("secondjetCHfrac"      , &secondjetCHfrac      , "secondjetCHfrac/D");
    tree->Branch("secondjetNHfrac"      , &secondjetNHfrac      , "secondjetNHfrac/D");
    tree->Branch("secondjetEMfrac"      , &secondjetEMfrac      , "secondjetEMfrac/D");
    tree->Branch("secondjetCEMfrac"     , &secondjetCEMfrac     , "secondjetCEMfrac/D");
    tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");
    tree->Branch("genjetjetdphi"        , &genjetjetdphi        , "genjetjetdphi/D");
    // Photon info
    tree->Branch("phpt"                 , &phpt                 , "phpt/D");
    tree->Branch("pheta"                , &pheta                , "pheta/D");
    tree->Branch("phphi"                , &phphi                , "phphi/D");
    // W/Z gen-level info
    tree->Branch("vid"                  , &vid                  , "vid/I");
    tree->Branch("vmass"                , &vmass                , "vmass/D");
    tree->Branch("vmt"                  , &vmt                  , "vmt/D");
    tree->Branch("vpt"                  , &vpt                  , "vpt/D");
    tree->Branch("veta"                 , &veta                 , "veta/D");
    tree->Branch("vphi"                 , &vphi                 , "vphi/D");
    tree->Branch("l1id"                 , &l1id                 , "l1id/I");
    tree->Branch("l1pt"                 , &l1pt                 , "l1pt/D");
    tree->Branch("l1eta"                , &l1eta                , "l1eta/D");
    tree->Branch("l1phi"                , &l1phi                , "l1phi/D");
    tree->Branch("l2id"                 , &l2id                 , "l2id/I");
    tree->Branch("l2pt"                 , &l2pt                 , "l2pt/D");
    tree->Branch("l2eta"                , &l2eta                , "l2eta/D");
    tree->Branch("l2phi"                , &l2phi                , "l2phi/D");
}

void ZGammaStudies::endJob() {
}

void ZGammaStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_NoiseCleaned");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned");
    triggerPathsVector.push_back("HLT_Photon165");
    triggerPathsVector.push_back("HLT_Photon175_HE10");

    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}

void ZGammaStudies::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ZGammaStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void ZGammaStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void ZGammaStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ZGammaStudies);

