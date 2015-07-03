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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <TH1F.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

class TreeDumper : public edm::EDAnalyzer {
    public:
        explicit TreeDumper(const edm::ParameterSet&);
        ~TreeDumper();
        
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
        edm::InputTag triggerResultsTag;
        edm::InputTag verticesTag;
        edm::InputTag muonsTag;
        edm::InputTag electronsTag;
        edm::InputTag photonsTag;
        edm::InputTag electronVetoIdMap;
        edm::InputTag electronLooseIdMap;
        edm::InputTag electronMediumIdMap;
        edm::InputTag electronTightIdMap;
        edm::InputTag photonLooseIdMap;
        edm::InputTag photonMediumIdMap;
        edm::InputTag photonTightIdMap;
        edm::InputTag photonSIEIEMap;
        edm::InputTag tausTag;
        edm::InputTag jetsTag;
        edm::InputTag fatjetsTag;
        edm::InputTag t1pfmetTag;
        edm::InputTag genevtInfoTag;
        edm::InputTag gensTag;

        // Tokens
        edm::EDGetTokenT<edm::TriggerResults>            triggerResultsToken;
        edm::EDGetTokenT<std::vector<reco::Vertex> >     verticesToken;
        edm::EDGetTokenT<edm::View<pat::Muon> >          muonsToken;
        edm::EDGetTokenT<edm::View<pat::Electron> >      electronsToken;
        edm::EDGetTokenT<edm::View<pat::Photon> >        photonsToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           electronVetoIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           electronLooseIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           electronMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           electronTightIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           photonLooseIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           photonMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >           photonTightIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<float> >          photonSIEIEMapToken;
        edm::EDGetTokenT<edm::View<pat::Tau> >           tausToken;
        edm::EDGetTokenT<edm::View<pat::Jet> >           jetsToken;
        edm::EDGetTokenT<edm::View<pat::Jet> >           fatjetsToken;
        edm::EDGetTokenT<edm::View<pat::MET> >           t1pfmetToken;
        edm::EDGetTokenT<GenEventInfoProduct>            genevtInfoToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle> >  gensToken;

        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        bool isVMCSample;   
        bool uselheweights;   
        bool applyHighMETFilter;
        TTree* tree;

        uint32_t event, run, lumi;
        uint32_t hltmet90, hltmet120, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel;
        uint32_t nvtx;

        double   pfmet, pfmetphi, t1pfmet, t1pfmetphi, mumet, mumetphi, t1mumet, t1mumetphi;
        uint32_t njets, nfatjets, ncntjets, ncntjets2p5, ncntjets3p0;
        double   jetpt[100], jeteta[100], jetphi[100], jetbtag[100], jetCHfrac[100], jetNHfrac[100], jetEMfrac[100], jetCEMfrac[100];
        int32_t  jetid[100], pujetid[100];
        double   fatjetpt[100], fatjeteta[100], fatjetphi[100], fatjetprmass[100], fatjetsdmass[100], fatjettrmass[100], fatjetftmass[100], fatjettau2[100], fatjettau1[100]; 
        double   fatjetCHfrac[100], fatjetNHfrac[100], fatjetEMfrac[100], fatjetCEMfrac[100], fatjetmetdphi[100];
        int32_t  fatjetid[100], pufatjetid[100];

        uint32_t nmuons, nelectrons, ntaus, nphotons, nvetomuons, nvetoelectrons, nvetotaus, nvetophotons;
        double   mupt[100], mueta[100], muphi[100], mupfpt[100], mupfeta[100], mupfphi[100], muiso[100];
        int32_t  muid[100], muidonly[100], mupid[100];
        double   elpt[100], eleta[100], elphi[100];
        int32_t  elidloose[100], elidmedium[100], elidtight[100], elpid[100];
        double   phpt[100], pheta[100], phphi[100], phsieie[100];
        int32_t  phidloose[100], phidmedium[100], phidtight[100];
        double   tapt[100], taeta[100], taphi[100], taiso[100];

        double   vmass, vmt, vpt, veta, vphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
        int32_t  vid, l1id, l2id; 
        double   xsec, wgt, kfact;

        struct PatJetPtSorter {
            bool operator() (const pat::Jet& i, const pat::Jet& j) {
                return (i.pt() > j.pt());
            }
        } jetsorter;
        
        struct PatMuonPtSorter {
            bool operator() (const pat::Muon& i, const pat::Muon& j) {
                return (i.pt() > j.pt());
            }
        } muonsorter;
        
        struct PatElectronPtSorter {
            bool operator() (const pat::Electron& i, const pat::Electron& j) {
                return (i.pt() > j.pt());
            }
        } electronsorter;

        struct PatPhotonPtSorter {
            bool operator() (const pat::Photon& i, const pat::Photon& j) {
                return (i.pt() > j.pt());
            }
        } photonsorter;

        struct PatTauPtSorter {
            bool operator() (const pat::Tau& i, pat::Tau& j) {
                return (i.pt() > j.pt());
            }
        } tausorter;

};

TreeDumper::TreeDumper(const edm::ParameterSet& iConfig): 
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
    muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
    electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
    photonsTag(iConfig.getParameter<edm::InputTag>("photons")),
    electronVetoIdMap(iConfig.getParameter<edm::InputTag>("electronidveto")),
    electronLooseIdMap(iConfig.getParameter<edm::InputTag>("electronidloose")),
    electronMediumIdMap(iConfig.getParameter<edm::InputTag>("electronidmedium")),
    electronTightIdMap(iConfig.getParameter<edm::InputTag>("electronidtight")),
    photonLooseIdMap(iConfig.getParameter<edm::InputTag>("photonidloose")),
    photonMediumIdMap(iConfig.getParameter<edm::InputTag>("photonidmedium")),
    photonTightIdMap(iConfig.getParameter<edm::InputTag>("photonidtight")),
    photonSIEIEMap(iConfig.getParameter<edm::InputTag>("photonsieie")),
    tausTag(iConfig.getParameter<edm::InputTag>("taus")),
    jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
    fatjetsTag(iConfig.getParameter<edm::InputTag>("fatjets")),
    t1pfmetTag(iConfig.getParameter<edm::InputTag>("t1pfmet")),
    genevtInfoTag((iConfig.existsAs<edm::InputTag>("genevt") ? iConfig.getParameter<edm::InputTag>("genevt") : edm::InputTag("generator"))),
    gensTag((iConfig.existsAs<edm::InputTag>("gens") ? iConfig.getParameter<edm::InputTag>("gens") : edm::InputTag("prunedGenParticles"))),
    isVMCSample(iConfig.existsAs<bool>("isVMCSample") ? iConfig.getParameter<bool>("isVMCSample") : false),
    uselheweights(iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
    applyHighMETFilter(iConfig.existsAs<bool>("applyHighMETFilter") ? iConfig.getParameter<bool>("applyHighMETFilter") : false),
    xsec((iConfig.existsAs<double>("xsec") ? iConfig.getParameter<double>("xsec") : 1.0)),
    kfact((iConfig.existsAs<double>("kfactor") ? iConfig.getParameter<double>("kfactor") : 1.0))
{
    // Token consumes instructions
    triggerResultsToken      = consumes<edm::TriggerResults>           (triggerResultsTag); 
    verticesToken            = consumes<std::vector<reco::Vertex> >    (verticesTag);
    muonsToken               = consumes<edm::View<pat::Muon> >         (muonsTag); 
    electronsToken           = consumes<edm::View<pat::Electron> >     (electronsTag); 
    photonsToken             = consumes<edm::View<pat::Photon> >       (photonsTag); 
    electronVetoIdMapToken   = consumes<edm::ValueMap<bool> >          (electronVetoIdMap);
    electronLooseIdMapToken  = consumes<edm::ValueMap<bool> >          (electronLooseIdMap);
    electronMediumIdMapToken = consumes<edm::ValueMap<bool> >          (electronMediumIdMap);
    electronTightIdMapToken  = consumes<edm::ValueMap<bool> >          (electronTightIdMap);
    photonLooseIdMapToken    = consumes<edm::ValueMap<bool> >          (photonLooseIdMap);
    photonMediumIdMapToken   = consumes<edm::ValueMap<bool> >          (photonMediumIdMap);
    photonTightIdMapToken    = consumes<edm::ValueMap<bool> >          (photonTightIdMap);
    photonSIEIEMapToken      = consumes<edm::ValueMap<float> >         (photonSIEIEMap);
    tausToken                = consumes<edm::View<pat::Tau> >          (tausTag); 
    jetsToken                = consumes<edm::View<pat::Jet> >          (jetsTag); 
    fatjetsToken             = consumes<edm::View<pat::Jet> >          (fatjetsTag); 
    t1pfmetToken             = consumes<edm::View<pat::MET> >          (t1pfmetTag); 
    genevtInfoToken          = consumes<GenEventInfoProduct>           (genevtInfoTag);
    gensToken                = consumes<edm::View<reco::GenParticle> > (gensTag); 
  
    // Scaling the cross-section to fb 
    xsec *= 1000.; 

}


TreeDumper::~TreeDumper() {
}

void TreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get handles to all the requisite collections
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByToken(triggerResultsToken, triggerResultsH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);

    Handle<View<pat::Muon> > muonsH;
    iEvent.getByToken(muonsToken, muonsH);

    Handle<View<pat::Electron> > electronsH;
    iEvent.getByToken(electronsToken, electronsH);

    Handle<View<pat::Photon> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    Handle<ValueMap<bool> > electronVetoIdH;
    iEvent.getByToken(electronVetoIdMapToken, electronVetoIdH);

    Handle<ValueMap<bool> > electronLooseIdH;
    iEvent.getByToken(electronLooseIdMapToken, electronLooseIdH);

    Handle<ValueMap<bool> > electronMediumIdH;
    iEvent.getByToken(electronMediumIdMapToken, electronMediumIdH);

    Handle<ValueMap<bool> > electronTightIdH;
    iEvent.getByToken(electronTightIdMapToken, electronTightIdH);

    Handle<ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdMapToken, photonLooseIdH);

    Handle<ValueMap<bool> > photonMediumIdH;
    iEvent.getByToken(photonMediumIdMapToken, photonMediumIdH);

    Handle<ValueMap<bool> > photonTightIdH;
    iEvent.getByToken(photonTightIdMapToken, photonTightIdH);

    Handle<ValueMap<float> > photonSIEIEH;
    iEvent.getByToken(photonSIEIEMapToken, photonSIEIEH);

    Handle<View<pat::Tau> > tausH;
    iEvent.getByToken(tausToken, tausH);

    Handle<View<pat::Jet> > jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<View<pat::Jet> > fatjetsH;
    iEvent.getByToken(fatjetsToken, fatjetsH);

    Handle<View<pat::MET> > t1pfmetH;
    iEvent.getByToken(t1pfmetToken, t1pfmetH);

    Handle<GenEventInfoProduct> genevtInfoH;
    if (uselheweights) iEvent.getByToken(genevtInfoToken, genevtInfoH);

    Handle<View<GenParticle> > gensH;
    if (isVMCSample) iEvent.getByToken(gensToken, gensH);

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
    hltdoublemu  = 0;
    hltsinglemu  = 0;
    hltdoubleel  = 0;
    hltsingleel  = 0;

    // Which triggers fired
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
        if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90     = 1; // MET trigger
        if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120    = 1; // MET trigger
        if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90  = 1; // Jet-MET trigger
        if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120 = 1; // Jet-MET trigger
        if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165 = 1; // Photon trigger
        if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175 = 1; // Photon trigger
        if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu  = 1; // Double muon trigger
        if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu  = 1; // Double muon trigger
        if (i == 8  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 9  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 10 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 11 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 12 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 13 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 14 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 15 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 16 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu  = 1; // Single muon trigger
        if (i == 17 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 18 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 19 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 20 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 21 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 22 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel  = 1; // Double muon trigger
        if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel  = 1; // Single electron trigger
        if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel  = 1; // Single electron trigger
        if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel  = 1; // Single electron trigger
        if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel  = 1; // Single electron trigger
        if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel  = 1; // Single electron trigger
    }

    // Pileup info -- For now just the number of vertices
    nvtx = verticesH->size();

    // Event weight -- Pertinent in the case of aMC@NLO samples where event weights can be negative
    if (uselheweights) wgt = genevtInfoH->weight();
    else wgt = 1.0;

    // MET information 
    pfmet           = t1pfmetH->front().uncorrectedPt();
    pfmetphi        = t1pfmetH->front().uncorrectedPhi();

    t1pfmet         = t1pfmetH->front().et();
    t1pfmetphi      = t1pfmetH->front().phi();

    mumet           = pfmet;
    mumetphi        = pfmetphi;

    t1mumet         = t1pfmet;
    t1mumetphi      = t1pfmetphi;

    double mumetx   = mumet * cos(mumetphi);
    double mumety   = mumet * sin(mumetphi);

    double t1mumetx = t1mumet * cos(t1mumetphi);
    double t1mumety = t1mumet * sin(t1mumetphi);

    // Muon information
    pat::MuonCollection muons;
    for (View<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        if (muons_iter->pt() > 10. && fabs(muons_iter->eta()) < 2.4 && muon::isLooseMuon(*muons_iter)) {
            pat::Muon muon = *muons_iter;
            muons.push_back(muon);
        }
    }
    sort(muons.begin(), muons.end(), muonsorter);

    nmuons     = muons.size();
    if (nmuons > 100) nmuons = 100;
    nvetomuons = 0;

    for (size_t i = 0; i < nmuons; i++) {
        mupt[i]  = muons[i].pt();
        mueta[i] = muons[i].eta();
        muphi[i] = muons[i].phi();
        mupid[i] = muons[i].pdgId();

        mupfpt[i]  = muons[i].pfP4().Pt();
        mupfeta[i] = muons[i].pfP4().Eta();
        mupfphi[i] = muons[i].pfP4().Phi();

        float isoval = muons[i].pfIsolationR04().sumNeutralHadronEt;
        isoval += muons[i].pfIsolationR04().sumPhotonEt;
        isoval -= 0.5*muons[i].pfIsolationR04().sumPUPt;
        if (isoval < 0.) isoval = 0.;
        isoval += muons[i].pfIsolationR04().sumChargedHadronPt;
        isoval /= muons[i].pt();
        muiso[i] = isoval;

        if (isoval < 0.2) {
            nvetomuons++;
            mumetx   += mupfpt[i] * cos(mupfphi[i]);
            mumety   += mupfpt[i] * sin(mupfphi[i]);
            t1mumetx += mupfpt[i] * cos(mupfphi[i]);
            t1mumety += mupfpt[i] * sin(mupfphi[i]);
        }

        muidonly[i] = 0;
        muid[i] = 0;
        if (nvtx > 0 && muon::isTightMuon(muons[i], *(verticesH->begin()))) {
            muidonly[i] = 1;
            if (isoval < 0.12) muid[i] = 1;
        }
    }
    
    mumet      = sqrt(mumetx*mumetx + mumety*mumety);
    mumetphi   = atan2(mumety, mumetx);

    t1mumet    = sqrt(t1mumetx*t1mumetx + t1mumety*t1mumety);
    t1mumetphi = atan2(t1mumety, t1mumetx);

    if (applyHighMETFilter && mumet < 200. && t1mumet < 200.) return;

    // Electron information
    pat::ElectronCollection electrons;
    vector<int> eidloose;
    vector<int> eidmedium;
    vector<int> eidtight;
    for (View<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        const Ptr<pat::Electron> electronPtr(electronsH, electrons_iter - electronsH->begin());
        if (electrons_iter->pt() > 10. && fabs(electrons_iter->superCluster()->eta()) < 2.5 && (*electronVetoIdH)[electronPtr]) {
            pat::Electron electron = *electrons_iter;
            electrons.push_back(electron);
            if ((*electronLooseIdH) [electronPtr]) eidloose .push_back(1);
            else                                   eidloose .push_back(0);
            if ((*electronMediumIdH)[electronPtr]) eidmedium.push_back(1);
            else                                   eidmedium.push_back(0);
            if ((*electronTightIdH) [electronPtr]) eidtight .push_back(1);
            else                                   eidtight .push_back(0);
        }
    }
    sort(electrons.begin(), electrons.end(), electronsorter);

    nelectrons     = electrons.size();
    nvetoelectrons = nelectrons;
    if (nelectrons > 100) nelectrons = 100;

    for (size_t i = 0; i < nelectrons; i++) {
        elpt[i]  = electrons[i].pt();
        eleta[i] = electrons[i].eta();
        elphi[i] = electrons[i].phi();
        elpid[i] = electrons[i].pdgId();

        elidloose [i] = eidloose[i];
        elidmedium[i] = eidmedium[i];
        elidtight [i] = eidtight[i];
    }

    // Photon information
    pat::PhotonCollection photons;
    vector<int> pidloose;
    vector<int> pidmedium;
    vector<int> pidtight;
    vector<double> psieie;
    for (View<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
        const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
        if (photons_iter->pt() > 15. && fabs(photons_iter->superCluster()->eta()) < 2.5) {
            if (photons_iter->pt() > 0.8 || photons_iter->chargedHadronIso() < 20. || photons_iter->chargedHadronIso() < 0.3*photons_iter->pt()) {
                pat::Photon photon = *photons_iter;
                photons.push_back(photon);
                if ((*photonLooseIdH) [photonPtr]) pidloose .push_back(1);
                else                               pidloose .push_back(0);
                if ((*photonMediumIdH)[photonPtr]) pidmedium.push_back(1);
                else                               pidmedium.push_back(0);
                if ((*photonTightIdH) [photonPtr]) pidtight .push_back(1);
                else                               pidtight .push_back(0);
                psieie.push_back((*photonSIEIEH)[photonPtr]);
            }
        }
    }
    sort(photons.begin(), photons.end(), photonsorter);

    nphotons     = photons.size();
    if (nphotons > 100) nphotons = 100;
    nvetophotons = 0;

    for (size_t i = 0; i < nphotons; i++) {
        phpt[i]  = photons[i].pt();
        pheta[i] = photons[i].eta();
        phphi[i] = photons[i].phi();

        phidloose [i] = pidloose[i];
        phidmedium[i] = pidmedium[i];
        phidtight [i] = pidtight[i];

        phsieie[i] = psieie[i];

        if (pidloose[i] == 1) nvetophotons++;
    }

    // Tau information
    pat::TauCollection taus;
    for (View<pat::Tau>::const_iterator taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
        const Ptr<pat::Tau> tauPtr(tausH, taus_iter - tausH->begin());
        if (taus_iter->pt() > 18. && fabs(taus_iter->eta()) < 2.3 && taus_iter->tauID("decayModeFinding")) {
            pat::Tau tau = *taus_iter;
            taus.push_back(tau);
        }
    }
    sort(taus.begin(), taus.end(), tausorter);

    ntaus     = taus.size();
    if (ntaus > 100) ntaus = 100.;
    nvetotaus = 0;

    for (size_t i = 0; i < ntaus; i++) {
        tapt[i]  = taus[i].pt();
        taeta[i] = taus[i].eta();
        taphi[i] = taus[i].phi();
        taiso[i] = taus[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        if (taiso[i] < 5.) nvetotaus++;
    }

    // AK4 Jets information
    pat::JetCollection jets;
    for (View<pat::Jet>::const_iterator jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
        if (fabs(jets_iter->pt()) < 10.) continue;
        pat::Jet jet = *jets_iter;
        jets.push_back(jet);
    }
    sort(jets.begin(), jets.end(), jetsorter);

    njets       = jets.size();
    if (njets   > 100) njets = 100;
    ncntjets    = 0;
    ncntjets2p5 = 0;
    ncntjets3p0 = 0;

    for (size_t i = 0; i < njets; i++) {
        jetpt[i]  = jets[i].pt();
        jeteta[i] = jets[i].eta();
        jetphi[i] = jets[i].phi();

        jetCHfrac[i]  = jets[i].chargedHadronEnergyFraction();
        jetNHfrac[i]  = jets[i].neutralHadronEnergyFraction();
        jetEMfrac[i]  = jets[i].neutralEmEnergyFraction();
        jetCEMfrac[i] = jets[i].chargedEmEnergyFraction();

        jetbtag[i]    = jets[0].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        double jetabseta = fabs(jeteta[i]);
        jetid[i] = 0;
        if (jetNHfrac[i] < 0.99 && jetEMfrac[i] < 0.99 && (jets[i].chargedMultiplicity() + jets[i].neutralMultiplicity()) > 1 && jets[i].muonEnergyFraction() < 0.8) {
            if (jetabseta > 2.4) jetid[i] = 1;
            else if (jetabseta <= 2.4 && jetCHfrac[i] > 0. && jetCEMfrac[i] < 0.99 && jets[i].chargedMultiplicity() > 0) jetid[i] = 1;
        }
        pujetid[i] = 0;
        double puidval = jets[i].userFloat("pileupJetId:fullDiscriminant");
        if (jetabseta >= 0.00 && jetabseta < 2.50 && puidval > -0.63) pujetid[i] = 1;
        if (jetabseta >= 2.50 && jetabseta < 2.75 && puidval > -0.60) pujetid[i] = 1;
        if (jetabseta >= 2.75 && jetabseta < 3.00 && puidval > -0.55) pujetid[i] = 1;
        if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.45) pujetid[i] = 1;

        if (jetpt[i] > 30 && jetid[i] == 1 && pujetid[i] == 1) {
            ncntjets++;
            if (jetabseta < 2.5) ncntjets2p5++;
            if (jetabseta < 3.0) ncntjets3p0++;
        }
    }

    // Fat Jets information
    pat::JetCollection fatjets;
    for (View<pat::Jet>::const_iterator fatjets_iter = fatjetsH->begin(); fatjets_iter != fatjetsH->end(); ++fatjets_iter) {
        if (fabs(fatjets_iter->pt()) < 200.) continue;
        pat::Jet fatjet = *fatjets_iter;
        fatjets.push_back(fatjet);
    }
    sort(fatjets.begin(), fatjets.end(), jetsorter);

    nfatjets = fatjets.size();
    if (nfatjets > 100) nfatjets = 100;
    
    for (size_t i = 0; i < nfatjets; i++) {
        fatjetpt[i]  = fatjets[i].pt();
        fatjeteta[i] = fatjets[i].eta();
        fatjetphi[i] = fatjets[i].phi();

        fatjetCHfrac[i]  = fatjets[i].chargedHadronEnergyFraction();
        fatjetNHfrac[i]  = fatjets[i].neutralHadronEnergyFraction();
        fatjetEMfrac[i]  = fatjets[i].neutralEmEnergyFraction();
        fatjetCEMfrac[i] = fatjets[i].chargedEmEnergyFraction();

        fatjetprmass[i]  = fatjets[i].userFloat("ak8PFJetsCHSPrunedMass");    
        fatjetsdmass[i]  = fatjets[i].userFloat("ak8PFJetsCHSSoftDropMass");    
        fatjettrmass[i]  = fatjets[i].userFloat("ak8PFJetsCHSTrimmedMass");    
        fatjetftmass[i]  = fatjets[i].userFloat("ak8PFJetsCHSFilteredMass");    

        double fatjetabseta = fabs(fatjeteta[i]);
        fatjetid[i] = 0;
        if (fatjetNHfrac[i] < 0.99 && fatjetEMfrac[i] < 0.99 && (fatjets[i].chargedMultiplicity() + fatjets[i].neutralMultiplicity()) > 1 && fatjets[i].muonEnergyFraction() < 0.8) {
            if (fatjetabseta > 2.4) fatjetid[i] = 1;
            else if (fatjetabseta <= 2.4 && fatjetCHfrac[i] > 0. && fatjetCEMfrac[i] < 0.99 && fatjets[i].chargedMultiplicity() > 0) fatjetid[i] = 1;
        }
        pufatjetid[i] = 0;
        double puidval = fatjets[i].userFloat("pileupJetId:fullDiscriminant");
        if (fatjetabseta >= 0.00 && fatjetabseta < 2.50 && puidval > -0.63) pufatjetid[i] = 1;
        if (fatjetabseta >= 2.50 && fatjetabseta < 2.75 && puidval > -0.60) pufatjetid[i] = 1;
        if (fatjetabseta >= 2.75 && fatjetabseta < 3.00 && puidval > -0.55) pufatjetid[i] = 1;
        if (fatjetabseta >= 3.00 && fatjetabseta < 5.00 && puidval > -0.45) pufatjetid[i] = 1;
    }

    // Generator-level information
    vid           = 0;
    vmass         = 0.0;
    vmt           = 0.0;
    vpt           = 0.0;
    veta          = 0.0;
    vphi          = 0.0;
    l1id          = 0;
    l1pt          = 0.0;
    l1eta         = 0.0;
    l1phi         = 0.0;
    l2id          = 0;
    l2pt          = 0.0;
    l2eta         = 0.0;
    l2phi         = 0.0;

    if (isVMCSample && gensH.isValid()) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() > 1 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) {
                vid   = gens_iter->pdgId();
                vmass = gens_iter->mass();
                vpt   = gens_iter->pt();
                veta  = gens_iter->eta();
                vphi  = gens_iter->phi();
                l1id   = gens_iter->daughter(0)->pdgId();
                l1pt   = gens_iter->daughter(0)->pt();
                l1eta  = gens_iter->daughter(0)->eta();
                l1phi  = gens_iter->daughter(0)->phi();
                l2id   = gens_iter->daughter(1)->pdgId();
                l2pt   = gens_iter->daughter(1)->pt();
                l2eta  = gens_iter->daughter(1)->eta();
                l2phi  = gens_iter->daughter(1)->phi();
                vmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi)))); 
            }
            if (isVMCSample && gens_iter->pdgId() == 22 && gens_iter->status() == 1 && gens_iter->pt() > 100.) {
                vid   = gens_iter->pdgId();
                vpt   = gens_iter->pt();
                veta  = gens_iter->eta();
                vphi  = gens_iter->phi();
            }
        }
    }

    if (isVMCSample && gensH.isValid() && vid == 0) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if (gens_iter->pdgId() == 22 && gens_iter->status() == 1 && gens_iter->pt() > 100. && gens_iter->pt() > vpt) {
                vid   = gens_iter->pdgId();
                vpt   = gens_iter->pt();
                veta  = gens_iter->eta();
                vphi  = gens_iter->phi();
            }
        }
    }

    tree->Fill();

}


void TreeDumper::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"       , "tree");
    // Run, Lumi, Event info
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
    // Event weights
    tree->Branch("xsec"                 , &xsec                 , "xsec/D");
    tree->Branch("wgt"                  , &wgt                  , "wgt/D");
    tree->Branch("kfact"                , &kfact                , "kfact/D");
    // Pileup info
    tree->Branch("nvtx"                 , &nvtx                 , "nvtx/i");
    // Triggers and event filters
    tree->Branch("hltmet90"             , &hltmet90             , "hltmet90/i");
    tree->Branch("hltmet120"            , &hltmet120            , "hltmet120/i");
    tree->Branch("hltjetmet90"          , &hltjetmet90          , "hltjetmet90/i");
    tree->Branch("hltjetmet120"         , &hltjetmet120         , "hltjetmet120/i");
    tree->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/i");
    tree->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/i");
    tree->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/i");
    tree->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/i");
    tree->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/i");
    tree->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/i");
    // Object counts
    tree->Branch("nmuons"               , &nmuons               , "nmuons/i");
    tree->Branch("nelectrons"           , &nelectrons           , "nelectrons/i");
    tree->Branch("nphotons"             , &nphotons             , "nphotons/i");
    tree->Branch("ntaus"                , &ntaus                , "ntaus/i");
    tree->Branch("njets"                , &njets                , "njets/i");
    tree->Branch("nfatjets"             , &nfatjets             , "nfatjets/i");
    tree->Branch("nvetomuons"           , &nvetomuons           , "nvetomuons/i");
    tree->Branch("nvetoelectrons"       , &nvetoelectrons       , "nvetoelectrons/i");
    tree->Branch("nvetophotons"         , &nvetophotons         , "nvetophotons/i");
    tree->Branch("nvetotaus"            , &nvetotaus            , "nvetotaus/i");
    // MET info
    tree->Branch("pfmet"                , &pfmet                , "pfmet/D");
    tree->Branch("pfmetphi"             , &pfmetphi             , "pfmetphi/D");
    tree->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/D");
    tree->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/D");
    tree->Branch("mumet"                , &mumet                , "mumet/D");
    tree->Branch("mumetphi"             , &mumetphi             , "mumetphi/D");
    tree->Branch("t1mumet"              , &t1mumet              , "t1mumet/D");
    tree->Branch("t1mumetphi"           , &t1mumetphi           , "t1mumetphi/D");
    // Jet info
    tree->Branch("jetpt"                , jetpt                 , "jetpt[njets]/D");
    tree->Branch("jeteta"               , jeteta                , "jeteta[njets]/D");
    tree->Branch("jetphi"               , jetphi                , "jetphi[njets]/D");
    tree->Branch("jetbtag"              , jetbtag               , "jetbtag[njets]/D");
    tree->Branch("jetCHfrac"            , jetCHfrac             , "jetCHfrac[njets]/D");
    tree->Branch("jetNHfrac"            , jetNHfrac             , "jetNHfrac[njets]/D");
    tree->Branch("jetEMfrac"            , jetEMfrac             , "jetEMfrac[njets]/D");
    tree->Branch("jetCEMfrac"           , jetCEMfrac            , "jetCEMfrac[njets]/D");
    tree->Branch("jetid"                , jetid                 , "jetid[njets]/I");
    tree->Branch("pujetid"              , pujetid               , "pujetid[njets]/I");
    tree->Branch("fatjetpt"             , fatjetpt              , "fatjetpt[nfatjets]/D");
    tree->Branch("fatjeteta"            , fatjeteta             , "fatjeteta[nfatjets]/D");
    tree->Branch("fatjetphi"            , fatjetphi             , "fatjetphi[nfatjets]/D");
    tree->Branch("fatjetprmass"         , fatjetprmass          , "fatjetprmass[nfatjets]/D");
    tree->Branch("fatjetsdmass"         , fatjetsdmass          , "fatjetsdmass[nfatjets]/D");
    tree->Branch("fatjettrmass"         , fatjettrmass          , "fatjettrmass[nfatjets]/D");
    tree->Branch("fatjetftmass"         , fatjetftmass          , "fatjetftmass[nfatjets]/D");
    tree->Branch("fatjettau2"           , fatjettau2            , "fatjettau2[nfatjets]/D");
    tree->Branch("fatjettau1"           , fatjettau1            , "fatjettau1[nfatjets]/D");
    tree->Branch("fatjetCHfrac"         , fatjetCHfrac          , "fatjetCHfrac[nfatjets]/D");
    tree->Branch("fatjetNHfrac"         , fatjetNHfrac          , "fatjetNHfrac[nfatjets]/D");
    tree->Branch("fatjetEMfrac"         , fatjetEMfrac          , "fatjetEMfrac[nfatjets]/D");
    tree->Branch("fatjetCEMfrac"        , fatjetCEMfrac         , "fatjetCEMfrac[nfatjets]/D");
    tree->Branch("fatjetid"             , fatjetid              , "fatjetid[nfatjets]/I");
    tree->Branch("pufatjetid"           , pufatjetid            , "pufatjetid[nfatjets]/I");
    // Lepton-Photon info
    tree->Branch("mupt"                 , mupt                  , "mupt[nmuons]/D");
    tree->Branch("mueta"                , mueta                 , "mueta[nmuons]/D");
    tree->Branch("muphi"                , muphi                 , "muphi[nmuons]/D");
    tree->Branch("muidonly"             , muidonly              , "muidonly[nmuons]/I");
    tree->Branch("muid"                 , muid                  , "muid[nmuons]/I");
    tree->Branch("mupid"                , mupid                 , "mupid[nmuons]/I");
    tree->Branch("muiso"                , muiso                 , "muiso[nmuons]/D");
    tree->Branch("elpt"                 , elpt                  , "elpt[nelectrons]/D");
    tree->Branch("eleta"                , eleta                 , "eleta[nelectrons]/D");
    tree->Branch("elphi"                , elphi                 , "elphi[nelectrons]/D");
    tree->Branch("elidloose"            , elidloose             , "elidloose[nelectrons]/I");
    tree->Branch("elidmedium"           , elidmedium            , "elidmedium[nelectrons]/I");
    tree->Branch("elidtight"            , elidtight             , "elidtight[nelectrons]/I");
    tree->Branch("elpid"                , elpid                 , "elpid[nelectrons]/I");
    tree->Branch("phpt"                 , phpt                  , "phpt[nphotons]/D");
    tree->Branch("pheta"                , pheta                 , "pheta[nphotons]/D");
    tree->Branch("phphi"                , phphi                 , "phphi[nphotons]/D");
    tree->Branch("phidloose"            , phidloose             , "phidloose[nphotons]/I");
    tree->Branch("phidmedium"           , phidmedium            , "phidmedium[nphotons]/I");
    tree->Branch("phidtight"            , phidtight             , "phidtight[nphotons]/I");
    tree->Branch("phsieie"              , phsieie               , "phsieie[nphotons]/D");
    tree->Branch("tapt"                 , tapt                  , "tapt[ntaus]/D");
    tree->Branch("taeta"                , taeta                 , "taeta[ntaus]/D");
    tree->Branch("taphi"                , taphi                 , "taphi[ntaus]/D");
    tree->Branch("taiso"                , taiso                 , "taiso[ntaus]/D");
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

void TreeDumper::endJob() {
}

void TreeDumper::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // Trigger Paths
    triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_NoiseCleaned");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_NoiseCleaned");
    triggerPathsVector.push_back("HLT_Photon165");
    triggerPathsVector.push_back("HLT_Photon175_HE10");
    triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
    triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
    triggerPathsVector.push_back("HLT_IsoMu17_eta2p1");
    triggerPathsVector.push_back("HLT_IsoMu20_eta2p1");
    triggerPathsVector.push_back("HLT_IsoMu24_eta2p1");
    triggerPathsVector.push_back("HLT_IsoMu20");
    triggerPathsVector.push_back("HLT_IsoMu27");
    triggerPathsVector.push_back("HLT_IsoTkMu20_eta2p1");
    triggerPathsVector.push_back("HLT_IsoTkMu24_eta2p1");
    triggerPathsVector.push_back("HLT_IsoTkMu20");
    triggerPathsVector.push_back("HLT_IsoTkMu27");
    triggerPathsVector.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf");
    triggerPathsVector.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW");
    triggerPathsVector.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL");
    triggerPathsVector.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300");
    triggerPathsVector.push_back("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL");
    triggerPathsVector.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    triggerPathsVector.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf");
    triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPTight_Gsf");
    triggerPathsVector.push_back("HLT_Ele32_eta2p1_WPLoose_Gsf");
    triggerPathsVector.push_back("HLT_Ele32_eta2p1_WPTight_Gsf");
    triggerPathsVector.push_back("HLT_Ele27_WPLoose_Gsf_WHbbBoost");

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

void TreeDumper::endRun(edm::Run const&, edm::EventSetup const&) {
}

void TreeDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void TreeDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void TreeDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TreeDumper);
