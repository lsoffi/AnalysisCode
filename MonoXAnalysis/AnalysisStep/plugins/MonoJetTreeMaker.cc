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
#include "FWCore/Common/interface/TriggerNames.h"
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
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
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

class MonoJetTreeMaker : public edm::EDAnalyzer {
    public:
        explicit MonoJetTreeMaker(const edm::ParameterSet&);
        ~MonoJetTreeMaker();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        void findMother(const reco::Candidate*, int &, double &, double &, double &);
        void findFirstNonPhotonMother(const reco::Candidate*, int &, double &, double &, double &);
        std::string concatenate(std::vector<std::string> vstring);
        void Init();

        // InputTags
        edm::InputTag triggerResultsTag;
        edm::InputTag _IT_trg_obj;
        edm::InputTag filterResultsTag;

        // Tokens
        edm::EDGetTokenT<edm::TriggerResults>              triggerResultsToken;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trgObjectsToken_;
        edm::EDGetTokenT<edm::TriggerResults>              filterResultsToken;
        edm::EDGetTokenT<HcalNoiseSummary>                 hcalnoiseToken;
        edm::EDGetTokenT<bool>                             hbhelooseToken;
        edm::EDGetTokenT<bool>                             hbhetightToken;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
        edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
        edm::EDGetTokenT<std::vector<reco::Vertex> >       verticesToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
        edm::EDGetTokenT<pat::MuonRefVector>               muonsToken;
        edm::EDGetTokenT<pat::ElectronRefVector>           electronsToken;
        edm::EDGetTokenT<pat::PhotonRefVector>             photonsToken;
        edm::EDGetTokenT<pat::MuonRefVector>               tightmuonsToken;
        edm::EDGetTokenT<pat::ElectronRefVector>           tightelectronsToken;
        edm::EDGetTokenT<pat::PhotonRefVector>             tightphotonsToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonMediumIdToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonTightIdToken;
        edm::EDGetTokenT<std::vector<pat::Tau> >           tausToken;
        edm::EDGetTokenT<std::vector<pat::Jet> >           jetsToken;
        edm::EDGetTokenT<std::vector<pat::Jet> >           fatjetsToken;
        edm::EDGetTokenT<edm::ValueMap<float> >            qglToken;
        edm::EDGetTokenT<edm::ValueMap<float> >            qgs2Token;
        edm::EDGetTokenT<edm::ValueMap<int> >              qgmultToken;
        edm::EDGetTokenT<edm::ValueMap<float> >            qgptdToken;
        edm::EDGetTokenT<edm::View<pat::MET> >             t1pfmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            partmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            pfmuptToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            mumetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            t1mumetToken;

        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        std::vector<std::string> filterPathsVector;
        std::map<std::string, int> filterPathsMap;
        bool applyHLTFilter;
        bool isWorZMCSample, isSignalSample;   
        bool cleanMuonJet, cleanElectronJet, cleanPhotonJet;   
        bool uselheweights;   
        TTree* tree;

        int32_t  puobs, putrue; 
        int32_t  wzid, l1id, l2id, mu1pid, mu2pid, mu1id, mu2id, el1pid, el2pid, el1id, el2id, phidm, phidt, parid, ancid; 
        uint32_t event, run, lumi;
        uint32_t nvtx, nmuons, nelectrons, ntaus, ntightmuons, ntightelectrons, nphotons, njets, nbjets, nfatjets;
        uint8_t  hltmet90, hltmet120, hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel;
        uint8_t  flagcsctight, flaghbhenoise, flaghbheloose, flaghbhetight, flaghcallaser, flagecaltrig, flageebadsc, flagecallaser, flagtrkfail, flagtrkpog, flaghnoiseloose, flaghnoisetight, flaghnoisehilvl;
        double   pfmet, pfmetphi, t1pfmet, t1pfmetphi, calomet, calometphi, pfmupt, pfmuphi, mumet, mumetphi, t1mumet, t1mumetphi;
        double   hmet, hmetphi, amet, ametphi, bmet, bmetphi, cmet, cmetphi, emet, emetphi, mmet, mmetphi, pmet, pmetphi, omet, ometphi;
        double   fatjetpt, fatjeteta, fatjetphi, fatjettau2, fatjettau1, fatjetCHfrac, fatjetNHfrac, fatjetEMfrac, fatjetCEMfrac, fatjetmetdphi, fatjetprmass, fatjetsdmass, fatjettrmass, fatjetftmass;
        double   signaljetpt, signaljeteta, signaljetphi, signaljetbtag, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi, signaljetqgl, signaljetqgs2, signaljetqgptd;
        double   secondjetpt, secondjeteta, secondjetphi, secondjetbtag, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetCEMfrac, secondjetmetdphi, secondjetqgl, secondjetqgs2, secondjetqgptd;
        double   thirdjetpt , thirdjeteta , thirdjetphi , thirdjetbtag , thirdjetCHfrac , thirdjetNHfrac , thirdjetEMfrac , thirdjetCEMfrac , thirdjetmetdphi , thirdjetqgl , thirdjetqgs2 , thirdjetqgptd ;
        int      signaljetqgmult, secondjetqgmult, thirdjetqgmult;
        double   jetjetdphi, jetmetdphimin, incjetmetdphimin;
        double   ht, dht, mht, alphat, apcjetmetmax, apcjetmetmin; 
        double   wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi, parpt, pareta, parphi, ancpt, anceta, ancphi;
        double   mu1pt, mu1eta, mu1phi, mu1pfpt, mu1pfeta, mu1pfphi, mu2pt, mu2eta, mu2phi, mu2pfpt, mu2pfeta, mu2pfphi, el1pt, el1eta, el1phi, el2pt, el2eta, el2phi, phpt, pheta, phphi;
        double   zmass, zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt;
        double   xsec, wgt, kfact, puwgt;
        int32_t  _verbose;

        // Trigger objects
        TString _trig_pass;
        Int_t _trig_n;
        //
        Int_t _trig_obj_n;
        std::vector< double > _trig_obj_pt, _trig_obj_eta, _trig_obj_phi;
        std::vector< std::string > _trig_obj_col, _trig_obj_lab;
        //std::vector< std::string > _trig_obj_path_FF, _trig_obj_path_FT, 
        //_trig_obj_path_TF, _trig_obj_path_TT ;
        //std::vector< std::vector<int> > _trig_obj_ids;

        // Jet informations
        const static UInt_t _nMaxJets=20;
        uint32_t _jet_n;
        double _jet_pt[_nMaxJets], _jet_eta[_nMaxJets], _jet_phi[_nMaxJets], 
	  _jet_CHfrac[_nMaxJets], _jet_NHfrac[_nMaxJets], _jet_EMfrac[_nMaxJets], _jet_CEMfrac[_nMaxJets];

        struct PatJetPtSorter {
            bool operator() (pat::JetRef i, pat::JetRef j) {
                return (i->pt() > j->pt());
            }
        } jetsorter;
        
        struct PatMuonPtSorter {
            bool operator() (pat::MuonRef i, pat::MuonRef j) {
                return (i->pt() > j->pt());
            }
        } muonsorter;
        
        struct PatElectronPtSorter {
            bool operator() (pat::ElectronRef i, pat::ElectronRef j) {
                return (i->pt() > j->pt());
            }
        } electronsorter;

        struct PatPhotonPtSorter {
            bool operator() (pat::PhotonRef i, pat::PhotonRef j) {
                return (i->pt() > j->pt());
            }
        } photonsorter;

};

MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    _IT_trg_obj(iConfig.getParameter<edm::InputTag>("objects")),
    filterResultsTag(iConfig.getParameter<edm::InputTag>("filterResults")),
    triggerResultsToken(consumes<edm::TriggerResults> (triggerResultsTag)),
    trgObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(_IT_trg_obj)),
    filterResultsToken(consumes<edm::TriggerResults> (filterResultsTag)),
    hcalnoiseToken(consumes<HcalNoiseSummary> (iConfig.getParameter<edm::InputTag>("hcalnoise"))),
    hbhelooseToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbheloose"))),
    hbhetightToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbhetight"))),
    pileupInfoToken(consumes<std::vector<PileupSummaryInfo> > ((iConfig.existsAs<edm::InputTag>("pileup") ? iConfig.getParameter<edm::InputTag>("pileup") : edm::InputTag("addPileupInfo")))),
    genevtInfoToken(consumes<GenEventInfoProduct> ((iConfig.existsAs<edm::InputTag>("genevt") ? iConfig.getParameter<edm::InputTag>("genevt") : edm::InputTag("generator")))),
    verticesToken(consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
    gensToken(consumes<edm::View<reco::GenParticle> > ((iConfig.existsAs<edm::InputTag>("gens") ? iConfig.getParameter<edm::InputTag>("gens") : edm::InputTag("prunedGenParticles")))),
    muonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("electrons"))),
    photonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photons"))),
    tightmuonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("tightmuons"))),
    tightelectronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("tightelectrons"))),
    tightphotonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("tightphotons"))),
    photonMediumIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonMediumId"))),
    photonTightIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonTightId"))),
    tausToken(consumes<std::vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("taus"))),
    jetsToken(consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetsToken(consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("fatjets"))),
    qglToken(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("qgl"))),
    qgs2Token(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("qgs2"))),
    qgmultToken(consumes<edm::ValueMap<int> > (iConfig.getParameter<edm::InputTag>("qgmult"))),
    qgptdToken(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("qgptd"))),
    t1pfmetToken(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("t1pfmet"))),
    partmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("partmet"))),
    pfmuptToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("pfmupt"))),
    mumetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("mumet"))),
    t1mumetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("t1mumet"))),
    applyHLTFilter(iConfig.existsAs<bool>("applyHLTFilter") ? iConfig.getParameter<bool>("applyHLTFilter") : false),
    isWorZMCSample(iConfig.existsAs<bool>("isWorZMCSample") ? iConfig.getParameter<bool>("isWorZMCSample") : false),
    isSignalSample(iConfig.existsAs<bool>("isSignalSample") ? iConfig.getParameter<bool>("isSignalSample") : false),
    cleanMuonJet(iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
    cleanElectronJet(iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
    cleanPhotonJet(iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
    uselheweights(iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
    xsec(iConfig.getParameter<double>("xsec") * 1000.0),
    kfact((iConfig.existsAs<double>("kfactor") ? iConfig.getParameter<double>("kfactor") : 1.0)),
    _verbose((iConfig.existsAs<int32_t>("verbose") ? iConfig.getParameter<int32_t>("verbose") : 0))
{
}


MonoJetTreeMaker::~MonoJetTreeMaker() {
}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    if(_verbose>3) std::cout << "- analyze(...) starts" << std::endl;
    Init();

    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get handles to all the requisite collections
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByToken(triggerResultsToken, triggerResultsH);

    edm::Handle<pat::TriggerObjectStandAloneCollection> H_trg_obj;
    iEvent.getByToken(trgObjectsToken_, H_trg_obj);

    Handle<TriggerResults> filterResultsH;
    iEvent.getByToken(filterResultsToken, filterResultsH);

    Handle<HcalNoiseSummary> hcalnoiseH;
    iEvent.getByToken(hcalnoiseToken, hcalnoiseH);

    Handle<bool> hbhelooseH;
    iEvent.getByToken(hbhelooseToken, hbhelooseH);

    Handle<bool> hbhetightH;
    iEvent.getByToken(hbhetightToken, hbhetightH);

    Handle<vector<PileupSummaryInfo> > pileupInfoH;
    iEvent.getByToken(pileupInfoToken, pileupInfoH);

    Handle<GenEventInfoProduct> genevtInfoH;
    if (uselheweights) iEvent.getByToken(genevtInfoToken, genevtInfoH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);

    Handle<View<GenParticle> > gensH;
    if (isWorZMCSample || isSignalSample) iEvent.getByToken(gensToken, gensH);

    Handle<pat::MuonRefVector> muonsH;
    iEvent.getByToken(muonsToken, muonsH);
    pat::MuonRefVector muons = *muonsH;

    Handle<pat::ElectronRefVector> electronsH;
    iEvent.getByToken(electronsToken, electronsH);
    pat::ElectronRefVector electrons = *electronsH;

    Handle<pat::PhotonRefVector> photonsH;
    iEvent.getByToken(photonsToken, photonsH);
    pat::PhotonRefVector photons = *photonsH;

    Handle<pat::MuonRefVector> tightmuonsH;
    iEvent.getByToken(tightmuonsToken, tightmuonsH);
    pat::MuonRefVector tightmuons = *tightmuonsH;

    Handle<pat::ElectronRefVector> tightelectronsH;
    iEvent.getByToken(tightelectronsToken, tightelectronsH);
    pat::ElectronRefVector tightelectrons = *tightelectronsH;

    Handle<pat::PhotonRefVector> tightphotonsH;
    iEvent.getByToken(tightphotonsToken, tightphotonsH);
    pat::PhotonRefVector tightphotons = *tightphotonsH;

    Handle<edm::ValueMap<bool> > photonMediumIdH;
    iEvent.getByToken(photonMediumIdToken, photonMediumIdH);

    Handle<edm::ValueMap<bool> > photonTightIdH;
    iEvent.getByToken(photonTightIdToken, photonTightIdH);

    Handle<std::vector<pat::Tau> > tausH;
    iEvent.getByToken(tausToken, tausH);

    Handle<std::vector<pat::Jet> > jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<std::vector<pat::Jet> > fatjetsH;
    iEvent.getByToken(fatjetsToken, fatjetsH);

    Handle<edm::ValueMap<float> > qglH;
    iEvent.getByToken(qglToken, qglH);

    Handle<edm::ValueMap<float> > qgs2H;
    iEvent.getByToken(qgs2Token, qgs2H);

    Handle<edm::ValueMap<int> > qgmultH;
    iEvent.getByToken(qgmultToken, qgmultH);

    Handle<edm::ValueMap<float> > qgptdH;
    iEvent.getByToken(qgptdToken, qgptdH);

    Handle<View<pat::MET> > t1pfmetH;
    iEvent.getByToken(t1pfmetToken, t1pfmetH);

    Handle<View<reco::MET> > partmetH;
    iEvent.getByToken(partmetToken, partmetH);

    Handle<View<reco::MET> > pfmuptH;
    iEvent.getByToken(pfmuptToken, pfmuptH);

    Handle<View<reco::MET> > mumetH;
    iEvent.getByToken(mumetToken, mumetH);

    Handle<View<reco::MET> > t1mumetH;
    iEvent.getByToken(t1mumetToken, t1mumetH);

    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();

    // Trigger info
    hltmet90        = 0;
    hltmet120       = 0;
    hltmetwithmu90  = 0;
    hltmetwithmu120 = 0;
    hltmetwithmu170 = 0;
    hltmetwithmu300 = 0;
    hltjetmet90     = 0;
    hltjetmet120    = 0;
    hltphoton165    = 0;
    hltphoton175    = 0;
    hltdoublemu     = 0;
    hltsinglemu     = 0;
    hltdoubleel     = 0;
    hltsingleel     = 0;

    // Which triggers fired
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
        if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu90  = 1; // MET trigger
        if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu120 = 1; // MET trigger
        if (i == 8  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 9  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 10 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 11 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 12 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 13 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 14 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 15 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 16 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 17 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 18 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 19 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 20 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 21 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165    = 1; // Photon trigger
        if (i == 22 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175    = 1; // Photon trigger
        if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 29 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 30 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 32 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 33 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 34 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 35 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 36 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 37 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 38 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 39 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 40 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 41 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 42 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 43 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 44 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 45 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 46 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 47 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 48 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 49 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 50 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 51 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 52 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 53 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 54 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 55 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
    }

    bool triggered = false;
    if (hltmet90        == 1) triggered = true;
    if (hltmet120       == 1) triggered = true;
    if (hltmetwithmu90  == 1) triggered = true;
    if (hltmetwithmu120 == 1) triggered = true;
    if (hltmetwithmu170 == 1) triggered = true;
    if (hltmetwithmu300 == 1) triggered = true;
    if (hltjetmet90     == 1) triggered = true;
    if (hltjetmet120    == 1) triggered = true;
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltdoublemu     == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltdoubleel     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (applyHLTFilter && !triggered) return;

    // MET filter info
    flagcsctight  = 0;
    flaghbhenoise = 0;
    flaghcallaser = 0;
    flagecaltrig  = 0;
    flageebadsc   = 0;
    flagecallaser = 0;
    flagtrkfail   = 0;
    flagtrkpog    = 0;

    // HCAL Noise info
    flaghnoiseloose = 0;
    flaghnoisetight = 0;
    flaghnoisehilvl = 0;
    if (hcalnoiseH->passLooseNoiseFilter()    ) flaghnoiseloose = 1; 
    if (hcalnoiseH->passTightNoiseFilter()    ) flaghnoisetight = 1; 
    if (hcalnoiseH->passHighLevelNoiseFilter()) flaghnoisehilvl = 1; 

    flaghbheloose = (*hbhelooseH ? 1 : 0);
    flaghbhetight = (*hbhetightH ? 1 : 0);

    // Which MET filters passed
    for (size_t i = 0; i < filterPathsVector.size(); i++) {
        if (filterPathsMap[filterPathsVector[i]] == -1) continue;
        if (i == 0  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagcsctight  = 1; // CSCTightHaloFilter
        if (i == 1  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbhenoise = 1; // HBHENoiseFilter
        if (i == 2  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghcallaser = 1; // hcalLaserEventFilter
        if (i == 3  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagecaltrig  = 1; // EcalDeadCellTriggerPrimitiveFilter
        if (i == 4  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flageebadsc   = 1; // eeBadScFilter
        if (i == 5  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagecallaser = 1; // ecalLaserCorrFilter
        if (i == 6  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagtrkfail   = 1; // trackingFailureFilter
        if (i == 7  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagtrkpog    = 1; // trkPOGFilters
    }

    // Trigger paths
    _trig_pass = "";
    _trig_n    = 0;
    TString path="";
    //
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsH);
    for (int iHLT = 0 ; iHLT<static_cast<int>(triggerResultsH->size()); ++iHLT) {	
      if (triggerResultsH->accept (iHLT)) {
	path = TString(triggerNames.triggerName(iHLT));
	_trig_pass += "_%_"+path ;
	_trig_n++ ;
      }
    }

    if(_verbose>4) cout << "- " << _trig_n << "paths : " << _trig_pass << endl;

    // Trigger objects
    if(!H_trg_obj.isValid()) {
      if(_verbose>0) cout << "Missing collection : " << _IT_trg_obj << " ... skip entry !" << endl;
      return;
    }

    string trgColl,trgFiltStr,trgPathsFFStr;
    vector<int> trgIds;
    vector<string> trgFilt, trgPathsFF, trgPathsFT, trgPathsTF, trgPathsTT;
    ///bool isNone, isL3, isLF, isBoth;
    int iObj=0;

    // Loop over trigger objects    
    for (pat::TriggerObjectStandAlone obj : *H_trg_obj) { // note: not "const &" since we want to call unpackPathNames
      
      if(_verbose>3) cout << "-- iteration #" << iObj << endl;
      iObj++ ;

      obj.unpackPathNames(triggerNames);
      
      // pt,eta,phi
      if(_verbose>1) cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << endl;
      //
      _trig_obj_pt.push_back( obj.pt());
      _trig_obj_eta.push_back(obj.eta());
      _trig_obj_phi.push_back(obj.phi());
      
      // Collection
      trgColl = obj.collection();
      if(_verbose>1) cout << "\t   Collection: " << trgColl << endl;
      //
      _trig_obj_col.push_back(trgColl);
      
      // Trigger Filter Ids
      if(_verbose>1) cout << "\t   Type IDs:   ";
      trgIds = obj.filterIds();
      //
      //_trig_obj_ids.push_back(trgIds);
      //
      for (unsigned h = 0; h < trgIds.size(); ++h) {
	if(_verbose>1) cout << " " << trgIds[h] ;
      }
      if(_verbose>1) cout << endl;
      
      // Trigger filters
      if(_verbose>1) cout << "\t   Filters:    ";
      trgFilt = obj.filterLabels();
      for (unsigned h = 0; h < trgFilt.size(); ++h) {
	if(_verbose>1) cout << " " << trgFilt[h];
	trgFiltStr += trgFilt[h]+"_%_" ;
      }
      if(_verbose>1) cout << endl;
      _trig_obj_lab.push_back(trgFiltStr);
      
      // Trigger paths
      /*
      trgPathsFF = obj.pathNames(false,false);
      trgPathsFT = obj.pathNames(false,true);
      trgPathsTF = obj.pathNames(true,false);
      trgPathsTT = obj.pathNames(true,true);
      //
      _trig_obj_path_FT.push_back(concatenate(trgPathsFT));
      _trig_obj_path_TF.push_back(concatenate(trgPathsTF));
      _trig_obj_path_TT.push_back(concatenate(trgPathsTT));
      //
      if(_verbose>1) cout << "\t   Paths (" << trgPathsFF.size()<<"/"<<trgPathsTT.size()<<"):    ";
      //

      // Loop over all associated paths

      for (unsigned h = 0, n = trgPathsFF.size(); h < n; ++h) {
	
	trgPathsFFStr += trgPathsFF[h]+"_%_";      
	
	isNone = obj.hasPathName( trgPathsFF[h], false, false ); 
	isL3   = obj.hasPathName( trgPathsFF[h], false, true ); 
	isLF   = obj.hasPathName( trgPathsFF[h], true, false ); 
	isBoth = obj.hasPathName( trgPathsFF[h], true, true ); 
	
	if(_verbose>1) {
	  cout << "   " << trgPathsFF[h];
	  if (isBoth)  cout << "(L,3)";
	  if (isL3 && !isBoth)  cout << "(*,3)";
	  if (isLF && !isBoth)  cout << "(L,*)";
	  if (isNone && !isBoth && !isL3 && !isLF)  cout << "(*,*)";
	}
      }
      if(_verbose>1) cout << endl;
      //
      _trig_obj_path_FF.push_back(trgPathsFFStr);
      */
      
      // Clear vectors
      trgIds.clear();
      trgFilt.clear(); 
      trgFiltStr="";
      trgColl="";
      trgPathsFFStr="";
      
      // Increment trigger object index
      _trig_obj_n++ ;
    }  
    //////////// END TRIGGER OBJECTS //////////////

    // Pileup info -- Will need to the updated to the Run-II specifications
    nvtx   = verticesH->size();
    puobs  = 0;
    putrue = 0;
    puwgt  = 1.;
    if (uselheweights && genevtInfoH.isValid()) wgt = genevtInfoH->weight();
    else wgt = 1.0;

    if (pileupInfoH.isValid()) {
        for (vector<PileupSummaryInfo>::const_iterator pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
            if (pileupInfo_iter->getBunchCrossing() == 0) {
                puobs  = pileupInfo_iter->getPU_NumInteractions();
                putrue = pileupInfo_iter->getTrueNumInteractions();
            }
        }
    }

    // MET information 
    hmet    = 0.;
    amet    = 0.;
    bmet    = 0.;
    cmet    = 0.;
    emet    = 0.;
    mmet    = 0.;
    pmet    = 0.;
    omet    = 0.;

    hmetphi = 0.;
    ametphi = 0.;
    bmetphi = 0.;
    cmetphi = 0.;
    emetphi = 0.;
    mmetphi = 0.;
    pmetphi = 0.;
    ometphi = 0.;

    t1pfmet        = t1pfmetH->front().et();
    t1pfmetphi     = t1pfmetH->front().phi();

    calomet        = t1pfmetH->front().caloMETPt();
    calometphi     = t1pfmetH->front().caloMETPhi();

    pfmet          = t1pfmetH->front().uncorrectedPt();
    pfmetphi       = t1pfmetH->front().uncorrectedPhi();

    if (partmetH.isValid() && partmetH->size() == 9) {
    hmet           = (*partmetH)[1].et();
    hmetphi        = (*partmetH)[1].phi();

    amet           = (*partmetH)[2].et();
    ametphi        = (*partmetH)[2].phi();

    bmet           = (*partmetH)[3].et();
    bmetphi        = (*partmetH)[3].phi();

    cmet           = (*partmetH)[4].et();
    cmetphi        = (*partmetH)[4].phi();

    emet           = (*partmetH)[5].et();
    emetphi        = (*partmetH)[5].phi();

    mmet           = (*partmetH)[6].et();
    mmetphi        = (*partmetH)[6].phi();

    pmet           = (*partmetH)[7].et();
    pmetphi        = (*partmetH)[7].phi();

    omet           = (*partmetH)[8].et();
    ometphi        = (*partmetH)[8].phi();
    }

    mumet          = mumetH->front().et();
    mumetphi       = mumetH->front().phi();
 
    t1mumet        = t1mumetH->front().et();
    t1mumetphi     = t1mumetH->front().phi();
 
    pfmupt         = pfmuptH->front().et();
    pfmuphi        = pfmuptH->front().phi();

    // Jet information
    int hardestPhotonIndex = -1;
    double hardestPhotonPt = 0.0;
    for (size_t i = 0; i < tightphotons.size(); i++) {
        if (tightphotons[i]->pt() > hardestPhotonPt) {
            hardestPhotonIndex = i;
            hardestPhotonPt = tightphotons[i]->pt();
        }
    }

    vector<pat::JetRef> incjets;
    vector<pat::JetRef> jets;
    for (vector<pat::Jet>::const_iterator jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
        if (fabs(jets_iter->eta()) > 4.5) continue;
        bool skipjet = false;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < photons.size(); j++) {
            if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
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
        pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());
        incjets.push_back(jetref);
    }

    vector<pat::JetRef> fatjets;
    for (vector<pat::Jet>::const_iterator fatjets_iter = fatjetsH->begin(); fatjets_iter != fatjetsH->end(); ++fatjets_iter) {
        if (fabs(fatjets_iter->eta()) > 2.5) continue;
        bool skipfatjet = false;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), fatjets_iter->eta(), fatjets_iter->phi()) < 0.4) skipfatjet = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), fatjets_iter->eta(), fatjets_iter->phi()) < 0.4) skipfatjet = true;
        }
        for (std::size_t j = 0; j < photons.size(); j++) {
            if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), fatjets_iter->eta(), fatjets_iter->phi()) < 0.4) skipfatjet = true;
        }
        if (skipfatjet) continue;
        bool passfatjetid = false;
        if (fatjets_iter->neutralHadronEnergyFraction()<0.99 && fatjets_iter->neutralEmEnergyFraction()<0.99 && (fatjets_iter->chargedMultiplicity() + fatjets_iter->neutralMultiplicity()) > 1 && fatjets_iter->muonEnergyFraction()<0.8) {
            if (fabs(fatjets_iter->eta()) > 2.4) passfatjetid = true;
            else if (fabs(fatjets_iter->eta()) <= 2.4 && fatjets_iter->chargedHadronEnergyFraction() > 0. && fatjets_iter->chargedEmEnergyFraction() < 0.99 && fatjets_iter->chargedMultiplicity() > 0) passfatjetid = true;
        }
        if (!passfatjetid) continue;
        bool passpuid = false;
        double puidval = fatjets_iter->userFloat("pileupJetId:fullDiscriminant");
        double fatjetabseta = fabs(fatjets_iter->eta());
        if (fatjetabseta >= 0.00 && fatjetabseta < 2.50 && puidval > -0.63) passpuid = true;
        if (fatjetabseta >= 2.50 && fatjetabseta < 2.75 && puidval > -0.60) passpuid = true;
        if (fatjetabseta >= 2.75 && fatjetabseta < 3.00 && puidval > -0.55) passpuid = true;
        if (fatjetabseta >= 3.00 && fatjetabseta < 5.00 && puidval > -0.45) passpuid = true;
        if (!passpuid) continue;
        pat::JetRef fatjetref(fatjetsH, fatjets_iter - fatjetsH->begin());
        fatjets.push_back(fatjetref);
    }


    for (size_t i = 0; i < incjets.size(); i++) {
        if (fabs(incjets[i]->eta()) <= 2.5) jets.push_back(incjets[i]);
    }        

    sort(jets.begin(), jets.end(), jetsorter);
    sort(incjets.begin(), incjets.end(), jetsorter);
    sort(fatjets.begin(), fatjets.end(), jetsorter);

    // AK4 Jets
    njets     = 0;
    nbjets    = 0;
    _jet_n = 0;
    for (size_t i = 0; i < jets.size(); i++) {

        if (jets[i]->pt() > 30) { 
	  njets++;
	  if(_jet_n <= _nMaxJets) {
	    _jet_n++ ;
	    _jet_pt     [i] = jets[i]->pt();
	    _jet_eta    [i] = jets[i]->eta();
	    _jet_phi    [i] = jets[i]->phi();
	    _jet_CHfrac [i] = jets[i]->chargedHadronEnergyFraction();
	    _jet_NHfrac [i] = jets[i]->neutralHadronEnergyFraction();
	    _jet_EMfrac [i] = jets[i]->neutralEmEnergyFraction();
	    _jet_CEMfrac[i] = jets[i]->chargedEmEnergyFraction();
	  }
	}

        if (jets[i]->pt() > 30 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.814) nbjets++;
    }

    signaljetpt      = 0.0;
    signaljeteta     = 0.0;
    signaljetphi     = 0.0;
    signaljetbtag    = 0.0;
    signaljetCHfrac  = 0.0;
    signaljetNHfrac  = 0.0;
    signaljetEMfrac  = 0.0;
    signaljetCEMfrac = 0.0;
    signaljetmetdphi = 0.0;
    signaljetqgl     = 0.0;
    signaljetqgs2    = 0.0;
    signaljetqgmult  = 0;
    signaljetqgptd   = 0.0;
    secondjetpt      = 0.0;
    secondjeteta     = 0.0;
    secondjetphi     = 0.0;
    secondjetbtag    = 0.0;
    secondjetCHfrac  = 0.0;
    secondjetNHfrac  = 0.0;
    secondjetEMfrac  = 0.0;
    secondjetCEMfrac = 0.0;
    secondjetmetdphi = 0.0;
    secondjetqgl     = 0.0;
    secondjetqgs2    = 0.0;
    secondjetqgmult  = 0;
    secondjetqgptd   = 0.0;
    thirdjetpt       = 0.0;
    thirdjeteta      = 0.0;
    thirdjetphi      = 0.0;
    thirdjetbtag     = 0.0;
    thirdjetCHfrac   = 0.0;
    thirdjetNHfrac   = 0.0;
    thirdjetEMfrac   = 0.0;
    thirdjetCEMfrac  = 0.0;
    thirdjetmetdphi  = 0.0;
    thirdjetqgl      = 0.0;
    thirdjetqgs2     = 0.0;
    thirdjetqgmult   = 0;
    thirdjetqgptd    = 0.0;
    jetjetdphi       = 0.0;
    jetmetdphimin    = 0.0;
    incjetmetdphimin = 0.0;


    if (njets > 0) {
        signaljetpt      = jets[0]->pt();
        signaljeteta     = jets[0]->eta();
        signaljetphi     = jets[0]->phi();
        signaljetbtag    = jets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        signaljetCHfrac  = jets[0]->chargedHadronEnergyFraction();
        signaljetNHfrac  = jets[0]->neutralHadronEnergyFraction();
        signaljetEMfrac  = jets[0]->neutralEmEnergyFraction();
        signaljetCEMfrac = jets[0]->chargedEmEnergyFraction();
        signaljetqgl     = (*qglH)   [jets[0]];
        signaljetqgs2    = (*qgs2H)  [jets[0]];
        signaljetqgmult  = (*qgmultH)[jets[0]];
        signaljetqgptd   = (*qgptdH) [jets[0]];
    }

    if (njets > 1) {
        secondjetpt      = jets[1]->pt();
        secondjeteta     = jets[1]->eta();
        secondjetphi     = jets[1]->phi();
        secondjetbtag    = jets[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        secondjetCHfrac  = jets[1]->chargedHadronEnergyFraction();
        secondjetNHfrac  = jets[1]->neutralHadronEnergyFraction();
        secondjetEMfrac  = jets[1]->neutralEmEnergyFraction();
        secondjetCEMfrac = jets[1]->chargedEmEnergyFraction();
        secondjetqgl     = (*qglH)   [jets[1]];
        secondjetqgs2    = (*qgs2H)  [jets[1]];
        secondjetqgmult  = (*qgmultH)[jets[1]];
        secondjetqgptd   = (*qgptdH) [jets[1]];
    }

    if (njets > 2) {
        thirdjetpt       = jets[2]->pt();
        thirdjeteta      = jets[2]->eta();
        thirdjetphi      = jets[2]->phi();
        thirdjetbtag     = jets[2]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        thirdjetCHfrac   = jets[2]->chargedHadronEnergyFraction();
        thirdjetNHfrac   = jets[2]->neutralHadronEnergyFraction();
        thirdjetEMfrac   = jets[2]->neutralEmEnergyFraction();
        thirdjetCEMfrac  = jets[2]->chargedEmEnergyFraction();
        thirdjetqgl      = (*qglH)   [jets[1]];
        thirdjetqgs2     = (*qgs2H)  [jets[1]];
        thirdjetqgmult   = (*qgmultH)[jets[1]];
        thirdjetqgptd    = (*qgptdH) [jets[1]];
    }

    if (signaljetpt > 0.0 && secondjetpt > 0.0) jetjetdphi = deltaPhi(signaljetphi, secondjetphi);
    if (signaljetpt > 0.0) signaljetmetdphi = deltaPhi(signaljetphi, mumetphi);
    if (secondjetpt > 0.0) secondjetmetdphi = deltaPhi(secondjetphi, mumetphi);
    if (thirdjetpt  > 0.0) thirdjetmetdphi  = deltaPhi(thirdjetphi , mumetphi);

    std::vector<double> jetmetdphiminvector;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > 30) {
            double jetphi = atan2(sin(jets[i]->phi()), cos(jets[i]->phi()));
            jetmetdphiminvector.push_back(fabs(deltaPhi(jetphi, mumetphi)));
        }
    }
    if (jetmetdphiminvector.size() > 0) jetmetdphimin = *min_element(jetmetdphiminvector.begin(), jetmetdphiminvector.end());

    std::vector<double> incjetmetdphiminvector;
    for (size_t i = 0; i < incjets.size(); i++) {
        if (incjets[i]->pt() > 30) {
            double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
            incjetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, mumetphi)));
        }
    }
    if (incjetmetdphiminvector.size() > 0) incjetmetdphimin = *min_element(incjetmetdphiminvector.begin(), incjetmetdphiminvector.end());

    // Fat jets
    nfatjets = 0;
    for (size_t i = 0; i < fatjets.size(); i++) {
        if (fatjets[i]->pt() > 100) nfatjets++;
    }

    fatjetpt         = 0.0;
    fatjeteta        = 0.0;
    fatjetphi        = 0.0;
    fatjettau2       = -1.0;
    fatjettau1       = 1.0;
    fatjetprmass     = 0.0;
    fatjetsdmass     = 0.0;
    fatjettrmass     = 0.0;
    fatjetftmass     = 0.0;
    fatjetCHfrac     = 0.0;
    fatjetNHfrac     = 0.0;
    fatjetEMfrac     = 0.0;
    fatjetCEMfrac    = 0.0;
    fatjetmetdphi    = 0.0;

    if (fatjets.size() > 0) {
        fatjetpt         = fatjets[0]->pt();
        fatjeteta        = fatjets[0]->eta();
        fatjetphi        = fatjets[0]->phi();
        fatjettau2       = fatjets[0]->userFloat("NjettinessAK8:tau2");
        fatjettau1       = fatjets[0]->userFloat("NjettinessAK8:tau1");
        fatjetprmass     = fatjets[0]->userFloat("ak8PFJetsCHSPrunedMass");
        fatjetsdmass     = fatjets[0]->userFloat("ak8PFJetsCHSSoftDropMass");
        fatjettrmass     = fatjets[0]->userFloat("ak8PFJetsCHSTrimmedMass");
        fatjetftmass     = fatjets[0]->userFloat("ak8PFJetsCHSFilteredMass");
        fatjetCHfrac     = fatjets[0]->chargedHadronEnergyFraction();
        fatjetNHfrac     = fatjets[0]->neutralHadronEnergyFraction();
        fatjetEMfrac     = fatjets[0]->neutralEmEnergyFraction();
        fatjetCEMfrac    = fatjets[0]->chargedEmEnergyFraction();
        fatjetmetdphi    = deltaPhi(fatjetphi, mumetphi);;
    }

    // QCD suppression handles
    ht     = 0.;
    dht    = -1.;
    mht    = 0.;
    alphat = -1.;

    double mhtx = 0.;
    double mhty = 0.;
    std::vector<double> jetEts;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > 30) {
            ht += jets[i]->pt();
            mhtx -= jets[i]->pt() * cos(jets[i]->phi());
            mhty -= jets[i]->pt() * sin(jets[i]->phi());
            jetEts.push_back(jets[i]->pt());
        }
    }

    mht = sqrt(mhtx*mhtx + mhty*mhty);

    /*
    if (jetEts.size() > 1 && jetEts.size() < 15) { // Memory consumption explodes with large number of jets -- this should be addressed
        // This code is ripped off from UserCode/SusyAnalysis/HadronicSUSYOverlapExercise/ANALYSIS/src 
        std::vector<double> diff( 1<<(jetEts.size()-1) , 0. );
        for(size_t i = 0; i < diff.size(); i++) {
            for(size_t j = 0; j < jetEts.size(); j++) diff[i] += jetEts[j] * ( 1 - 2 * (int(i>>j)&1) );
        }        
        for(size_t i = 0; i < diff.size(); i++) {
            diff[i] = fabs(diff[i]);
        }        
        dht = *min_element(diff.begin(), diff.end());
        alphat = 0.5 * (ht - dht) / sqrt(ht*ht - mht*mht);
    }
    else alphat = 0.0;
    */

    apcjetmetmax = 0.0;
    apcjetmetmin = 0.0;

    std::vector<double> apcjetmetvector;
    for (size_t j = 0; j < jets.size(); j++) {
        if (jets[j]->pt() > 30) {
            apcjetmetvector.push_back(0.);
            for (size_t i = 0; i < jets.size(); i++) {
                if (jets[i]->pt() > 30) {
                    double dphijet = fabs(deltaPhi(jets[i]->phi(), jets[j]->phi()));
                    double jetphi  = atan2(sin(jets[i]->phi()), cos(jets[i]->phi()));
                    double dphimet = fabs(deltaPhi(jetphi, mumetphi));

                    apcjetmetvector.back() += jets[i]->pt() * cos(dphijet/2.0) * sin(dphimet/2.0);
                }
            }
        }
    }
    if (apcjetmetvector.size() > 0) {
        apcjetmetmax = *max_element(apcjetmetvector.begin(), apcjetmetvector.end());
        apcjetmetmin = *min_element(apcjetmetvector.begin(), apcjetmetvector.end());
    }    

    if (ht != 0) {
        apcjetmetmax /= ht;
        apcjetmetmin /= ht;
    }

    // Lepton counts
    nmuons          = muonsH->size();
    nelectrons      = electronsH->size();
    ntightmuons     = tightmuonsH->size();
    ntightelectrons = tightelectronsH->size();
    ntaus           = 0;

    for (vector<pat::Tau>::const_iterator taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
        bool skiptau = false;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
        }
        for (std::size_t j = 0; j < photons.size(); j++) {
            if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
        }
        if (taus_iter->pt() > 18 && fabs(taus_iter->eta()) < 2.3 && taus_iter->tauID("decayModeFinding") > 0.5 && taus_iter->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 5 && !skiptau) ntaus++;
    }

    vector<pat::MuonRef> muonvector;
    for (size_t i = 0; i < muons.size(); i++) muonvector.push_back(muons[i]);

    vector<pat::ElectronRef> electronvector;
    for (size_t i = 0; i < electrons.size(); i++) electronvector.push_back(electrons[i]);

    // W, Z control sample information
    zmass       = 0.0; 
    zpt         = 0.0;
    zeta        = 0.0; 
    zphi        = 0.0;
    zeemass     = 0.0; 
    zeept       = 0.0;
    zeeeta      = 0.0; 
    zeephi      = 0.0;
    wmt         = 0.0;
    wemt        = 0.0;
    emumass     = 0.0;
    emupt       = 0.0;
    emueta      = 0.0;
    emuphi      = 0.0;
    mu1pid      = 0;
    mu1pt       = 0.0;
    mu1eta      = 0.0;
    mu1phi      = 0.0;
    mu1pfpt     = 0.0;
    mu1pfeta    = 0.0;
    mu1pfphi    = 0.0;
    mu1id       = 0;
    mu2pid      = 0;
    mu2pt       = 0.0;
    mu2eta      = 0.0; 
    mu2phi      = 0.0;
    mu2pfpt     = 0.0;
    mu2pfeta    = 0.0;
    mu2pfphi    = 0.0;
    mu2id       = 0;
    el1pid      = 0;
    el1pt       = 0.0;
    el1eta      = 0.0;
    el1phi      = 0.0;
    el1id       = 0;
    el2pid      = 0;
    el2pt       = 0.0;
    el2eta      = 0.0;
    el2phi      = 0.0;
    el2id       = 0;

    sort(muonvector.begin(), muonvector.end(), muonsorter);
    sort(electronvector.begin(), electronvector.end(), electronsorter);

    if (nmuons == 1 || nmuons == 2) {
        pat::MuonRef muon = muons[0];
        mu1pid   = muon->pdgId(); 
        mu1pt    = muon->pt(); 
        mu1eta   = muon->eta(); 
        mu1phi   = muon->phi();
        mu1pfpt  = muon->pfP4().Pt();
        mu1pfeta = muon->pfP4().Eta();
        mu1pfphi = muon->pfP4().Phi();

        for (std::size_t i = 0; i < tightmuons.size(); i++) {
            if (muon == tightmuons[i]) mu1id = 1;
        }

        if (nmuons == 1) wmt = sqrt(2.0 * mu1pt * t1pfmet * (1.0 - cos(deltaPhi(mu1phi, t1pfmetphi))));
    }
   

 
    if (nmuons == 2) {        
        pat::MuonRef muon = muons[1];
        mu2pid = muon->pdgId(); 
        mu2pt  = muon->pt(); 
        mu2eta = muon->eta(); 
        mu2phi = muon->phi();
        mu2pfpt  = muon->pfP4().Pt();
        mu2pfeta = muon->pfP4().Eta();
        mu2pfphi = muon->pfP4().Phi();
    
        for (std::size_t i = 0; i < tightmuons.size(); i++) {
            if (muon == tightmuons[i]) mu2id = 1;
        }
    
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
        TLorentzVector mu2vec; mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muon->p());
    
        TLorentzVector zvec(mu1vec);
        zvec += mu2vec;
    
        zmass = zvec.M();
        zpt   = zvec.Pt();
        zeta  = zvec.Eta();            
        zphi  = zvec.Phi();
    }

    if (nelectrons == 1 || nelectrons == 2) {
        pat::ElectronRef electron = electrons[0];
        el1pid = electron->pdgId();
        el1pt  = electron->pt();
        el1eta = electron->eta();
        el1phi = electron->phi();
        
        for (std::size_t i = 0; i < tightelectrons.size(); i++) {
            if (electron == tightelectrons[i]) el1id = 1;
        }
        
        if (electrons.size() == 1) wemt = sqrt(2.0 * el1pt * t1pfmet * (1.0 - cos(deltaPhi(el1phi, t1pfmetphi))));
    }

    if (nelectrons == 2) {
        pat::ElectronRef electron = electrons[1];
        el2pid = electron->pdgId();
        el2pt  = electron->pt();
        el2eta = electron->eta();
        el2phi = electron->phi();

        for (std::size_t i = 0; i < tightelectrons.size(); i++) {
            if (electron == tightelectrons[i]) el2id = 1;
        }

        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
        TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electron->p());

        TLorentzVector zvec(el1vec);
        zvec += el2vec;

        zeemass = zvec.M();
        zeept   = zvec.Pt();
        zeeeta  = zvec.Eta();
        zeephi  = zvec.Phi();
    }

    if (nmuons == 1 && nelectrons == 1) {
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
        
        TLorentzVector emuvec(mu1vec);
        emuvec += el1vec;
        
        emumass = emuvec.M();
        emupt   = emuvec.Pt();
        emueta  = emuvec.Eta();
        emuphi  = emuvec.Phi();
    } 

    // Photon information
    phidm    = 0;
    phidt    = 0;
    phpt     = 0.0;
    pheta    = 0.0;
    phphi    = 0.0;
    nphotons = photonsH->size();

    if (hardestPhotonIndex >= 0) {
        phidm   = ((*photonMediumIdH)[tightphotons[hardestPhotonIndex]] ? 1 : 0);
        phidt   = ((*photonTightIdH )[tightphotons[hardestPhotonIndex]] ? 1 : 0);
        phpt    = tightphotons[hardestPhotonIndex]->pt();
        pheta   = tightphotons[hardestPhotonIndex]->eta();
        phphi   = tightphotons[hardestPhotonIndex]->phi();
    }

    // Generator-level information
    wzid          = 0;
    wzmass        = 0.0;
    wzpt          = 0.0;
    wzeta         = 0.0;
    wzphi         = 0.0;
    l1id          = 0;
    l1pt          = 0.0;
    l1eta         = 0.0;
    l1phi         = 0.0;
    l2id          = 0;
    l2pt          = 0.0;
    l2eta         = 0.0;
    l2phi         = 0.0;
    parid         = 0;
    parpt         = 0.0;
    pareta        = 0.0;
    parphi        = 0.0;
    ancid         = 0;
    ancpt         = 0.0;
    anceta        = 0.0;
    ancphi        = 0.0;

    if (isWorZMCSample && gensH.isValid()) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() > 1 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) {
                wzid   = gens_iter->pdgId();
                wzmass = gens_iter->mass();
                wzpt   = gens_iter->pt();
                wzeta  = gens_iter->eta();
                wzphi  = gens_iter->phi();
                l1id   = gens_iter->daughter(0)->pdgId();
                l1pt   = gens_iter->daughter(0)->pt();
                l1eta  = gens_iter->daughter(0)->eta();
                l1phi  = gens_iter->daughter(0)->phi();
                l2id   = gens_iter->daughter(1)->pdgId();
                l2pt   = gens_iter->daughter(1)->pt();
                l2eta  = gens_iter->daughter(1)->eta();
                l2phi  = gens_iter->daughter(1)->phi();
                wzmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
            }
        }

        if (wzid == 0) {
            for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
                if (gens_iter->pdgId() == 22 && gens_iter->status() == 1 && gens_iter->isPromptFinalState() && gens_iter->pt() > wzpt) {
                    wzid   = gens_iter->pdgId();
                    wzpt   = gens_iter->pt();
                    wzeta  = gens_iter->eta();
                    wzphi  = gens_iter->phi();

                    findFirstNonPhotonMother(&(*gens_iter), ancid, ancpt, anceta, ancphi);
                    findMother(&(*gens_iter), parid, parpt, pareta, parphi);
                }
            }
        }
    }

    tree->Fill();

    if(_verbose>3) std::cout << "-- tree is filled" << std::endl;
}


void MonoJetTreeMaker::beginJob() {

    if(_verbose>3) std::cout << "- beginJob() starts" << std::endl;

    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"       , "tree");
    int buffersize = 1000; // trig_obj vectors
    // Run, Lumi, Event info
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
    // Event weights
    tree->Branch("xsec"                 , &xsec                 , "xsec/D");
    tree->Branch("wgt"                  , &wgt                  , "wgt/D");
    tree->Branch("kfact"                , &kfact                , "kfact/D");
    tree->Branch("puwgt"                , &puwgt                , "puwgt/D");
    // Pileup info
    tree->Branch("puobs"                , &puobs                , "puobs/I");
    tree->Branch("putrue"               , &putrue               , "putrue/I");
    tree->Branch("nvtx"                 , &nvtx                 , "nvtx/i");

    // Triggers
    tree->Branch("hltmet90"             , &hltmet90             , "hltmet90/b");
    tree->Branch("hltmet120"            , &hltmet120            , "hltmet120/b");
    tree->Branch("hltmetwithmu90"       , &hltmetwithmu90       , "hltmetwithmu90/b");
    tree->Branch("hltmetwithmu120"      , &hltmetwithmu120      , "hltmetwithmu120/b");
    tree->Branch("hltmetwithmu170"      , &hltmetwithmu170      , "hltmetwithmu170/b");
    tree->Branch("hltmetwithmu300"      , &hltmetwithmu300      , "hltmetwithmu300/b");
    tree->Branch("hltjetmet90"          , &hltjetmet90          , "hltjetmet90/b");
    tree->Branch("hltjetmet120"         , &hltjetmet120         , "hltjetmet120/b");
    tree->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/b");
    tree->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/b");
    tree->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/b");
    tree->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/b");
    tree->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/b");
    tree->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/b");

    // Trigger objects
    //tree->Branch("trig_pass",&_trig_pass);
    tree->Branch("trig_n",&_trig_n,"trig_n/I");  
    //
    tree->Branch("trig_obj_n",&_trig_obj_n,"trig_obj_n/I");
    tree->Branch("trig_obj_pt","std::vector<double>",&_trig_obj_pt,buffersize);
    tree->Branch("trig_obj_eta","std::vector<double>",&_trig_obj_eta,buffersize);
    tree->Branch("trig_obj_phi","std::vector<double>",&_trig_obj_phi,buffersize);
    tree->Branch("trig_obj_col","std::vector<std::string>",&_trig_obj_col,buffersize);
    tree->Branch("trig_obj_lab", "std::vector<std::string>",&_trig_obj_lab, buffersize);
    //tree->Branch("trig_obj_ids","std::vector<std::vector<std::int>>",&_trig_obj_ids,buffersize);
    //tree->Branch("trig_obj_path_FF","std::vector<std::string>",&_trig_obj_path_FF,buffersize);
    //tree->Branch("trig_obj_path_FT","std::vector<std::string>",&_trig_obj_path_FT,buffersize);
    //tree->Branch("trig_obj_path_TF","std::vector<std::string>",&_trig_obj_path_TF,buffersize);
    //tree->Branch("trig_obj_path_TT","std::vector<std::string>",&_trig_obj_path_TT,buffersize);

    // Jet vectors
    tree->Branch("jet_n" ,&_jet_n );
    tree->Branch("jet_pt" ,&_jet_pt , "jet_pt[jet_n]/D");
    tree->Branch("jet_eta",&_jet_eta, "jet_eta[jet_n]/D");
    tree->Branch("jet_phi",&_jet_phi,"jet_phi[jet_n]/D");
    tree->Branch("jet_CHfrac",&_jet_CHfrac,"jet_CHfrac[jet_n]/D");
    tree->Branch("jet_NHfrac",&_jet_NHfrac,"jet_NHfrac[jet_n]/D");
    tree->Branch("jet_EMfrac",&_jet_EMfrac,"jet_EMfrac[jet_n]/D");
    tree->Branch("jet_CEMfrac",&_jet_CEMfrac,"jet_CEMfrac[jet_n]/D");

    // MET filters
    tree->Branch("flagcsctight"         , &flagcsctight         , "flagcsctight/b");
    tree->Branch("flaghbhenoise"        , &flaghbhenoise        , "flaghbhenoise/b");
    tree->Branch("flaghbheloose"        , &flaghbheloose        , "flaghbheloose/b");
    tree->Branch("flaghbhetight"        , &flaghbhetight        , "flaghbhetight/b");
    tree->Branch("flaghcallaser"        , &flaghcallaser        , "flaghcallaser/b");
    tree->Branch("flagecaltrig"         , &flagecaltrig         , "flagecaltrig/b");
    tree->Branch("flageebadsc"          , &flageebadsc          , "flageebadsc/b");
    tree->Branch("flagecallaser"        , &flagecallaser        , "flagecallaser/b");
    tree->Branch("flagtrkfail"          , &flagtrkfail          , "flagtrkfail/b");
    tree->Branch("flagtrkpog"           , &flagtrkpog           , "flagtrkpog/b");
    tree->Branch("flaghnoiseloose"      , &flaghnoiseloose      , "flaghnoiseloose/b");
    tree->Branch("flaghnoisetight"      , &flaghnoisetight      , "flaghnoisetight/b");
    tree->Branch("flaghnoisehilvl"      , &flaghnoisehilvl      , "flaghnoisehilvl/b");

    // Object counts
    tree->Branch("nmuons"               , &nmuons               , "nmuons/i");
    tree->Branch("nelectrons"           , &nelectrons           , "nelectrons/i");
    tree->Branch("ntightmuons"          , &ntightmuons          , "ntightmuons/i");
    tree->Branch("ntightelectrons"      , &ntightelectrons      , "ntightelectrons/i");
    tree->Branch("ntaus"                , &ntaus                , "ntaus/i");
    tree->Branch("njets"                , &njets                , "njets/i");
    tree->Branch("nbjets"               , &nbjets               , "nbjets/i");
    tree->Branch("nfatjets"             , &nfatjets             , "nfatjets/i");
    tree->Branch("nphotons"             , &nphotons             , "nphotons/i");
    // MET info
    tree->Branch("pfmet"                , &pfmet                , "pfmet/D");
    tree->Branch("pfmetphi"             , &pfmetphi             , "pfmetphi/D");
    tree->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/D");
    tree->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/D");
    tree->Branch("calomet"              , &calomet              , "calomet/D");
    tree->Branch("calometphi"           , &calometphi           , "calometphi/D");
    tree->Branch("pfmupt"               , &pfmupt               , "pfmupt/D");
    tree->Branch("pfmuphi"              , &pfmuphi              , "pfmuphi/D");
    tree->Branch("mumet"                , &mumet                , "mumet/D");
    tree->Branch("mumetphi"             , &mumetphi             , "mumetphi/D");
    tree->Branch("t1mumet"              , &t1mumet              , "t1mumet/D");
    tree->Branch("t1mumetphi"           , &t1mumetphi           , "t1mumetphi/D");
    tree->Branch("hmet"                 , &hmet                 , "hmet/D");
    tree->Branch("hmetphi"              , &hmetphi              , "hmetphi/D");
    tree->Branch("amet"                 , &amet                 , "amet/D");
    tree->Branch("ametphi"              , &ametphi              , "ametphi/D");
    tree->Branch("bmet"                 , &bmet                 , "bmet/D");
    tree->Branch("bmetphi"              , &bmetphi              , "bmetphi/D");
    tree->Branch("cmet"                 , &cmet                 , "cmet/D");
    tree->Branch("cmetphi"              , &cmetphi              , "cmetphi/D");
    tree->Branch("emet"                 , &emet                 , "emet/D");
    tree->Branch("emetphi"              , &emetphi              , "emetphi/D");
    tree->Branch("mmet"                 , &mmet                 , "mmet/D");
    tree->Branch("mmetphi"              , &mmetphi              , "mmetphi/D");
    tree->Branch("pmet"                 , &pmet                 , "pmet/D");
    tree->Branch("pmetphi"              , &pmetphi              , "pmetphi/D");
    tree->Branch("omet"                 , &omet                 , "omet/D");
    tree->Branch("ometphi"              , &ometphi              , "ometphi/D");
    // Jet info
    tree->Branch("fatjetpt"             , &fatjetpt             , "fatjetpt/D");
    tree->Branch("fatjeteta"            , &fatjeteta            , "fatjeteta/D");
    tree->Branch("fatjetphi"            , &fatjetphi            , "fatjetphi/D");
    tree->Branch("fatjetprmass"         , &fatjetprmass         , "fatjetprmass/D");
    tree->Branch("fatjetsdmass"         , &fatjetsdmass         , "fatjetsdmass/D");
    tree->Branch("fatjettrmass"         , &fatjettrmass         , "fatjettrmass/D");
    tree->Branch("fatjetftmass"         , &fatjetftmass         , "fatjetftmass/D");
    tree->Branch("fatjettau2"           , &fatjettau2           , "fatjettau2/D");
    tree->Branch("fatjettau1"           , &fatjettau1           , "fatjettau1/D");
    tree->Branch("fatjetCHfrac"         , &fatjetCHfrac         , "fatjetCHfrac/D");
    tree->Branch("fatjetNHfrac"         , &fatjetNHfrac         , "fatjetNHfrac/D");
    tree->Branch("fatjetEMfrac"         , &fatjetEMfrac         , "fatjetEMfrac/D");
    tree->Branch("fatjetCEMfrac"        , &fatjetCEMfrac        , "fatjetCEMfrac/D");
    tree->Branch("fatjetmetdphi"        , &fatjetmetdphi        , "fatjetmetdphi/D");
    tree->Branch("signaljetpt"          , &signaljetpt          , "signaljetpt/D");
    tree->Branch("signaljeteta"         , &signaljeteta         , "signaljeteta/D");
    tree->Branch("signaljetphi"         , &signaljetphi         , "signaljetphi/D");
    tree->Branch("signaljetbtag"        , &signaljetbtag        , "signaljetbtag/D");
    tree->Branch("signaljetCHfrac"      , &signaljetCHfrac      , "signaljetCHfrac/D");
    tree->Branch("signaljetNHfrac"      , &signaljetNHfrac      , "signaljetNHfrac/D");
    tree->Branch("signaljetEMfrac"      , &signaljetEMfrac      , "signaljetEMfrac/D");
    tree->Branch("signaljetCEMfrac"     , &signaljetCEMfrac     , "signaljetCEMfrac/D");
    tree->Branch("signaljetmetdphi"     , &signaljetmetdphi     , "signaljetmetdphi/D");
    tree->Branch("signaljetqgl"         , &signaljetqgl         , "signaljetqgl/D");
    tree->Branch("signaljetqgs2"        , &signaljetqgs2        , "signaljetqgs2/D");
    tree->Branch("signaljetqgmult"      , &signaljetqgmult      , "signaljetqgmult/I");
    tree->Branch("signaljetqgptd"       , &signaljetqgptd       , "signaljetqgptd/D");
    tree->Branch("secondjetpt"          , &secondjetpt          , "secondjetpt/D");
    tree->Branch("secondjeteta"         , &secondjeteta         , "secondjeteta/D");
    tree->Branch("secondjetphi"         , &secondjetphi         , "secondjetphi/D");
    tree->Branch("secondjetbtag"        , &secondjetbtag        , "secondjetbtag/D");
    tree->Branch("secondjetCHfrac"      , &secondjetCHfrac      , "secondjetCHfrac/D");
    tree->Branch("secondjetNHfrac"      , &secondjetNHfrac      , "secondjetNHfrac/D");
    tree->Branch("secondjetEMfrac"      , &secondjetEMfrac      , "secondjetEMfrac/D");
    tree->Branch("secondjetCEMfrac"     , &secondjetCEMfrac     , "secondjetCEMfrac/D");
    tree->Branch("secondjetmetdphi"     , &secondjetmetdphi     , "secondjetmetdphi/D");
    tree->Branch("secondjetqgl"         , &secondjetqgl         , "secondjetqgl/D");
    tree->Branch("secondjetqgs2"        , &secondjetqgs2        , "secondjetqgs2/D");
    tree->Branch("secondjetqgmult"      , &secondjetqgmult      , "secondjetqgmult/I");
    tree->Branch("secondjetqgptd"       , &secondjetqgptd       , "secondjetqgptd/D");
    tree->Branch("thirdjetpt"           , &thirdjetpt           , "thirdjetpt/D");
    tree->Branch("thirdjeteta"          , &thirdjeteta          , "thirdjeteta/D");
    tree->Branch("thirdjetphi"          , &thirdjetphi          , "thirdjetphi/D");
    tree->Branch("thirdjetbtag"         , &thirdjetbtag         , "thirdjetbtag/D");
    tree->Branch("thirdjetCHfrac"       , &thirdjetCHfrac       , "thirdjetCHfrac/D");
    tree->Branch("thirdjetNHfrac"       , &thirdjetNHfrac       , "thirdjetNHfrac/D");
    tree->Branch("thirdjetEMfrac"       , &thirdjetEMfrac       , "thirdjetEMfrac/D");
    tree->Branch("thirdjetCEMfrac"      , &thirdjetCEMfrac      , "thirdjetCEMfrac/D");
    tree->Branch("thirdjetmetdphi"      , &thirdjetmetdphi      , "thirdjetmetdphi/D");
    tree->Branch("thirdjetqgl"          , &thirdjetqgl          , "thirdjetqgl/D");
    tree->Branch("thirdjetqgs2"         , &thirdjetqgs2         , "thirdjetqgs2/D");
    tree->Branch("thirdjetqgmult"       , &thirdjetqgmult       , "thirdjetqgmult/I");
    tree->Branch("thirdjetqgptd"        , &thirdjetqgptd        , "thirdjetqgptd/D");
    tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");
    tree->Branch("jetmetdphimin"        , &jetmetdphimin        , "jetmetdphimin/D");
    tree->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/D");
    // QCD suppression
    tree->Branch("ht"                   , &ht                   , "ht/D");
    tree->Branch("dht"                  , &dht                  , "dht/D");
    tree->Branch("mht"                  , &mht                  , "mht/D");
    tree->Branch("alphat"               , &alphat               , "alphat/D");
    tree->Branch("apcjetmetmax"         , &apcjetmetmax         , "apcjetmetmax/D");
    tree->Branch("apcjetmetmin"         , &apcjetmetmin         , "apcjetmetmin/D");
    // Lepton info
    tree->Branch("mu1pid"               , &mu1pid               , "mu1pid/I");
    tree->Branch("mu1pt"                , &mu1pt                , "mu1pt/D");
    tree->Branch("mu1eta"               , &mu1eta               , "mu1eta/D");
    tree->Branch("mu1phi"               , &mu1phi               , "mu1phi/D");
    tree->Branch("mu1pfpt"              , &mu1pfpt              , "mu1pfpt/D");
    tree->Branch("mu1pfeta"             , &mu1pfeta             , "mu1pfeta/D");
    tree->Branch("mu1pfphi"             , &mu1pfphi             , "mu1pfphi/D");
    tree->Branch("mu1id"                , &mu1id                , "mu1id/I");
    tree->Branch("mu2pid"               , &mu2pid               , "mu2pid/I");
    tree->Branch("mu2pt"                , &mu2pt                , "mu2pt/D");
    tree->Branch("mu2eta"               , &mu2eta               , "mu2eta/D");
    tree->Branch("mu2phi"               , &mu2phi               , "mu2phi/D");
    tree->Branch("mu2pfpt"              , &mu2pfpt              , "mu2pfpt/D");
    tree->Branch("mu2pfeta"             , &mu2pfeta             , "mu2pfeta/D");
    tree->Branch("mu2pfphi"             , &mu2pfphi             , "mu2pfphi/D");
    tree->Branch("mu2id"                , &mu2id                , "mu2id/I");
    tree->Branch("el1pid"               , &el1pid               , "el1pid/I");
    tree->Branch("el1pt"                , &el1pt                , "el1pt/D");
    tree->Branch("el1eta"               , &el1eta               , "el1eta/D");
    tree->Branch("el1phi"               , &el1phi               , "el1phi/D");
    tree->Branch("el1id"                , &el1id                , "el1id/I");
    tree->Branch("el2pid"               , &el2pid               , "el2pid/I");
    tree->Branch("el2pt"                , &el2pt                , "el2pt/D");
    tree->Branch("el2eta"               , &el2eta               , "el2eta/D");
    tree->Branch("el2phi"               , &el2phi               , "el2phi/D");
    tree->Branch("el2id"                , &el2id                , "el2id/I");
    // Dilepton info
    tree->Branch("zmass"                , &zmass                , "zmass/D");
    tree->Branch("zpt"                  , &zpt                  , "zpt/D");
    tree->Branch("zeta"                 , &zeta                 , "zeta/D");
    tree->Branch("zphi"                 , &zphi                 , "zphi/D");
    tree->Branch("wmt"                  , &wmt                  , "wmt/D");
    tree->Branch("emumass"              , &emumass              , "emumass/D");
    tree->Branch("emupt"                , &emupt                , "emupt/D");
    tree->Branch("emueta"               , &emueta               , "emueta/D");
    tree->Branch("emuphi"               , &emuphi               , "emuphi/D");
    tree->Branch("zeemass"              , &zeemass              , "zeemass/D");
    tree->Branch("zeept"                , &zeept                , "zeeept/D");
    tree->Branch("zeeeta"               , &zeeeta               , "zeeeta/D");
    tree->Branch("zeephi"               , &zeephi               , "zeephi/D");
    tree->Branch("wemt"                 , &wemt                 , "wemt/D");
    // Photon info
    tree->Branch("phidm"                , &phidm                , "phidm/I");
    tree->Branch("phidt"                , &phidt                , "phidt/I");
    tree->Branch("phpt"                 , &phpt                 , "phpt/D");
    tree->Branch("pheta"                , &pheta                , "pheta/D");
    tree->Branch("phphi"                , &phphi                , "phphi/D");
    // W/Z gen-level info
    tree->Branch("wzid"                 , &wzid                 , "wzid/I");
    tree->Branch("wzmass"               , &wzmass               , "wzmass/D");
    tree->Branch("wzmt"                 , &wzmt                 , "wzmt/D");
    tree->Branch("wzpt"                 , &wzpt                 , "wzpt/D");
    tree->Branch("wzeta"                , &wzeta                , "wzeta/D");
    tree->Branch("wzphi"                , &wzphi                , "wzphi/D");
    tree->Branch("l1id"                 , &l1id                 , "l1id/I");
    tree->Branch("l1pt"                 , &l1pt                 , "l1pt/D");
    tree->Branch("l1eta"                , &l1eta                , "l1eta/D");
    tree->Branch("l1phi"                , &l1phi                , "l1phi/D");
    tree->Branch("l2id"                 , &l2id                 , "l2id/I");
    tree->Branch("l2pt"                 , &l2pt                 , "l2pt/D");
    tree->Branch("l2eta"                , &l2eta                , "l2eta/D");
    tree->Branch("l2phi"                , &l2phi                , "l2phi/D");
    tree->Branch("parid"                , &parid                , "parid/I");
    tree->Branch("parpt"                , &parpt                , "parpt/D");
    tree->Branch("pareta"               , &pareta               , "pareta/D");
    tree->Branch("parphi"               , &parphi               , "parphi/D");
    tree->Branch("ancid"                , &ancid                , "ancid/I");
    tree->Branch("ancpt"                , &ancpt                , "ancpt/D");
    tree->Branch("anceta"               , &anceta               , "anceta/D");
    tree->Branch("ancphi"               , &ancphi               , "ancphi/D");

    if(_verbose>3) std::cout << "- beginJob() ends" << std::endl;
}

void MonoJetTreeMaker::endJob() {
}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_PFMET90_PFMHT90_IDTight");
    triggerPathsVector.push_back("HLT_PFMET120_PFMHT120_IDTight");
    triggerPathsVector.push_back("HLT_PFMET170_NoiseCleaned");
    triggerPathsVector.push_back("HLT_PFMET170_JetIdCleaned");
    triggerPathsVector.push_back("HLT_PFMET170_HBHECleaned");
    triggerPathsVector.push_back("HLT_PFMET170_v");
    triggerPathsVector.push_back("HLT_PFMET300_NoiseCleaned");
    triggerPathsVector.push_back("HLT_PFMET300_JetIdCleaned");
    triggerPathsVector.push_back("HLT_PFMET300_v");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight");
    triggerPathsVector.push_back("HLT_Photon165_HE10");
    triggerPathsVector.push_back("HLT_Photon175");
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
    triggerPathsVector.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"); //ND: Spring15 MC single ele
    triggerPathsVector.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    triggerPathsVector.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    triggerPathsVector.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    triggerPathsVector.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v");
    triggerPathsVector.push_back("HLT_Ele27_WP85_Gsf_v");
    triggerPathsVector.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v");
    triggerPathsVector.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v");
    triggerPathsVector.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
    triggerPathsVector.push_back("HLT_Ele25WP60_SC4_Mass55_v");

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

    // Trigger objects
    _trig_pass = "";
    _trig_n    = 0;
    //
    _trig_obj_n = 0;
    _trig_obj_pt.clear();
    _trig_obj_eta.clear();
    _trig_obj_phi.clear();
    _trig_obj_col.clear();
    _trig_obj_lab.clear();
    //_trig_obj_ids.clear();
    //_trig_obj_path_FF.clear();
    //_trig_obj_path_FT.clear();
    //_trig_obj_path_TF.clear();
    //_trig_obj_path_TT.clear();

    // jets
    for (size_t i = 0; i < _nMaxJets; i++) {
      _jet_pt     [i] = 0; 
      _jet_eta    [i] = 0; 
      _jet_phi    [i] = 0; 
      _jet_CHfrac [i] = 0; 
      _jet_NHfrac [i] = 0; 
      _jet_EMfrac [i] = 0; 
      _jet_CEMfrac[i] = 0; 
    }

    // MET filter Paths
    filterPathsVector.push_back("Flag_CSCTightHaloFilter");
    filterPathsVector.push_back("Flag_HBHENoiseFilter");
    filterPathsVector.push_back("Flag_hcalLaserEventFilter");
    filterPathsVector.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
    filterPathsVector.push_back("Flag_eeBadScFilter");
    filterPathsVector.push_back("Flag_ecalLaserCorrFilter");
    filterPathsVector.push_back("Flag_trackingFailureFilter");
    filterPathsVector.push_back("Flag_trkPOGFilters");

    HLTConfigProvider fltrConfig;
    fltrConfig.init(iRun, iSetup, filterResultsTag.process(), changedConfig);

    for (size_t i = 0; i < filterPathsVector.size(); i++) {
        filterPathsMap[filterPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < filterPathsVector.size(); i++){
        TPRegexp pattern(filterPathsVector[i]);
        for(size_t j = 0; j < fltrConfig.triggerNames().size(); j++){
            std::string pathName = fltrConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                filterPathsMap[filterPathsVector[i]] = j;
            }
        }
    }

}

void MonoJetTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

/*
This code is ripped off from https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/PhotonNtuplerMiniAOD.cc
*/
void MonoJetTreeMaker::findFirstNonPhotonMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {
    if (particle == 0) {
        return;
    }
    if (abs(particle->pdgId()) == 22) {
        findFirstNonPhotonMother(particle->mother(0), ancestorid, ancestorpt, ancestoreta, ancestorphi);
    }
    else {
        ancestorid  = particle->pdgId();
        ancestorpt  = particle->pt();
        ancestoreta = particle->eta();
        ancestorphi = particle->phi();
    }
    return;
}

void MonoJetTreeMaker::findMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {
    if (particle == 0) {
        return;
    }
    if (abs(particle->pdgId()) == 22) {
        ancestorid  = particle->pdgId();
        ancestorpt  = particle->pt();
        ancestoreta = particle->eta();
        ancestorphi = particle->phi();
    }
    return;
}

std::string MonoJetTreeMaker::concatenate(std::vector<std::string> vstring)
{
    std::string result;
    for(UInt_t i=0 ; i<vstring.size() ; i++) {
      result += vstring[i]+"_%_";
    }
    return result;
}


void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

void MonoJetTreeMaker::Init()
{

  _trig_obj_n = 0;
  _trig_obj_pt.clear();
  _trig_obj_eta.clear();
  _trig_obj_phi.clear();
  _trig_obj_col.clear();
  _trig_obj_lab.clear();
  //_trig_obj_ids.clear();
  //_trig_obj_path_FF.clear();
  //_trig_obj_path_FT.clear();
  //_trig_obj_path_TF.clear();
  //_trig_obj_path_TT.clear();

  for (size_t i = 0; i < _nMaxJets; i++) {
    _jet_pt     [i] = 0; 
    _jet_eta    [i] = 0; 
    _jet_phi    [i] = 0; 
    _jet_CHfrac [i] = 0; 
    _jet_NHfrac [i] = 0; 
    _jet_EMfrac [i] = 0; 
    _jet_CEMfrac[i] = 0; 
  }

}

DEFINE_FWK_MODULE(MonoJetTreeMaker);

