#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class METBreakDownProducer : public edm::EDProducer {
    public:
        explicit METBreakDownProducer(const edm::ParameterSet&);
        ~METBreakDownProducer();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        edm::InputTag candsTag;
        edm::EDGetTokenT<edm::View<reco::Candidate> > candsToken;
};

METBreakDownProducer::METBreakDownProducer(const edm::ParameterSet& iConfig): 
    candsTag(iConfig.getParameter<edm::InputTag>("pfcands"))
{
    produces<reco::METCollection>();

    candsToken = consumes<edm::View<reco::Candidate> >(candsTag);
}


METBreakDownProducer::~METBreakDownProducer() {
}

void METBreakDownProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<View<Candidate> > candsH;
    iEvent.getByToken(candsToken, candsH);

    std::auto_ptr<METCollection> output(new METCollection);

    double metx  = 0.0;
    double mety  = 0.0;

    double hmetx = 0.0;
    double hmety = 0.0;

    double ametx = 0.0;
    double amety = 0.0;

    double bmetx = 0.0;
    double bmety = 0.0;

    double cmetx = 0.0;
    double cmety = 0.0;

    double emetx = 0.0;
    double emety = 0.0;

    double mmetx = 0.0;
    double mmety = 0.0;

    double pmetx = 0.0;
    double pmety = 0.0;

    double ometx = 0.0;
    double omety = 0.0;

    for (View<Candidate>::const_iterator cands_iter = candsH->begin(); cands_iter != candsH->end(); ++cands_iter) {
        int absid = abs(cands_iter->pdgId());
        if (absid == 130) {
            hmetx -= cands_iter->pt() * cos(cands_iter->phi());
            hmety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 1) {
            ametx -= cands_iter->pt() * cos(cands_iter->phi());
            amety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 2) {
            bmetx -= cands_iter->pt() * cos(cands_iter->phi());
            bmety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 211) {
            cmetx -= cands_iter->pt() * cos(cands_iter->phi());
            cmety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 11) {
            emetx -= cands_iter->pt() * cos(cands_iter->phi());
            emety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 13) {
            mmetx -= cands_iter->pt() * cos(cands_iter->phi());
            mmety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else if (absid == 22) {
            pmetx -= cands_iter->pt() * cos(cands_iter->phi());
            pmety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        else {
            ometx -= cands_iter->pt() * cos(cands_iter->phi());
            omety -= cands_iter->pt() * sin(cands_iter->phi());
        }
        metx -= cands_iter->pt() * cos(cands_iter->phi());
        mety -= cands_iter->pt() * sin(cands_iter->phi());
    }            
    double  met = sqrt( metx* metx +  mety* mety);
    double hmet = sqrt(hmetx*hmetx + hmety*hmety);
    double amet = sqrt(ametx*ametx + amety*amety);
    double bmet = sqrt(bmetx*bmetx + bmety*bmety);
    double cmet = sqrt(cmetx*cmetx + cmety*cmety);
    double emet = sqrt(emetx*emetx + emety*emety);
    double mmet = sqrt(mmetx*mmetx + mmety*mmety);
    double pmet = sqrt(pmetx*pmetx + pmety*pmety);
    double omet = sqrt(ometx*ometx + omety*omety);

    MET  metcand;
    MET hmetcand;
    MET ametcand;
    MET bmetcand;
    MET cmetcand;
    MET emetcand;
    MET mmetcand;
    MET pmetcand;
    MET ometcand;

    metcand .setP4(reco::Candidate::LorentzVector( metx,  mety, 0.,  met));
    hmetcand.setP4(reco::Candidate::LorentzVector(hmetx, hmety, 0., hmet));
    ametcand.setP4(reco::Candidate::LorentzVector(ametx, amety, 0., amet));
    bmetcand.setP4(reco::Candidate::LorentzVector(bmetx, bmety, 0., bmet));
    cmetcand.setP4(reco::Candidate::LorentzVector(cmetx, cmety, 0., cmet));
    emetcand.setP4(reco::Candidate::LorentzVector(emetx, emety, 0., emet));
    mmetcand.setP4(reco::Candidate::LorentzVector(mmetx, mmety, 0., mmet));
    pmetcand.setP4(reco::Candidate::LorentzVector(pmetx, pmety, 0., pmet));
    ometcand.setP4(reco::Candidate::LorentzVector(ometx, omety, 0., omet));

    output->push_back( metcand);
    output->push_back(hmetcand);
    output->push_back(ametcand);
    output->push_back(bmetcand);
    output->push_back(cmetcand);
    output->push_back(emetcand);
    output->push_back(mmetcand);
    output->push_back(pmetcand);
    output->push_back(ometcand);

    iEvent.put(output);
}

void METBreakDownProducer::beginJob() {
}

void METBreakDownProducer::endJob() {
}

void METBreakDownProducer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void METBreakDownProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void METBreakDownProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void METBreakDownProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void METBreakDownProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(METBreakDownProducer);
