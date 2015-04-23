import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("TREE")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

# Is this a simulation or real data
isMC = True

# Filter on high MET events
filterHighMETEvents = True

# Is this a dielectron control sample?
isPhotonSample = False

# Define the global tag depending on the sample type
if isMC:
    process.GlobalTag.globaltag = 'START53_V27::All'     # for Simulation
else:
    process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'   # for Jan 22 2013 ReReco

# Define the input source
process.source = cms.Source("PoolSource", 
    #fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/G_Pt-120to170_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/007CC955-3C17-E211-B882-003048C69272.root')
    fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/0012AB48-AA02-E211-ACF7-001E67398043.root')
)

# Define the output -- Needed for some of the PAT sequences
process.out = cms.OutputModule("PoolOutputModule", 
    outputCommands =  cms.untracked.vstring(), 
    fileName = cms.untracked.string('monoxSkim.root') 
)

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# Define rho to be used for electron and photon isolation
from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJets25 = kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)

# Setup the particle flow isolation for electrons and photons (the default values in the GsfElectron object correspond to dR = 0.4)
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import *
from CommonTools.ParticleFlow.Isolation.tools_cfi import *
from RecoParticleFlow.PFProducer.electronPFIsolationValues_cff import *
from RecoParticleFlow.PFProducer.photonPFIsolationValues_cff import *

process.pfAllChargedHadrons = pfAllChargedHadrons.clone()
process.pfAllNeutralHadrons = pfAllNeutralHadrons.clone()
process.pfAllPhotons = pfAllPhotons.clone()

process.elPFIsoDepositCharged = isoDepositReplace('gsfElectrons', 'pfAllChargedHadrons')
process.elPFIsoDepositNeutral = isoDepositReplace('gsfElectrons', 'pfAllNeutralHadrons')
process.elPFIsoDepositGamma   = isoDepositReplace('gsfElectrons', 'pfAllPhotons')
process.elPFIsoDepositPU      = isoDepositReplace('gsfElectrons', 'pfPileUpAllChargedParticles')
process.elPFIsoValueCharged03 = elPFIsoValueCharged03PFId.clone()
process.elPFIsoValueNeutral03 = elPFIsoValueNeutral03PFId.clone()
process.elPFIsoValueGamma03   = elPFIsoValueGamma03PFId.clone()
process.elPFIsoValuePU03      = elPFIsoValuePU03PFId.clone()

process.phPFIsoDepositCharged = isoDepositReplace(cms.InputTag("pfPhotonTranslator", "pfphot"), 'pfAllChargedHadrons')
process.phPFIsoDepositNeutral = isoDepositReplace(cms.InputTag("pfPhotonTranslator", "pfphot"), 'pfAllNeutralHadrons')
process.phPFIsoDepositGamma   = isoDepositReplace(cms.InputTag("pfPhotonTranslator", "pfphot"), 'pfAllPhotons')
process.phPFIsoValueCharged03 = phPFIsoValueCharged03PFId.clone()
process.phPFIsoValueNeutral03 = phPFIsoValueNeutral03PFId.clone()
process.phPFIsoValueGamma03   = phPFIsoValueGamma03PFId.clone()

process.phPFIsoDepositCharged.ExtractorPSet.DR_Veto = cms.double(0)
process.phPFIsoDepositNeutral.ExtractorPSet.DR_Veto = cms.double(0)
process.phPFIsoDepositGamma.ExtractorPSet.ComponentName = cms.string('PFCandWithSuperClusterExtractor')
process.phPFIsoDepositGamma.ExtractorPSet.SCMatch_Veto = cms.bool(True)
process.phPFIsoDepositGamma.ExtractorPSet.MissHitSCMatch_Veto = cms.bool(False)

# Load the PAT sequence
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Define ak5PFJets without leptons and pileup
process.pfPileUp = cms.EDProducer("PFPileUp",
    checkClosestZVertex = cms.bool(True),
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("goodVertices")
)

process.particleFlowClean = cms.EDProducer("PFCleaner",
    src = cms.InputTag("particleFlow"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("goodVertices"),
    pfpileup = cms.InputTag("pfPileUp"),
    rho = cms.InputTag("kt6PFJets", "rho"),
    conversions = cms.InputTag("allConversions"),
    electrons = cms.InputTag("gsfElectrons"),
    photons = cms.InputTag("photons"),
    electronPFIsoCH = cms.InputTag("elPFIsoValueCharged03"),
    electronPFIsoNH = cms.InputTag("elPFIsoValueNeutral03"),
    electronPFIsoPH = cms.InputTag("elPFIsoValueGamma03"),
    electronPFIsoPU = cms.InputTag("elPFIsoValuePU03"),
    photonPFIsoCH = cms.InputTag("phPFIsoValueCharged03"),
    photonPFIsoNH = cms.InputTag("phPFIsoValueNeutral03"),
    photonPFIsoPH = cms.InputTag("phPFIsoValueGamma03"),
    d0cut = cms.double(0.2),
    dzcut = cms.double(0.5),
    vetophotons = cms.bool(False)
)

from RecoJets.JetProducers.ak5PFJets_cfi import *
process.ak5PFJetsClean = ak5PFJets.clone(src = cms.InputTag("particleFlowClean", "pfcands"))

process.ca8PFJetsClean = ak5PFJets.clone(
    src = cms.InputTag("particleFlowClean", "pfcands"),
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
)

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsCleanPruned = ak5PFJetsPruned.clone(
    src = cms.InputTag("particleFlowClean", "pfcands"),
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.Configuration.RecoGenJets_cff import ak7GenJetsNoNu
process.ca8GenJetsNoNu = ak7GenJetsNoNu.clone()
process.ca8GenJetsNoNu.rParam = 0.8
process.ca8GenJetsNoNu.jetAlgorithm = "CambridgeAachen"

process.ca8PFJetsCleanVMaps = cms.EDProducer("JetSubstructureValueMapsProducer",
    src = cms.InputTag("ca8PFJetsClean"),
    jetRadius = cms.double(0.8)
)

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    cms.InputTag('ak5PFJetsClean'),
    doJTA            = False,
    doBTagging       = True,
    jetCorrLabel     = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    doType1MET       = False,
    genJetCollection = cms.InputTag('ak5GenJets'),
    doJetID          = True
)

addJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    algoLabel        = 'AK5',
    typeLabel        = 'PF',
    doJTA            = False,
    doBTagging       = False,
    jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    doType1MET       = False,
    genJetCollection = cms.InputTag('ak5GenJets'),
    doJetID          = True
)   

addJetCollection(
    process,
    cms.InputTag('ca8PFJetsClean'),
    algoLabel        = 'CA8',
    typeLabel        = 'PF',
    doJTA            = False,
    doBTagging       = False,
    jetCorrLabel     = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    doType1MET       = False,
    genJetCollection = cms.InputTag('ca8GenJetsNoNu'),
    doJetID          = True
)

addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCleanPruned'),
    algoLabel        = 'CA8Pruned',
    typeLabel        = 'PF',
    doJTA            = False,
    doBTagging       = False,
    jetCorrLabel     = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    doType1MET       = False,
    genJetCollection = cms.InputTag('ca8GenJetsNoNu'),
    doJetID          = True
)

process.patJetsCA8PF.userData.userFloats.src = cms.VInputTag("ca8PFJetsCleanVMaps:tau1", "ca8PFJetsCleanVMaps:tau2", "ca8PFJetsCleanVMaps:tau3")

process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
)
if isMC:
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
else:
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

# Use HPS taus
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.selectedPatTaus.cut = "pt > 20 && abs(eta) < 2.3 && tauID('decayModeFinding') > 0.5 && tauID('byLooseCombinedIsolationDeltaBetaCorr') > 0.5 && tauID('againstMuonTight2') > 0.5 && tauID('againstElectronLoose') > 0.5"

# Define all the METs corrected for lepton/photon momenta
process.mumet = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("pfMet"),
    muons = cms.InputTag("particleFlowClean", "muons")
)

process.phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("pfMet"),
    cands = cms.VInputTag(cms.InputTag("particleFlowClean", "leadingphoton")) 
)

process.t1mumet = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("pfType1CorrectedMet"),
    muons = cms.InputTag("particleFlowClean", "muons")
)

process.t1phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("pfType1CorrectedMet"),
    cands = cms.VInputTag(cms.InputTag("particleFlowClean", "leadingphoton"))
)

process.pfmupt = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("pfMet"),
    muons = cms.InputTag("particleFlowClean", "muons"),
    muptonly = cms.bool(True)
)

# Remove MC matching when running on data
if not isMC:
    runOnData(process)

# Use good primary vertices in the PAT sequence
massSearchReplaceAnyInputTag(process.patDefaultSequence, cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodVertices"), True)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    vertices = cms.InputTag("goodVertices"),
    muons = cms.InputTag("particleFlowClean", "muons"),
    electrons = cms.InputTag("particleFlowClean", "electrons"),
    photons = cms.InputTag("particleFlowClean", "photons"),
    tightmuons = cms.InputTag("particleFlowClean", "tightmuons"),
    tightelectrons = cms.InputTag("particleFlowClean", "tightelectrons"),
    tightphotons = cms.InputTag("particleFlowClean", "tightphotons"),
    taus = cms.InputTag("selectedPatTaus"),
    jets = cms.InputTag("patJets"),
    nochsjets = cms.InputTag("patJetsAK5PF"),
    fatjets = cms.InputTag("patJetsCA8PF"),
    prunedfatjets = cms.InputTag("patJetsCA8PrunedPF"),
    gens = cms.InputTag("genParticles"),
    pfmet = cms.InputTag("pfMet"),
    t1pfmet = cms.InputTag("pfType1CorrectedMet"),
    calomet = cms.InputTag("met"),
    pfmupt = cms.InputTag("pfmupt"),
    mumet = cms.InputTag("mumet"),
    phmet = cms.InputTag("phmet"),
    t1mumet = cms.InputTag("t1mumet"),
    t1phmet = cms.InputTag("t1phmet"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    weight = cms.double(1000.0*225.2/6923750.0),
    isWorZMCSample = cms.bool(False),
    cleanPhotonJet = cms.bool(False)
)

# Select events passing the monojet triggers and having MET > 200 GeV
if isPhotonSample :
    process.triggerfilter = cms.EDFilter("HLTCheckFilter",
        triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        triggerPaths   = cms.vstring('HLT_Photon135', 'HLT_Photon150')   
    )

    process.metfilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("phmet"),
        cut = cms.string("et > 200"),
        filter = cms.bool(True)
    )
else :
    process.triggerfilter = cms.EDFilter("HLTCheckFilter",
        triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        triggerPaths   = cms.vstring('HLT_MET120_HBHENoiseCleaned', 'HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95', 'HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95', 'HLT_Mu17_Mu8', 'HLT_Mu17_TkMu8')   
        #triggerPaths   = cms.vstring('HLT_Mu17_Mu8', 'HLT_Mu17_TkMu8')   
    )

    process.metfilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("mumet"),
        cut = cms.string("et > 200"),
        filter = cms.bool(True)
    )


if filterHighMETEvents: 
    #process.treePath = cms.Path(process.goodVertices + process.triggerfilter + process.metfilter + process.tree)
    process.treePath = cms.Path(process.goodVertices + process.triggerfilter + process.tree)
else :
    process.treePath = cms.Path(process.goodVertices + process.tree)

