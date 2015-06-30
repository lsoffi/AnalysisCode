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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

# Is this a simulation or real data
isMC = True

# Filter on high MET events
filterHighMETEvents = False

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/04DC7D13-1E0C-E511-847C-00A0D1EE923C.root')
)

# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

# Set the global tag depending on the sample type
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if isMC:
    process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'   # for Simulation
else:
    process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'   # for Data (for now set to the MC tag)

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True)
)

# Electron and Photon ValueMaps for identification
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']
ph_id_modules  = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']
for idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.selectedObjects = cms.EDProducer("PFCleaner",
    vertices = cms.InputTag("goodVertices"),
    pfcands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    electronidveto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    electronidmedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    photonidloose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose")
)

# Define all the METs corrected for lepton/photon momenta
process.mumet = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    muons = cms.InputTag("selectedObjects", "muons"),
    useuncorrmet = cms.bool(True)
)
process.pfmupt = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    muons = cms.InputTag("selectedObjects", "muons"),
    muptonly = cms.bool(True)
)
process.phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
    useuncorrmet = cms.bool(True)
)
process.t1mumet = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    muons = cms.InputTag("selectedObjects", "muons")
)
process.t1phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")) 
)

# Remove MC matching when running on data
if not isMC:
    runOnData(process)

process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    genevt = cms.InputTag("generator"),
    vertices = cms.InputTag("goodVertices"),
    gens = cms.InputTag("prunedGenParticles"),
    muons = cms.InputTag("selectedObjects", "muons"),
    electrons = cms.InputTag("selectedObjects", "electrons"),
    photons = cms.InputTag("selectedObjects", "photons"),
    tightmuons = cms.InputTag("selectedObjects", "tightmuons"),
    tightelectrons = cms.InputTag("selectedObjects", "tightelectrons"),
    tightphotons = cms.InputTag("selectedObjects", "tightphotons"),
    loosephotons = cms.InputTag("selectedObjects", "loosephotons"),
    rndgammaiso = cms.InputTag("selectedObjects", "rndgammaiso"),
    photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    pfmet = cms.InputTag("pfMet"),
    t1pfmet = cms.InputTag("slimmedMETs"),
    pfmupt = cms.InputTag("pfmupt"),
    mumet = cms.InputTag("mumet"),
    phmet = cms.InputTag("phmet"),
    t1mumet = cms.InputTag("t1mumet"),
    t1phmet = cms.InputTag("t1phmet"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    xsec = cms.double(32.293),
    cleanMuonJet = cms.bool(True),
    cleanElectronJet = cms.bool(True),
    cleanPhotonJet = cms.bool(True),
    uselheweights = cms.bool(False),
    isWorZMCSample = cms.bool(False)
)

process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheinfo = cms.InputTag("externalLHEProducer"),
    geninfo = cms.InputTag("generator"),
    uselheweights = cms.bool(False),
    addqcdpdfweights = cms.bool(False)
)

process.metfilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("t1mumet"),
    cut = cms.string("et > 200"),
    filter = cms.bool(True)
)


if filterHighMETEvents: 
    process.treePath = cms.Path(process.gentree + process.goodVertices + process.metfilter + process.tree)
else :
    process.treePath = cms.Path(process.gentree + process.goodVertices + process.tree)

