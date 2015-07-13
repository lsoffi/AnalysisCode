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

# Filter on triggered events
filterOnHLT = False

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/80CF5456-B9EC-E411-93DA-002618FDA248.root')
    #fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValProdTTbar_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/0A9E2CED-C9EC-E411-A8E4-003048FFCBA8.root')
    #fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValZMM_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/1C258197-BEEC-E411-A70A-002618943901.root')
    #fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValWM_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/6AF5B73D-C8EC-E411-8559-002618943865.root')
    #fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValRSGravitonToGaGa_13TeV/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/189277BA-DCEC-E411-B3B8-0025905B859E.root')
    #fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/00000/6408230F-9F08-E511-A1A6-D4AE526A023A.root')
    #fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/60000/04963444-D107-E511-B245-02163E00F339.root')
)

# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("treedump.root"))

# Set the global tag depending on the sample type
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if isMC:
    process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'   # for Simulation
else:
    process.GlobalTag.globaltag = 'GR_P_V54::All'       # for Data

# Set the process for reading MET filter flags from TriggerResults
if isMC:
    metFilterProcess = "PAT"
else :
    metFilterProcess = "RECO"

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(False)
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

# The tree maker
process.tree = cms.EDAnalyzer("TreeDumper",
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    filterResults = cms.InputTag("TriggerResults", "", metFilterProcess),
    hcalnoise = cms.InputTag("hcalnoise"),
    vertices = cms.InputTag("goodVertices"),
    genevt = cms.InputTag("generator"),
    gens = cms.InputTag("prunedGenParticles"),
    genjets = cms.InputTag("slimmedGenJets"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    electronidveto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    electronidloose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    electronidmedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    electronidtight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    photonidloose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
    photonidmedium = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
    photonidtight = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
    photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    t1pfmet = cms.InputTag("slimmedMETs"),
    isVMCSample = cms.bool(False),
    uselheweights = cms.bool(False),
    xsec = cms.double(1.0),
    kfactor = cms.double(1.0),
    applyHLTFilter = cms.bool(filterOnHLT),
    applyHighMETFilter = cms.bool(filterHighMETEvents)
)

# Gen Tree -- For event weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    addqcdpdfweights = cms.bool(False),
    uselheweights = cms.bool(False),
    lheinfo = cms.InputTag("externalLHEProducer"),
    geninfo = cms.InputTag("generator")
)

# Process Path
if isMC:
    process.treePath = cms.Path(process.gentree + process.tree)
else :
    process.treePath = cms.Path(process.tree)

