import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("TREE")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
isMC = False

# Filter on high MET events
filterHighMETEvents = False

# Filter on triggered events
filterOnHLT = True

# Redo jets and MET with updated JEC
redoJetsMET = False

# Use private JECs since the GTs are not updated
usePrivateSQlite = False

# Apply L2L3 residual corrections
applyL2L3Residuals = False

# Process name used in MiniAOD -- needed to get the correct trigger results, and also for redoing the MET
miniAODProcess = "PAT"

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring([
        '/store/data/Run2015D/DoubleMuon/MINIAOD/05Oct2015-v1/30000/B0E177F6-8A6F-E511-A81F-0025905B85D6.root'
    ])
)

# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
if isMC:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'   # for Simulation
else:
    process.GlobalTag.globaltag = '74X_dataRun2_v4'            # for Data

# Setup the private SQLite -- Ripped from PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
if usePrivateSQlite:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    era = "Summer15_25nsV5"
    if isMC : 
        era += "_MC"
    else :
        era += "_DATA"
    dBFile = os.path.expandvars(era+".db")
    process.jec = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBFile),
        toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
            ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
            ),
        )
    )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')


# JEC levels when redoing jets and MET
if isMC:
    JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else :
    if not applyL2L3Residuals : 
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    else : 
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

# Re-run the HBHE Noise filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

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

ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
ph_id_modules  = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff']
for idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


# Re-running of jets and MET
jetCollName = "slimmedJets"

if redoJetsMET :  
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

    runMetCorAndUncFromMiniAOD(process,
        isData = (not isMC),
    )
    if miniAODProcess != "PAT" :
        if isMC : 
            process.genMetExtractor.metSource= cms.InputTag("slimmedMETs", "", miniAODProcess)   
        process.slimmedMETs.t01Variation = cms.InputTag("slimmedMETs", "", miniAODProcess) 

    if not applyL2L3Residuals : 
        process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJets"),
        levels = JECLevels,
        payload = 'AK4PFchs' 
    )

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.slimmedJetsRecorrected = patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
    )
    jetCollName = "slimmedJetsRecorrected"

# Create a set of objects to read from
process.selectedObjects = cms.EDProducer("PFCleaner",
    vertices = cms.InputTag("goodVertices"),
    pfcands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag(jetCollName),
    electronidveto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronidmedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    photonidloose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose")
)

# Define all the METs corrected for lepton/photon momenta
process.partMet = cms.EDProducer("METBreakDownProducer",
    pfcands = cms.InputTag("packedPFCandidates") 
)

process.noHFCands = cms.EDFilter("CandPtrSelector",
    src=cms.InputTag("packedPFCandidates"),
    cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
)

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
process.t1mumet = cms.EDProducer("MuonCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    muons = cms.InputTag("selectedObjects", "muons")
)

process.elmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
    useuncorrmet = cms.bool(True)
)
process.t1elmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
)
process.phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
    useuncorrmet = cms.bool(True)
)
process.t1phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
)

# Quark-Gluon Discriminant
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = jetCollName
process.QGTagger.srcVertexCollection = "goodVertices"

# Make the tree 
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
    photonMediumId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
    photonTightId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag(jetCollName),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    qgl = cms.InputTag("QGTagger", "qgLikelihood"),
    qgs2 = cms.InputTag("QGTagger", "axis2"),
    qgmult = cms.InputTag("QGTagger", "mult"),
    qgptd = cms.InputTag("QGTagger", "ptD"),
    partmet = cms.InputTag("partMet"),
    t1pfmet = cms.InputTag("slimmedMETs"),
    pfmupt = cms.InputTag("pfmupt"),
    mumet = cms.InputTag("mumet"),
    t1mumet = cms.InputTag("t1mumet"),
    elmet = cms.InputTag("elmet"),
    t1elmet = cms.InputTag("t1elmet"),
    phmet = cms.InputTag("phmet"),
    t1phmet = cms.InputTag("t1phmet"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    filterResults = cms.InputTag("TriggerResults", "", miniAODProcess),
    hcalnoise = cms.InputTag("hcalnoise"),
    hbheloose = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
    hbhetight = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    xsec = cms.double(0.001),
    cleanMuonJet = cms.bool(True),
    cleanElectronJet = cms.bool(True),
    cleanPhotonJet = cms.bool(True),
    applyHLTFilter = cms.bool(filterOnHLT),
    uselheweights = cms.bool(False),
    isWorZMCSample = cms.bool(False)
)

# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheinfo = cms.InputTag("externalLHEProducer"),
    geninfo = cms.InputTag("generator"),
    uselheweights = cms.bool(False),
    addqcdpdfweights = cms.bool(False)
)

# MET filter
process.metfilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("t1mumet"),
    cut = cms.string("et > 200"),
    filter = cms.bool(True)
)

# Set up the path
if filterHighMETEvents: 
    if (isMC):
        process.treePath = cms.Path(process.gentree + process.goodVertices + process.metfilter + process.tree)
    else :
        process.treePath = cms.Path(                  process.goodVertices + process.metfilter + process.tree)
else :
    if (isMC):
        process.treePath = cms.Path(process.gentree + process.goodVertices                     + process.tree)
    else :
        process.treePath = cms.Path(                  process.goodVertices                     + process.tree)
