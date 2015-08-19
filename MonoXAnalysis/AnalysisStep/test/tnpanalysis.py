import FWCore.ParameterSet.Config as cms

process = cms.Process("TNP")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


isMC = False

InputFileName = "tnptree.root"
if isMC : 
    OutputFilePrefix = "efficiency-mc"
else :
    OutputFilePrefix = "efficiency-data"
PDFName = "pdfSignalPlusBackground"

EfficiencyBins = cms.PSet(
    pt = cms.vdouble( 10., 20., 30., 40., 50., 70., 100.),
    abseta = cms.vdouble( 1.2, 2.4 )
)
EfficiencyBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
        EfficiencyBins,
    ),
    BinToPDFmap = cms.vstring([
        "pdfSignalPlusBackground", 
        "*pt_bin0*", "pdfSignalPlusBackgroundb0",
        "*pt_bin1*", "pdfSignalPlusBackgroundb1",
        "*pt_bin2*", "pdfSignalPlusBackgroundb2",
        "*pt_bin3*", "pdfSignalPlusBackgroundb3",
        "*pt_bin4*", "pdfSignalPlusBackgroundb4",
        "*pt_bin5*", "pdfSignalPlusBackgroundb5",
    ])  
)

IDModule = cms.PSet(
    Id = cms.PSet(
        EfficiencyBinningSpecification,
        EfficiencyCategoryAndState = cms.vstring("tightid","pass"),
    )
)

process.muonIdTnP = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    InputFileNames = cms.vstring(InputFileName),
    InputDirectoryName = cms.string("muontnptree"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(OutputFilePrefix+".root"),
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    binsForMassPlots = cms.uint32(25),
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "65.0", "115.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("Probe #eta", "-2.4", "2.4", ""),                
    ),

    Categories = cms.PSet(
        mcTrue = cms.vstring("MC Truth", "dummy[true=1,false=0]"),
        tightid = cms.vstring("Tight ID", "dummy[pass=1,fail=0]"), 
    ),
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
            "BreitWigner::signalPass(mass, mP[91, 89, 93], sP[3.0, 0.0, 7.0])",
            "BreitWigner::signalFail(mass, mF[91, 89, 93], sF[3.0, 0.0, 7.0])",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb0 = cms.vstring(
            "#import templates0.root:w:signalPassMC",
            "#import templates0.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            #"Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
            #"Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb1 = cms.vstring(
            "#import templates1.root:w:signalPassMC",
            "#import templates1.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb2 = cms.vstring(
            "#import templates2.root:w:signalPassMC",
            "#import templates2.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb3 = cms.vstring(
            "#import templates3.root:w:signalPassMC",
            "#import templates3.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb4 = cms.vstring(
            "#import templates4.root:w:signalPassMC",
            "#import templates4.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
        pdfSignalPlusBackgroundb5 = cms.vstring(
            "#import templates5.root:w:signalPassMC",
            "#import templates5.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
            "efficiency[0.8,0,1]"
        ),
    ),

    Efficiencies = cms.PSet(
        IDModule
    )
)

process.fit = cms.Path(process.muonIdTnP)
