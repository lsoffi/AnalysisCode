from WMCore.Configuration import Configuration
config = Configuration()

name = 'EffAN_7413Update_Spring15'
proc = 'WLNu_HT200to400_V11'
dataset = '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'

# GENERAL
config.section_("General")
config.General.requestName = name+"_"+proc
config.General.workArea    = '/user/ndaci/CRABBY/StudyHLT/Monojet/'+name+'/'+proc+'/'
config.General.transferLogs    = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nadirMC.py'
config.JobType.pyCfgParams = ['theXSec=359600']
##config.JobType.inputFiles = '../../data/dttf_config.db'
#config.JobType.inputFiles = ['Summer15_50nsV4_MC.db']
config.JobType.allowUndistributedCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits  = 1976
config.Data.publication = False
config.Data.publishDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = name+"_"+proc
config.Data.ignoreLocality = False # allows to process inputs on CE != site hosting inputs
#config.Data.lumiMask = 
#config.Data.runRange = 

#A custom string to insert in the output file name inside the CRAB-created directory path to allow organizing groups of tasks.
#config.Data.prefix =  

# USER
config.section_("User")
#config.User.email = 'nadir.daci@cern.ch'
#config.User.voRole = 
config.User.voGroup = 'becms'

# GRID
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.whitelist = ['T2_BE_IIHE']
#config.Site.blacklist = ['T1_US_FNAL']
