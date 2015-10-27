from WMCore.Configuration import Configuration
config = Configuration()

name = 'EffAN_7413Update_Run2015D_PRV4'
proc = 'DoubleEG_V2'
dataset = '/DoubleEG/Run2015D-PromptReco-v4/MINIAOD'

# GENERAL
config.section_("General")
config.General.requestName = name+"_"+proc
config.General.workArea    = '/user/ndaci/CRABBY/StudyHLT/Monojet/'+name+'/'+proc+'/'
config.General.transferLogs    = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nadirDataPrompt.py'
#config.JobType.pyCfgParams = ['reEmulation=True','reEmulMuons=True','reEmulCalos=True','patchNtuple=True','force2012Config=True','customDTTF=True','dttfLutsFile=sqlite:src/L1TriggerDPG/L1Menu/data/dttf_config.db','useUct2015=True','globalTag=POSTLS162_V2::All','runOnMC=True','runOnPostLS1=True','whichPU=40']
##config.JobType.inputFiles = '../../data/dttf_config.db'
config.JobType.inputFiles = ['Summer15_25nsV5_DATA.db']
config.JobType.allowUndistributedCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
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
#config.Site.whitelist = 
config.Site.blacklist = ['T1_US_FNAL','T2_UA_KIPT','T3_US_UCR']
