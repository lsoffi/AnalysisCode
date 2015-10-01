from CRABClient.UserUtilities import config
config = config()

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'tree.py'

config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000

config.Site.storageSite = 'T2_US_UCSD'

