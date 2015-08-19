from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MonoJet_dmAVM10'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'tree.py'
config.JobType.inputFiles = ['Summer15_50nsV2_MC.db']

config.section_("Data")
config.Data.inputDataset = '/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'Monojet_dmAVM10'

config.section_("Site")
config.Site.storageSite = 'T2_US_UCSD'

