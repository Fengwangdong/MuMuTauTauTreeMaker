from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'YOURS'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../../runDiMuDiTau_cfg.py'
config.JobType.pyCfgParams = ["isMC=1","numThreads=4","era=YOURS"]
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = 'YOURS'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/YOURS/Treelization/'
config.Data.publication = False
config.Data.outputDatasetTag = 'YOURS'

config.Site.storageSite = 'YOURS'
