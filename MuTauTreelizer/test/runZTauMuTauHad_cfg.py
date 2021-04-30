import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('tauCluster', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "different tau clusters")
options.register('numThreads', 4, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Set number of CPU cores")
options.parseArguments()

process = cms.Process("ZMuonTauTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runZTauMuTauHad_cfg.py isMC=1 ###############
##########################################################################
if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.ZTauMuTauHadSelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.ZTauMuTauHadSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

######### embed 2017v2 tauID into the miniAOD ###############
# reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_M_AN1
if options.tauCluster <= 0:
    print " ====== use slimmedTaus cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTaus as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 1:
    print " ====== use slimmedTausMuonCleaned cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausMuonCleaned as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))

elif options.tauCluster == 2:
    print " ====== use slimmedTausElectronCleaned cluster ======"
    updatedTauName = "slimmedTausNewID"
    import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausElectronCleaned as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms,
            debug = True,
            updatedTauName = updatedTauName,
            toKeep = ["deepTau2017v2p1","2017v2"]
            )
    tauIdEmbedder.runTauID()
    process.rerunTauIDSequence = cms.Sequence(process.rerunMvaIsolationSequence * getattr(process,updatedTauName))
############################################################

if options.isMC == 1:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.rerunTauIDSequence*
            process.TauHadSelector*
            process.JetSelector*
            process.GenMuonCandSelector*
            process.GenTauMuCandSelector*
            process.GenTauHadCandSelector*
            process.ZTauMuTauHadAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('ZTauMuTauHadTreelization_mc.root')
    )

else:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.rerunTauIDSequence*
            process.TauHadSelector*
            process.JetSelector*
            process.ZTauMuTauHadAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('ZTauMuTauHadTreelization_data.root')
    )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.options.numberOfThreads = cms.untracked.uint32(options.numThreads)
process.p = cms.Path(process.treelizer)
