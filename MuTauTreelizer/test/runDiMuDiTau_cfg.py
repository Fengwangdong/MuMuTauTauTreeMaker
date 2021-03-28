import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('numThreads', 8, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Set number of CPU cores")
options.parseArguments()

process = cms.Process("DiMuonDiTauTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runDiMuDiTau_cfg.py isMC=1 ###############
##########################################################################
if options.isMC == 1:
    print " ****** we will run on sample of: MC ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuDiTauSelectorMC_cfi")

else:
    print " ****** we will run on sample of: data ******"
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.DiMuDiTauSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

######### embed 2017v2 tauID into the miniAOD ###############
# reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Rerunning_of_the_tau_ID_on_M_AN1
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

print " ====== use slimmedTausMuonCleaned cluster ======"
updatedTauNameMuonCleaned = "slimmedTausMuonCleanedNewID"
import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausMuonCleaned as tauIdConfigMuonCleaned
tauIdEmbedderMuonCleaned = tauIdConfigMuonCleaned.TauIDEmbedderMuonCleaned(process, cms,
        debug = True,
        updatedTauName = updatedTauNameMuonCleaned,
        toKeep = ["deepTau2017v2p1MuonCleaned","2017v2MuonCleaned"]
        )
tauIdEmbedderMuonCleaned.runTauID()
process.rerunTauMuonCleanedIDSequence = cms.Sequence(process.rerunMvaIsolationSequenceMuonCleaned * getattr(process,updatedTauNameMuonCleaned))

print " ====== use slimmedTausElectronCleaned cluster ======"
updatedTauNameElectronCleaned = "slimmedTausElectronCleanedNewID"
import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausElectronCleaned as tauIdConfigElectronCleaned
tauIdEmbedderElectronCleaned = tauIdConfigElectronCleaned.TauIDEmbedderElectronCleaned(process, cms,
        debug = True,
        updatedTauName = updatedTauNameElectronCleaned,
        toKeep = ["deepTau2017v2p1ElectronCleaned","2017v2ElectronCleaned"]
        )
tauIdEmbedderElectronCleaned.runTauID()
process.rerunTauElectronCleanedIDSequence = cms.Sequence(process.rerunMvaIsolationSequenceElectronCleaned * getattr(process,updatedTauNameElectronCleaned))

print " ====== use slimmedTausBoosted cluster ======"
updatedTauNameBoosted = "slimmedTausBoostedNewID"
import MuMuTauTauTreeMaker.MuTauTreelizer.TauIdDeep_slimmedTausBoosted as tauIdConfigBoosted
tauIdEmbedderBoosted = tauIdConfigBoosted.TauIDEmbedderBoosted(process, cms,
        debug = False,
        updatedTauName = updatedTauNameBoosted,
        PATTauProducer = cms.InputTag('cleanedSlimmedTausBoosted'),
        srcChargedIsoPtSum = cms.string('chargedIsoPtSumNoOverLap'),
        srcNeutralIsoPtSum = cms.string('neutralIsoPtSumNoOverLap'),
        toKeep = ["deepTau2017v2p1Boosted","2017v2Boosted"]
        )
tauIdEmbedderBoosted.runTauID()
process.rerunTauBoostedIDSequence = cms.Sequence(process.ca8PFJetsCHSprunedForBoostedTausPAT * getattr(process, "cleanedSlimmedTausBoosted") * process.rerunMvaIsolationSequenceBoosted * getattr(process,updatedTauNameBoosted))
############################################################

if options.isMC == 1:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.rerunTauMuonCleanedIDSequence*
            process.TauMuonCleanedCandSelector*
            process.rerunTauElectronCleanedIDSequence*
            process.TauElectronCleanedCandSelector*
            process.rerunTauBoostedIDSequence*
            process.TauBoostedCandSelector*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.GenMuonCandSelector*
            process.GenElectronCandSelector*
            process.GenTauMuCandSelector*
            process.GenTauEleCandSelector*
            process.GenTauHadCandSelector*
            process.DiMuDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_mc.root')
    )

else:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonID*
            process.MuonSelector*
            process.TrigMuMatcher*
            process.ElectronCandSelector*
            process.rerunTauIDSequence*
            process.TauCandSelector*
            process.rerunTauMuonCleanedIDSequence*
            process.TauMuonCleanedCandSelector*
            process.rerunTauElectronCleanedIDSequence*
            process.TauElectronCleanedCandSelector*
            process.rerunTauBoostedIDSequence*
            process.TauBoostedCandSelector*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.DiMuDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_data.root')
    )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.options.numberOfThreads = cms.untracked.uint32(options.numThreads)
process.p = cms.Path(process.treelizer)
