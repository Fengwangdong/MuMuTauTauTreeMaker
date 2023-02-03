import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# -------- input files. Can be changed on the command line with the option inputFiles=... ---------
options.inputFiles = ['/store/group/phys_higgs/HiggsExo/fengwang/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-19_TuneCUETP8M1_13TeV_madgraph_pythia8/MiniAOD_H125AA19_DiMuDiTau_Fall17DRPremix_v1/190515_140053/0000/mumutautau_1.root']
options.register('isMC', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Sample is MC")
options.register('numThreads', 4, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Set number of CPU cores")
options.register('era', '2017', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Year of sample")
options.parseArguments()

process = cms.Process("JetDiTauTreelizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff') # need for RecoEgamma recipe
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') # globaltag is also needed for RecoEgamma recipe inclusion
from Configuration.AlCa.GlobalTag import GlobalTag

########## Please specify if you are running on data (0) or MC (1) in the command line: #########################
########### eg: cmsRun runJetDiTau_cfg.py isMC=1 ###############
##########################################################################
if options.isMC == 1:
    print (" ****** we will run on sample of: MC ******")
    if options.era == '2016preVFP':
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_preVFP_v11'

    elif options.era == '2016postVFP':
        process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v17'

    elif options.era == '2017':
        process.GlobalTag.globaltag = '106X_mc2017_realistic_v9'

    else:
        process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v16_L1v1'

    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.JetDiTauSelectorMC_cfi")

else:
    print (" ****** we will run on sample of: data ******")
    process.GlobalTag.globaltag = '106X_dataRun2_v35'
    process.load("MuMuTauTauTreeMaker.MuTauTreelizer.JetDiTauSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# --- please specify the sample that you need to run for local test ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

######## embed the egamma ID into the nanoAOD ############
# reference: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
from MuMuTauTauTreeMaker.MuTauTreelizer.EgammaPostRecoTools import setupEgammaPostRecoSeq

if options.era == '2016preVFP':
    setupEgammaPostRecoSeq(process, era='2016preVFP-UL')

elif options.era == '2016postVFP':
    setupEgammaPostRecoSeq(process, era='2016postVFP-UL')

elif options.era == '2017':
    setupEgammaPostRecoSeq(process, era='2017-UL')

else:
    setupEgammaPostRecoSeq(process, era='2018-UL')
###########################################################

if options.isMC == 1:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonCandSelector*
            process.egammaPostRecoSeq*
            process.ElectronCandSelector*
            process.TrigJetMatcher*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.GenMuonCandSelector*
            process.GenElectronCandSelector*
            process.GenTauMuCandSelector*
            process.GenTauEleCandSelector*
            process.GenTauHadCandSelector*
            process.JetDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_mc.root')
    )

else:
    process.treelizer = cms.Sequence(
            process.lumiTree*
            process.HLTEle*
            process.MuonCandSelector*
            process.egammaPostRecoSeq*
            process.ElectronCandSelector*
            process.TrigJetMatcher*
            process.DeepDiTauProducer*
            process.JetIdEmbedder*
            process.JetDiTauAnalyzer
    )

    process.TFileService = cms.Service("TFileService",
            fileName =  cms.string('MuMuTauTauTreelization_data.root')
    )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
)

process.options.numberOfThreads = cms.untracked.uint32(options.numThreads)
process.p = cms.Path(process.treelizer)