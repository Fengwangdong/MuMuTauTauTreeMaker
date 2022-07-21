import FWCore.ParameterSet.Config as cms

lumiTree = cms.EDAnalyzer("LumiTree",
        genEventInfo = cms.InputTag("generator"),
        nevents = cms.InputTag('lumiSummary','numberOfEvents'),
        summedWeights = cms.InputTag('lumiSummary','sumOfWeightedEvents'),
)

HLTEle = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_PFJet500_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False), # throw exception on unknown path names
)

MuonCandSelector = cms.EDFilter("MuonCandSelector",
        muonTag = cms.InputTag('slimmedMuons'),
        muonID = cms.string('loose'),
        relIsoCutVal = cms.double(-1), # positive number for iso threshold, -1 for ignoring iso
        normalRelIso = cms.bool(True), #True = Iso-mu; False = inverted Iso-mu
        Eta = cms.double(2.5),
        Pt = cms.double(3.0),
)

ElectronCandSelector = cms.EDFilter("ElectronCandSelector",
        electronTag = cms.InputTag('slimmedElectrons'),
        # --- need the two parameters below for electron isolation computation ---
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        # ========================================================================
        relIdName = cms.string("veto"), # customize electron ID options: Veto, Loose, Medium, Tight
        passRelIso = cms.bool(False),
        etaCut = cms.double(2.5),
        ptCut = cms.double(3.0),
)

TrigJetMatcher = cms.EDFilter("TrigJetMatcher",
        jetsTag = cms.InputTag('slimmedJets'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_PFJet500_v"),
        dRCut = cms.double(0.15),
        jetPtCut = cms.double(502.0),
)


DeepDiTauProducer = cms.EDProducer("DeepDiTauProducer",
        slimmedJetTag = cms.InputTag('TrigJetMatcher'),
        DeepDiTauConfiguration = cms.PSet(
            memmapped = cms.bool(False),
            graphDefinitions = cms.VPSet(
                cms.PSet(
                    name = cms.string('ditau2017v1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_v1.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_v1_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('ditau2017MDv1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_md_v1.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_md_v1_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('ditau2017v2'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_v2.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_v2_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('ditau2017MDv2'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_md_v2.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/ditau_2017_md_v2_means_sigmas.txt'),
            ),
        ),
    ),
)

JetIdEmbedder = cms.EDProducer("JetIdEmbedder",
        slimmedJetTag = cms.InputTag('TrigJetMatcher'),
        discriminator = cms.string('pileupJetId:fullDiscriminant'),
        ditau2017v1 = cms.InputTag("DeepDiTauProducer","ditau2017v1"),
        ditau2017MDv1 = cms.InputTag("DeepDiTauProducer","ditau2017MDv1"),
        ditau2017v2 = cms.InputTag("DeepDiTauProducer","ditau2017v2"),
        ditau2017MDv2 = cms.InputTag("DeepDiTauProducer","ditau2017MDv2"),
)

GenMuonCandSelector = cms.EDFilter("GenMuonCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenElectronCandSelector = cms.EDFilter("GenElectronCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauMuCandSelector = cms.EDFilter("GenTauMuCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauEleCandSelector = cms.EDFilter("GenTauEleCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

GenTauHadCandSelector = cms.EDFilter("GenTauHadCandSelector",
        genParticlesTag = cms.InputTag('prunedGenParticles'),
        etaCut = cms.double(2.6),
        ptCut = cms.double(2.5),
)

JetDiTauAnalyzer = cms.EDAnalyzer('JetDiTauAnalyzer',
        MuTag = cms.InputTag("MuonCandSelector"),
        EleTag = cms.InputTag("ElectronCandSelector"),
        JetTag = cms.InputTag("JetIdEmbedder"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        isMC = cms.bool(True),
        GenMuTag = cms.InputTag('GenMuonCandSelector'),
        GenEleTag = cms.InputTag('GenElectronCandSelector'),
        GenTauMuTag = cms.InputTag('GenTauMuCandSelector'),
        GenTauEleTag = cms.InputTag('GenTauEleCandSelector'),
        GenTauHadTag = cms.InputTag('GenTauHadCandSelector'),
        PileupTag = cms.InputTag("slimmedAddPileupInfo"),
        Generator = cms.InputTag("generator"),
)
