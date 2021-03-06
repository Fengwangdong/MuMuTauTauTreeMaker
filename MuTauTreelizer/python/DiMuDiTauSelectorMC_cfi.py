import FWCore.ParameterSet.Config as cms

lumiTree = cms.EDAnalyzer("LumiTree",
        genEventInfo = cms.InputTag("generator"),
        nevents = cms.InputTag('lumiSummary','numberOfEvents'),
        summedWeights = cms.InputTag('lumiSummary','sumOfWeightedEvents'),
)

HLTEle = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*","HLT_IsoMu27_v*","HLT_IsoTkMu27_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False), # throw exception on unknown path names
)

HLTEleDiMu = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
        HLTPaths = cms.vstring("HLT_Dimuon20_Jpsi_v*", "HLT_Dimuon25_Jpsi_v*","HLT_DoubleMu4_JpsiTrk_Displaced_v*"),
        eventSetupPathsKey = cms.string(''),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False), # throw exception on unknown path names
)

MuonID = cms.EDFilter("MuonID",
        muonTag = cms.InputTag('slimmedMuons'),
        muonID = cms.string('loose'),
        minNumObjsToPassFilter = cms.int32(2),
)

MuonSelector = cms.EDFilter("MuonSelector",
        muonTag = cms.InputTag("MuonID"),
        relIsoCutVal = cms.double(-1), # positive number for iso threshold, -1 for ignoring iso
        normalRelIso = cms.bool(True), #True = Iso-mu; False = inverted Iso-mu
        Eta = cms.double(2.5),
        Pt = cms.double(3.0),
)

TrigMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('MuonSelector'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v","HLT_IsoTkMu27_v"),
        dRCut = cms.double(0.15),
        numberOfTrigMus = cms.int32(1),
        muPtCut = cms.double(26.0),
)

TrigDiMuMatcher = cms.EDFilter("TrigMuMatcher",
        muonsTag = cms.InputTag('MuonSelector'),
        bits = cms.InputTag("TriggerResults","","HLT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        trigNames = cms.vstring("HLT_Dimuon20_Jpsi_v", "HLT_Dimuon25_Jpsi_v","HLT_DoubleMu4_JpsiTrk_Displaced_v"),
        dRCut = cms.double(0.15),
        numberOfTrigMus = cms.int32(2),
        muPtCut = cms.double(3.0),
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

TauCandSelector = cms.EDFilter("TauCandSelector",
        tauTag = cms.InputTag('slimmedTausNewID'),
        tauDiscriminatorTag = cms.vstring('decayModeFindingNewDMs'),
        passDiscriminator = cms.bool(False),
        pTMin = cms.double(8.0),
        etaMax = cms.double(2.4),
)

DeepDiTauProducer = cms.EDProducer("DeepDiTauProducer",
        slimmedJetTag = cms.InputTag('slimmedJets'),
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
            ),
        ),
)

JetIdEmbedder = cms.EDProducer("JetIdEmbedder",
        slimmedJetTag = cms.InputTag('slimmedJets'),
        discriminator = cms.string('pileupJetId:fullDiscriminant'),
        ditau2017v1 = cms.InputTag("DeepDiTauProducer","ditau2017v1"),
        ditau2017MDv1 = cms.InputTag("DeepDiTauProducer","ditau2017MDv1"),
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

DiMuDiTauAnalyzer = cms.EDAnalyzer('DiMuDiTauAnalyzer',
        MuTag = cms.InputTag("TrigMuMatcher"),
        EleTag = cms.InputTag("ElectronCandSelector"),
        TauTag = cms.InputTag("TauCandSelector"),
        JetTag = cms.InputTag("JetIdEmbedder"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        isMC = cms.bool(True),
        numberOfTrigMus = cms.int32(1),
        GenMuTag = cms.InputTag('GenMuonCandSelector'),
        GenEleTag = cms.InputTag('GenElectronCandSelector'),
        GenTauMuTag = cms.InputTag('GenTauMuCandSelector'),
        GenTauEleTag = cms.InputTag('GenTauEleCandSelector'),
        GenTauHadTag = cms.InputTag('GenTauHadCandSelector'),
        PileupTag = cms.InputTag("slimmedAddPileupInfo"),
        Generator = cms.InputTag("generator"),
)

DiMuDiTauAnalyzerDiMuTrig = cms.EDAnalyzer('DiMuDiTauAnalyzer',
        MuTag = cms.InputTag("TrigDiMuMatcher"),
        EleTag = cms.InputTag("ElectronCandSelector"),
        TauTag = cms.InputTag("TauCandSelector"),
        JetTag = cms.InputTag("JetIdEmbedder"),
        MetTag = cms.InputTag("slimmedMETs"),
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        rhoTag = cms.InputTag("fixedGridRhoAll"),
        effAreasConfigFile = cms.FileInPath("MuMuTauTauTreeMaker/MuTauTreelizer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
        isMC = cms.bool(True),
        numberOfTrigMus = cms.int32(2),
        GenMuTag = cms.InputTag('GenMuonCandSelector'),
        GenEleTag = cms.InputTag('GenElectronCandSelector'),
        GenTauMuTag = cms.InputTag('GenTauMuCandSelector'),
        GenTauEleTag = cms.InputTag('GenTauEleCandSelector'),
        GenTauHadTag = cms.InputTag('GenTauHadCandSelector'),
        PileupTag = cms.InputTag("slimmedAddPileupInfo"),
        Generator = cms.InputTag("generator"),
)
