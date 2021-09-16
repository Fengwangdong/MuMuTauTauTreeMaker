// -*- C++ -*-
//
// Package:    MuMuChannel/DiMuDiTauAnalyzer
// Class:      DiMuDiTauAnalyzer
// 
/**\class DiMuDiTauAnalyzer DiMuDiTauAnalyzer.cc MuMuChannel/DiMuDiTauAnalyzer/plugins/DiMuDiTauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Tue, 09 Apr 2019 16:30:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <math.h>
#include <string>
#include <iostream>

using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DiMuDiTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiMuDiTauAnalyzer(const edm::ParameterSet&);
      ~DiMuDiTauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      std::vector<const reco::Candidate*> findTauMuEleVisDaughters(const reco::Candidate*);
      std::vector<const reco::Candidate*> findTauHadVisDaughters(const reco::Candidate*);
      int findTauPiZeros(const reco::Candidate*);
      int findTauChargedHadrons(const reco::Candidate*);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> MuTag;
      edm::EDGetTokenT<edm::View<pat::Electron>> EleTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauMuonCleanedTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauElectronCleanedTag;
      edm::EDGetTokenT<edm::View<pat::Tau>> TauBoostedTag;
      edm::EDGetTokenT<edm::View<pat::Jet>> JetTag;
      edm::EDGetTokenT<edm::View<pat::MET>> MetTag;
      edm::EDGetTokenT<edm::View<reco::Vertex>> VertexTag;
      edm::EDGetTokenT<double> rhoTag;
      EffectiveAreas effectiveAreas;
      bool isMC;
      int numberOfTrigMus;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo>> PileupTag;
      edm::EDGetTokenT<GenEventInfoProduct> generator;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenMuTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenEleTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauMuTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauEleTag;
      edm::EDGetTokenT<edm::View<reco::GenParticle>> GenTauHadTag;

      TTree *objectTree;
      // --- below is the vectors of object variables ---
      
      // --- reconstructed muons ---
      vector<float> recoMuonPt;
      vector<float> recoMuonEta;
      vector<float> recoMuonPhi;
      vector<float> recoMuonEnergy;
      vector<int> recoMuonPDGId;
      vector<float> recoMuonIsolation;
      vector<float> recoMuonDXY;
      vector<float> recoMuonDZ;
      vector<int> recoMuonNTrackerLayers;
      vector<int> recoMuonTriggerFlag;
      vector<int> recoMuonRefToElectron;
      vector<int> recoMuonRefToTau;
      vector<int> recoMuonRefToTauMuonCleaned;
      vector<int> recoMuonIdLoose;
      vector<int> recoMuonIdMedium;
      vector<int> recoMuonIdTight;

      // --- reconstructed electrons ---
      vector<float> recoElectronPt;
      vector<float> recoElectronEta;
      vector<float> recoElectronPhi;
      vector<float> recoElectronEnergy;
      vector<int> recoElectronPDGId;
      vector<float> recoElectronIsolation;
      vector<int> recoElectronIdLoose;
      vector<int> recoElectronIdMedium;
      vector<int> recoElectronIdTight;
      vector<int> recoElectronIdLooseNoIso;
      vector<int> recoElectronIdMediumNoIso;
      vector<int> recoElectronIdTightNoIso;
      vector<float> recoElectronEcalTrkEnergyPostCorr;
      vector<float> recoElectronEcalTrkEnergyErrPostCorr;
      vector<float> recoElectronEnergyScaleValue;
      vector<float> recoElectronEnergyScaleUp;
      vector<float> recoElectronEnergyScaleDown;
      vector<float> recoElectronEnergySigmaValue;
      vector<float> recoElectronEnergySigmaUp;
      vector<float> recoElectronEnergySigmaDown;
      vector<int> recoElectronRefToMuon;
      vector<int> recoElectronRefToTau;
      vector<int> recoElectronRefToTauElectronCleaned;

      // --- reconstructed taus ---
      vector<float> recoTauPt;
      vector<float> recoTauEta;
      vector<float> recoTauPhi;
      vector<float> recoTauEnergy;
      vector<int> recoTauPDGId;
      vector<float> recoTauDecayMode;
      vector<float> recoTauDecayModeFinding;
      vector<float> recoTauDecayModeFindingNewDMs;
      vector<float> recoTauIsoMVArawValue;
      vector<float> recoTauIsoMVAVVLoose;
      vector<float> recoTauIsoMVAVLoose;
      vector<float> recoTauIsoMVALoose;
      vector<float> recoTauIsoMVAMedium;
      vector<float> recoTauIsoMVATight;
      vector<float> recoTauIsoMVAVTight;
      vector<float> recoTauIsoMVAVVTight;
      
      vector<float> recoTauAntiMuMVALoose;
      vector<float> recoTauAntiMuMVATight;

      vector<float> recoTauAntiEleMVArawValue;
      vector<float> recoTauAntiEleMVAVLoose;
      vector<float> recoTauAntiEleMVALoose;
      vector<float> recoTauAntiEleMVAMedium;
      vector<float> recoTauAntiEleMVATight;
      vector<float> recoTauAntiEleMVAVTight;

      // *** deep tau ID variables ***
      vector<float> recoTauDeepVSeraw;
      vector<float> recoTauDeepVSjetraw;
      vector<float> recoTauDeepVSmuraw;
 
      vector<float> recoTauDeepVSeLoose;
      vector<float> recoTauDeepVSjetLoose;
      vector<float> recoTauDeepVSmuLoose;

      vector<float> recoTauDeepVSeMedium;
      vector<float> recoTauDeepVSjetMedium;
      vector<float> recoTauDeepVSmuMedium;

      vector<float> recoTauDeepVSeTight;
      vector<float> recoTauDeepVSjetTight;
      vector<float> recoTauDeepVSmuTight;

      vector<float> recoTauDeepVSeVLoose;
      vector<float> recoTauDeepVSjetVLoose;
      vector<float> recoTauDeepVSmuVLoose;

      vector<float> recoTauDeepVSeVTight;
      vector<float> recoTauDeepVSjetVTight;

      vector<float> recoTauDeepVSeVVLoose;
      vector<float> recoTauDeepVSjetVVLoose;

      vector<float> recoTauDeepVSeVVTight;
      vector<float> recoTauDeepVSjetVVTight;

      vector<float> recoTauDeepVSeVVVLoose;
      vector<float> recoTauDeepVSjetVVVLoose;

      vector<int> recoTauRefToMuon;
      vector<int> recoTauRefToElectron;

      // --- reconstructed taus muonCleaned ---
      vector<float> recoTauMuonCleanedPt;
      vector<float> recoTauMuonCleanedEta;
      vector<float> recoTauMuonCleanedPhi;
      vector<float> recoTauMuonCleanedEnergy;
      vector<int> recoTauMuonCleanedPDGId;
      vector<float> recoTauMuonCleanedDecayMode;
      vector<float> recoTauMuonCleanedDecayModeFinding;
      vector<float> recoTauMuonCleanedDecayModeFindingNewDMs;
      vector<float> recoTauMuonCleanedIsoMVArawValue;
      vector<float> recoTauMuonCleanedIsoMVAVVLoose;
      vector<float> recoTauMuonCleanedIsoMVAVLoose;
      vector<float> recoTauMuonCleanedIsoMVALoose;
      vector<float> recoTauMuonCleanedIsoMVAMedium;
      vector<float> recoTauMuonCleanedIsoMVATight;
      vector<float> recoTauMuonCleanedIsoMVAVTight;
      vector<float> recoTauMuonCleanedIsoMVAVVTight;
      
      vector<float> recoTauMuonCleanedAntiMuMVALoose;
      vector<float> recoTauMuonCleanedAntiMuMVATight;

      vector<float> recoTauMuonCleanedAntiEleMVArawValue;
      vector<float> recoTauMuonCleanedAntiEleMVAVLoose;
      vector<float> recoTauMuonCleanedAntiEleMVALoose;
      vector<float> recoTauMuonCleanedAntiEleMVAMedium;
      vector<float> recoTauMuonCleanedAntiEleMVATight;
      vector<float> recoTauMuonCleanedAntiEleMVAVTight;

      // *** deep tau ID variables ***
      vector<float> recoTauMuonCleanedDeepVSeraw;
      vector<float> recoTauMuonCleanedDeepVSjetraw;
      vector<float> recoTauMuonCleanedDeepVSmuraw;
 
      vector<float> recoTauMuonCleanedDeepVSeLoose;
      vector<float> recoTauMuonCleanedDeepVSjetLoose;
      vector<float> recoTauMuonCleanedDeepVSmuLoose;

      vector<float> recoTauMuonCleanedDeepVSeMedium;
      vector<float> recoTauMuonCleanedDeepVSjetMedium;
      vector<float> recoTauMuonCleanedDeepVSmuMedium;

      vector<float> recoTauMuonCleanedDeepVSeTight;
      vector<float> recoTauMuonCleanedDeepVSjetTight;
      vector<float> recoTauMuonCleanedDeepVSmuTight;

      vector<float> recoTauMuonCleanedDeepVSeVLoose;
      vector<float> recoTauMuonCleanedDeepVSjetVLoose;
      vector<float> recoTauMuonCleanedDeepVSmuVLoose;

      vector<float> recoTauMuonCleanedDeepVSeVTight;
      vector<float> recoTauMuonCleanedDeepVSjetVTight;

      vector<float> recoTauMuonCleanedDeepVSeVVLoose;
      vector<float> recoTauMuonCleanedDeepVSjetVVLoose;

      vector<float> recoTauMuonCleanedDeepVSeVVTight;
      vector<float> recoTauMuonCleanedDeepVSjetVVTight;

      vector<float> recoTauMuonCleanedDeepVSeVVVLoose;
      vector<float> recoTauMuonCleanedDeepVSjetVVVLoose;

      vector<int> recoTauMuonCleanedRefToMuon;
      vector<int> recoTauMuonCleanedRefToElectron;

      // --- reconstructed taus electronCleaned ---
      vector<float> recoTauElectronCleanedPt;
      vector<float> recoTauElectronCleanedEta;
      vector<float> recoTauElectronCleanedPhi;
      vector<float> recoTauElectronCleanedEnergy;
      vector<int> recoTauElectronCleanedPDGId;
      vector<float> recoTauElectronCleanedDecayMode;
      vector<float> recoTauElectronCleanedDecayModeFinding;
      vector<float> recoTauElectronCleanedDecayModeFindingNewDMs;
      vector<float> recoTauElectronCleanedIsoMVArawValue;
      vector<float> recoTauElectronCleanedIsoMVAVVLoose;
      vector<float> recoTauElectronCleanedIsoMVAVLoose;
      vector<float> recoTauElectronCleanedIsoMVALoose;
      vector<float> recoTauElectronCleanedIsoMVAMedium;
      vector<float> recoTauElectronCleanedIsoMVATight;
      vector<float> recoTauElectronCleanedIsoMVAVTight;
      vector<float> recoTauElectronCleanedIsoMVAVVTight;
      
      vector<float> recoTauElectronCleanedAntiMuMVALoose;
      vector<float> recoTauElectronCleanedAntiMuMVATight;

      vector<float> recoTauElectronCleanedAntiEleMVArawValue;
      vector<float> recoTauElectronCleanedAntiEleMVAVLoose;
      vector<float> recoTauElectronCleanedAntiEleMVALoose;
      vector<float> recoTauElectronCleanedAntiEleMVAMedium;
      vector<float> recoTauElectronCleanedAntiEleMVATight;
      vector<float> recoTauElectronCleanedAntiEleMVAVTight;

      // *** deep tau ID variables ***
      vector<float> recoTauElectronCleanedDeepVSeraw;
      vector<float> recoTauElectronCleanedDeepVSjetraw;
      vector<float> recoTauElectronCleanedDeepVSmuraw;
 
      vector<float> recoTauElectronCleanedDeepVSeLoose;
      vector<float> recoTauElectronCleanedDeepVSjetLoose;
      vector<float> recoTauElectronCleanedDeepVSmuLoose;

      vector<float> recoTauElectronCleanedDeepVSeMedium;
      vector<float> recoTauElectronCleanedDeepVSjetMedium;
      vector<float> recoTauElectronCleanedDeepVSmuMedium;

      vector<float> recoTauElectronCleanedDeepVSeTight;
      vector<float> recoTauElectronCleanedDeepVSjetTight;
      vector<float> recoTauElectronCleanedDeepVSmuTight;

      vector<float> recoTauElectronCleanedDeepVSeVLoose;
      vector<float> recoTauElectronCleanedDeepVSjetVLoose;
      vector<float> recoTauElectronCleanedDeepVSmuVLoose;

      vector<float> recoTauElectronCleanedDeepVSeVTight;
      vector<float> recoTauElectronCleanedDeepVSjetVTight;

      vector<float> recoTauElectronCleanedDeepVSeVVLoose;
      vector<float> recoTauElectronCleanedDeepVSjetVVLoose;

      vector<float> recoTauElectronCleanedDeepVSeVVTight;
      vector<float> recoTauElectronCleanedDeepVSjetVVTight;

      vector<float> recoTauElectronCleanedDeepVSeVVVLoose;
      vector<float> recoTauElectronCleanedDeepVSjetVVVLoose;

      vector<int> recoTauElectronCleanedRefToMuon;
      vector<int> recoTauElectronCleanedRefToElectron;

      // --- reconstructed taus boosted ---
      vector<float> recoTauBoostedPt;
      vector<float> recoTauBoostedEta;
      vector<float> recoTauBoostedPhi;
      vector<float> recoTauBoostedEnergy;
      vector<int> recoTauBoostedPDGId;
      vector<float> recoTauBoostedDecayMode;
      vector<float> recoTauBoostedDecayModeFinding;
      vector<float> recoTauBoostedDecayModeFindingNewDMs;
      vector<float> recoTauBoostedIsoMVArawValue;
      vector<float> recoTauBoostedIsoMVAVVLoose;
      vector<float> recoTauBoostedIsoMVAVLoose;
      vector<float> recoTauBoostedIsoMVALoose;
      vector<float> recoTauBoostedIsoMVAMedium;
      vector<float> recoTauBoostedIsoMVATight;
      vector<float> recoTauBoostedIsoMVAVTight;
      vector<float> recoTauBoostedIsoMVAVVTight;
      
      vector<float> recoTauBoostedAntiMuMVALoose;
      vector<float> recoTauBoostedAntiMuMVATight;

      vector<float> recoTauBoostedAntiEleMVArawValue;
      vector<float> recoTauBoostedAntiEleMVAVLoose;
      vector<float> recoTauBoostedAntiEleMVALoose;
      vector<float> recoTauBoostedAntiEleMVAMedium;
      vector<float> recoTauBoostedAntiEleMVATight;
      vector<float> recoTauBoostedAntiEleMVAVTight;

      // *** deep tau ID variables ***
      vector<float> recoTauBoostedDeepVSeraw;
      vector<float> recoTauBoostedDeepVSjetraw;
      vector<float> recoTauBoostedDeepVSmuraw;
 
      vector<float> recoTauBoostedDeepVSeLoose;
      vector<float> recoTauBoostedDeepVSjetLoose;
      vector<float> recoTauBoostedDeepVSmuLoose;

      vector<float> recoTauBoostedDeepVSeMedium;
      vector<float> recoTauBoostedDeepVSjetMedium;
      vector<float> recoTauBoostedDeepVSmuMedium;

      vector<float> recoTauBoostedDeepVSeTight;
      vector<float> recoTauBoostedDeepVSjetTight;
      vector<float> recoTauBoostedDeepVSmuTight;

      vector<float> recoTauBoostedDeepVSeVLoose;
      vector<float> recoTauBoostedDeepVSjetVLoose;
      vector<float> recoTauBoostedDeepVSmuVLoose;

      vector<float> recoTauBoostedDeepVSeVTight;
      vector<float> recoTauBoostedDeepVSjetVTight;

      vector<float> recoTauBoostedDeepVSeVVLoose;
      vector<float> recoTauBoostedDeepVSjetVVLoose;

      vector<float> recoTauBoostedDeepVSeVVTight;
      vector<float> recoTauBoostedDeepVSjetVVTight;

      vector<float> recoTauBoostedDeepVSeVVVLoose;
      vector<float> recoTauBoostedDeepVSjetVVVLoose;

      vector<int> recoTauBoostedRefToMuon;
      vector<int> recoTauBoostedRefToElectron;

      // --- reconstructed jets ---
      vector<float> recoJetPt;
      vector<float> recoJetEta;
      vector<float> recoJetPhi;
      vector<float> recoJetEnergy;
      vector<float> recoJetCSV;
      vector<float> recoJetDeepDiTauValue;
      vector<float> recoJetDeepDiTauValueMD;
      vector<int> recoJetIdLoose;
      vector<int> recoJetIdTight;
      vector<int> recoJetIdTightLepVeto;
      vector<int> recoJetIdPileUp;
      
      // --- reconstructed MET ---
      vector<float> recoMET;
      vector<float> recoMETPhi;
      vector<float> recoMETPx;
      vector<float> recoMETPy;

      // --- pileup and reconstructed vertices ---
      int recoNPrimaryVertex;
      int recoNPU;
      int trueNInteraction;
      int eventID;

      // --- gen muons ----
      vector<float> genMuonPt;
      vector<float> genMuonEta;
      vector<float> genMuonPhi;
      vector<float> genMuonMass;
      vector<int> genMuonPDGId;
      vector<int> genMuonMotherPDGId;

      // --- gen electrons ----
      vector<float> genElectronPt;
      vector<float> genElectronEta;
      vector<float> genElectronPhi;
      vector<float> genElectronMass;
      vector<int> genElectronPDGId;
      vector<int> genElectronMotherPDGId;

      // --- gen tau_mu ----
      vector<float> genTauMuPt;
      vector<float> genTauMuEta;
      vector<float> genTauMuPhi;
      vector<float> genTauMuMass;
      vector<int> genTauMuPDGId;
      vector<int> genTauMuMotherPDGId;
      vector<float> genTauMuVisPt;
      vector<float> genTauMuVisMass;

      // --- gen tau_e ----
      vector<float> genTauElePt;
      vector<float> genTauEleEta;
      vector<float> genTauElePhi;
      vector<float> genTauEleMass;
      vector<int> genTauElePDGId;
      vector<int> genTauEleMotherPDGId;
      vector<float> genTauEleVisPt;
      vector<float> genTauEleVisMass;

      // --- gen tau_h ----
      vector<float> genTauHadPt;
      vector<float> genTauHadEta;
      vector<float> genTauHadPhi;
      vector<float> genTauHadMass;
      vector<int> genTauHadPDGId;
      vector<int> genTauHadMotherPDGId;
      vector<float> genTauHadVisPt;
      vector<float> genTauHadVisMass;
      vector<int> genTauHadNPionZero;
      vector<int> genTauHadNChargedHadrons;

      // --- event weight for MC ---
      float genEventWeight; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiMuDiTauAnalyzer::DiMuDiTauAnalyzer(const edm::ParameterSet& iConfig):
    MuTag(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuTag"))),
    EleTag(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("EleTag"))),
    TauTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauTag"))),
    TauMuonCleanedTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauMuonCleanedTag"))),
    TauElectronCleanedTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauElectronCleanedTag"))),
    TauBoostedTag(consumes<edm::View<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauBoostedTag"))),
    JetTag(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    MetTag(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("MetTag"))),
    VertexTag(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag"))),
    rhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
    effectiveAreas((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
    PileupTag(consumes<edm::View<PileupSummaryInfo>>(iConfig.existsAs<edm::InputTag>("PileupTag") ? iConfig.getParameter<edm::InputTag>("PileupTag") : edm::InputTag())),
    generator(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ? iConfig.getParameter<edm::InputTag>("Generator") : edm::InputTag())),
    GenMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenMuTag") ? iConfig.getParameter<edm::InputTag>("GenMuTag") : edm::InputTag())),
    GenEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenEleTag") ? iConfig.getParameter<edm::InputTag>("GenEleTag") : edm::InputTag())),
    GenTauMuTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauMuTag") ? iConfig.getParameter<edm::InputTag>("GenTauMuTag") : edm::InputTag())),
    GenTauEleTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauEleTag") ? iConfig.getParameter<edm::InputTag>("GenTauEleTag") : edm::InputTag())),
    GenTauHadTag(consumes<edm::View<reco::GenParticle>>(iConfig.existsAs<edm::InputTag>("GenTauHadTag") ? iConfig.getParameter<edm::InputTag>("GenTauHadTag") : edm::InputTag()))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   isMC = iConfig.getParameter<bool>("isMC");
   numberOfTrigMus = iConfig.getParameter<int>("numberOfTrigMus");
}


DiMuDiTauAnalyzer::~DiMuDiTauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiMuDiTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu;
   iEvent.getByToken(MuTag, pMu);

   edm::Handle<edm::View<pat::Electron>> pElectron;
   iEvent.getByToken(EleTag, pElectron);

   edm::Handle<edm::View<pat::Tau>> pTau;
   iEvent.getByToken(TauTag, pTau);

   edm::Handle<edm::View<pat::Tau>> pTauMuonCleaned;
   iEvent.getByToken(TauMuonCleanedTag, pTauMuonCleaned);

   edm::Handle<edm::View<pat::Tau>> pTauElectronCleaned;
   iEvent.getByToken(TauElectronCleanedTag, pTauElectronCleaned);

   edm::Handle<edm::View<pat::Tau>> pTauBoosted;
   iEvent.getByToken(TauBoostedTag, pTauBoosted);

   edm::Handle<edm::View<pat::Jet>> pJet;
   iEvent.getByToken(JetTag, pJet);

   edm::Handle<edm::View<pat::MET>> pMet;
   iEvent.getByToken(MetTag, pMet);

   edm::Handle<edm::View<reco::Vertex>> pVertex;
   iEvent.getByToken(VertexTag, pVertex);

   edm::Handle<double> pRho;
   iEvent.getByToken(rhoTag, pRho);

   if (isMC)
   {
       edm::Handle<edm::View<reco::GenParticle>> pGenMu;
       iEvent.getByToken(GenMuTag, pGenMu);

       edm::Handle<edm::View<reco::GenParticle>> pGenEle;
       iEvent.getByToken(GenEleTag, pGenEle);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauMu;
       iEvent.getByToken(GenTauMuTag, pGenTauMu);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauEle;
       iEvent.getByToken(GenTauEleTag, pGenTauEle);

       edm::Handle<edm::View<reco::GenParticle>> pGenTauHad;
       iEvent.getByToken(GenTauHadTag, pGenTauHad);

       edm::Handle<GenEventInfoProduct> gen_ev_info;
       iEvent.getByToken(generator, gen_ev_info);

       edm::Handle<edm::View<PileupSummaryInfo>> pileup_info;
       iEvent.getByToken(PileupTag, pileup_info);

       if (pGenMu->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iMuon=pGenMu->begin(); iMuon!=pGenMu->end(); iMuon++)
           {
               genMuonPt.push_back(iMuon->pt());
               genMuonEta.push_back(iMuon->eta());
               genMuonPhi.push_back(iMuon->phi());
               genMuonMass.push_back(iMuon->mass());
               genMuonPDGId.push_back(iMuon->pdgId());
               genMuonMotherPDGId.push_back(iMuon->mother()->pdgId());
           } // end for loop on gen muons
       } // end if pGenMu->size() > 0

       if (pGenEle->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iElectron=pGenEle->begin(); iElectron!=pGenEle->end(); iElectron++)
           {
               genElectronPt.push_back(iElectron->pt());
               genElectronEta.push_back(iElectron->eta());
               genElectronPhi.push_back(iElectron->phi());
               genElectronMass.push_back(iElectron->mass());
               genElectronPDGId.push_back(iElectron->pdgId());
               genElectronMotherPDGId.push_back(iElectron->mother()->pdgId());
           } // end for loop on gen electrons
       } // end if pGenEle->size() > 0

       if (pGenTauMu->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauMu->begin(); iTau!=pGenTauMu->end(); iTau++)
           {
               genTauMuPt.push_back(iTau->pt());
               genTauMuEta.push_back(iTau->eta());
               genTauMuPhi.push_back(iTau->phi());
               genTauMuMass.push_back(iTau->mass());
               genTauMuPDGId.push_back(iTau->pdgId());
               genTauMuMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauMuVis;
               std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauMuEleVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauMuVis = sumTauMuVis + p4Daughter;
                   } // end for loop on all generations of visible daughter particles of tau_mu
               } // end for loop on tau_mu direct daughter particles

               genTauMuVisPt.push_back(sumTauMuVis.Pt());
               genTauMuVisMass.push_back(sumTauMuVis.M());
           } // end for loop on gen tau_mu
       } // end if pGenTauMu->size() > 0

       if (pGenTauEle->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauEle->begin(); iTau!=pGenTauEle->end(); iTau++)
           {
               genTauElePt.push_back(iTau->pt());
               genTauEleEta.push_back(iTau->eta());
               genTauElePhi.push_back(iTau->phi());
               genTauEleMass.push_back(iTau->mass());
               genTauElePDGId.push_back(iTau->pdgId());
               genTauEleMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauEleVis;
               std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauMuEleVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauEleVis = sumTauEleVis + p4Daughter;
                   } // end for loop on all generations of visible daughter particles of tau_e
               } // end for loop on tau_e direct daughter particles

               genTauEleVisPt.push_back(sumTauEleVis.Pt());
               genTauEleVisMass.push_back(sumTauEleVis.M());
           } // end for loop on gen tau_e
       } // end if pGenTauEle->size() > 0

       if (pGenTauHad->size() > 0)
       {
           for (edm::View<reco::GenParticle>::const_iterator iTau=pGenTauHad->begin(); iTau!=pGenTauHad->end(); iTau++)
           {
               genTauHadPt.push_back(iTau->pt());
               genTauHadEta.push_back(iTau->eta());
               genTauHadPhi.push_back(iTau->phi());
               genTauHadMass.push_back(iTau->mass());
               genTauHadPDGId.push_back(iTau->pdgId());
               genTauHadMotherPDGId.push_back(iTau->mother()->pdgId());

               TLorentzVector sumTauHadVis;
               std::vector <const reco::Candidate*> daughters;
               daughters.clear();

               int nPiZeros = 0;
               int nChargedHadrons = 0;

               for (unsigned int iDau = 0; iDau < iTau->numberOfDaughters(); iDau++)
               {
                   const reco::Candidate* directDaughter = iTau->daughter(iDau);
                   daughters = findTauHadVisDaughters(directDaughter); // collect all the current daughter (if status == 1) or together its daughters (if status != 1) 

                   nPiZeros += findTauPiZeros(directDaughter);
                   nChargedHadrons += findTauChargedHadrons(directDaughter);

                   for (unsigned int jDau = 0; jDau < daughters.size(); jDau++)
                   {
                       TLorentzVector p4Daughter;
                       p4Daughter.SetPtEtaPhiM(daughters[jDau]->pt(), daughters[jDau]->eta(), daughters[jDau]->phi(), daughters[jDau]->mass());
                       sumTauHadVis = sumTauHadVis + p4Daughter;
                   } // end for loop on all generations of visible daughter particles of tau_h
               } // end for loop on tau_h direct daughter particles

               genTauHadVisPt.push_back(sumTauHadVis.Pt());
               genTauHadVisMass.push_back(sumTauHadVis.M());
               genTauHadNPionZero.push_back(nPiZeros);
               genTauHadNChargedHadrons.push_back(nChargedHadrons);
           } // end for loop on gen tau_h
       } // end if pGenTauHad->size() > 0

       if (gen_ev_info.isValid())
       {
           genEventWeight = gen_ev_info->weight();
       } // end if gen_ev_info.isValid()

       if (pileup_info.isValid())
       {
           for(edm::View<PileupSummaryInfo>::const_iterator iPileup=pileup_info->begin(); iPileup!=pileup_info->end(); iPileup++)
           {
               if (iPileup->getBunchCrossing() == 0)
               {
                   trueNInteraction = iPileup->getTrueNumInteractions();
                   recoNPU = iPileup->getPU_NumInteractions();
               } // end if iPileup->getBunchCrossing() == 0
           } // end for loop on pileup_info
       } // end if pileup_info.isValid()
   } // end if isMC == true

   // --- prepare for offline primary vertices ---
   recoNPrimaryVertex = 0; 
   if (pVertex.isValid())
   {
       for(edm::View<reco::Vertex>::const_iterator iPV=pVertex->begin(); iPV!=pVertex->end(); iPV++)
       {
           recoNPrimaryVertex++;
       } // end for loop on pVertex
   } // end if pVertex.isValid()

   eventID = iEvent.eventAuxiliary().event();

   // --- prepare muon vector ---
   if (pMu->size()>0)
   {
       int muonCounter = 0;
       for(edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
       {
           recoMuonPt.push_back(iMuon->pt());
           recoMuonEta.push_back(iMuon->eta());
           recoMuonPhi.push_back(iMuon->phi());
           recoMuonEnergy.push_back(iMuon->energy());
           recoMuonPDGId.push_back(iMuon->pdgId());
           reco::MuonPFIsolation iso = iMuon->pfIsolationR04();
           double reliso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / iMuon->pt();
           recoMuonIsolation.push_back(reliso);
           recoMuonDXY.push_back(iMuon->muonBestTrack()->dxy());
           recoMuonDZ.push_back(iMuon->muonBestTrack()->dz());
           recoMuonNTrackerLayers.push_back(iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
           recoMuonRefToElectron.push_back(0);
           recoMuonRefToTau.push_back(0);
           recoMuonRefToTauMuonCleaned.push_back(0);

           bool goodGlob = iMuon->isGlobalMuon() && iMuon->globalTrack()->normalizedChi2() < 3 && iMuon->combinedQuality().chi2LocalPosition < 12 && iMuon->combinedQuality().trkKink < 20;
           int isMedium = muon::isLooseMuon(*iMuon) && iMuon->innerTrack()->validFraction() > 0.8 && muon::segmentCompatibility(*iMuon) > (goodGlob ? 0.303 : 0.451);
           int isLoose = iMuon->isPFMuon() && (iMuon->isGlobalMuon() || iMuon->isTrackerMuon());
           int isTight = iMuon->isGlobalMuon() && iMuon->isPFMuon() && iMuon->globalTrack()->normalizedChi2() < 10 && iMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && iMuon->numberOfMatchedStations() > 1 && fabs(iMuon->muonBestTrack()->dxy()) < 0.2 && fabs(iMuon->muonBestTrack()->dz()) < 0.5 && iMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;

           recoMuonIdLoose.push_back(isLoose);
           recoMuonIdMedium.push_back(isMedium);
           recoMuonIdTight.push_back(isTight);

           if (muonCounter == 0)
           {
               recoMuonTriggerFlag.push_back(1);
               muonCounter++;
           } // end if muonCounter == 0

           else if (muonCounter == 1 && numberOfTrigMus == 2)
           {
               recoMuonTriggerFlag.push_back(1);
               muonCounter++;
           } // end if muonCounter == 1 && numberOfTrigMus == 2

           else{
               recoMuonTriggerFlag.push_back(0);
               muonCounter++;
           } // end else if muonCounter != 0
       } // end for loop on muons
   } // end if pMu->size()>0

   // --- prepare electron vector ---
   if (pElectron->size()>0)
   {
       for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
       {
           int isLoose = 0;
           int isMedium = 0;
           int isTight = 0;

           int isLooseNoIso = 0;
           int isMediumNoIso = 0;
           int isTightNoIso = 0;

           // ---  full5x5_sigmaIetaIeta ---
           float sigmaIetaIeta = iElectron->full5x5_sigmaIetaIeta();

           // --- fabs(dEtaSeed) ---
           float dEtaSeed = fabs(iElectron->superCluster().isNonnull() && iElectron->superCluster()->seed().isNonnull() ? iElectron->deltaEtaSuperClusterTrackAtVtx() - iElectron->superCluster()->eta() + iElectron->superCluster()->seed()->eta() : std::numeric_limits<float>::max()); 
           
           // --- fabs(dPhiIn) ---
           float dPhiIn = fabs(iElectron->deltaPhiSuperClusterTrackAtVtx());
           
           // --- variables for H/E cuts ---
           float HoE = iElectron->hadronicOverEm();
           float rho = pRho.isValid() ? (*pRho) : 0; 
           float energy = iElectron->superCluster()->energy();

           // --- variables for relIsoWithEffectiveArea ---
           float chad = iElectron->pfIsolationVariables().sumChargedHadronPt;
           float nhad = iElectron->pfIsolationVariables().sumNeutralHadronEt;
           float pho = iElectron->pfIsolationVariables().sumPhotonEt;
           float elePt = iElectron->pt();
           float eleEta = iElectron->superCluster()->eta();
           float eArea = effectiveAreas.getEffectiveArea(fabs(eleEta));
           float relIsoWithEffectiveArea = (chad + std::max(0.0f, nhad + pho - rho*eArea)) / elePt;

           // --- variables for fabs(1/E-1/p) ---
           float eInverseMinusPInverse = fabs(1.0 - iElectron->eSuperClusterOverP())*(1.0/iElectron->ecalEnergy());

           // --- expected missing inner hits ---
           int mHits = iElectron->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

           // --- pass conversion veto ---
           bool isPassConVeto = iElectron->passConversionVeto();

           // ========= select electrons in different cut-based ID accordingly ==========
           if (fabs(eleEta) <= 1.479)
           {
               isLoose = (sigmaIetaIeta < 0.0112) &&
                         (dEtaSeed < 0.00377) &&
                         (dPhiIn < 0.0884) &&
                         (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.112 + 0.506/elePt) &&
                         (eInverseMinusPInverse < 0.193) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMedium = (sigmaIetaIeta < 0.0106) &&
                          (dEtaSeed < 0.0032) &&
                          (dPhiIn < 0.0547) &&
                          (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
                          (relIsoWithEffectiveArea < 0.0478 + 0.506/elePt) &&
                          (eInverseMinusPInverse < 0.184) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTight = (sigmaIetaIeta < 0.0104) &&
                         (dEtaSeed < 0.00255) &&
                         (dPhiIn < 0.022) &&
                         (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.0287 + 0.506/elePt) &&
                         (eInverseMinusPInverse < 0.159) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isLooseNoIso = (sigmaIetaIeta < 0.0112) &&
                         (dEtaSeed < 0.00377) &&
                         (dPhiIn < 0.0884) &&
                         (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
                         (eInverseMinusPInverse < 0.193) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMediumNoIso = (sigmaIetaIeta < 0.0106) &&
                          (dEtaSeed < 0.0032) &&
                          (dPhiIn < 0.0547) &&
                          (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
                          (eInverseMinusPInverse < 0.184) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTightNoIso = (sigmaIetaIeta < 0.0104) &&
                         (dEtaSeed < 0.00255) &&
                         (dPhiIn < 0.022) &&
                         (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
                         (eInverseMinusPInverse < 0.159) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);
           }// endif (fabs(eleEta) <= 1.479)

           else{
               isLoose = (sigmaIetaIeta < 0.0425) &&
                         (dEtaSeed < 0.00674) &&
                         (dPhiIn < 0.169) &&
                         (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.108 + 0.963/elePt) &&
                         (eInverseMinusPInverse < 0.111) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMedium = (sigmaIetaIeta < 0.0387) &&
                          (dEtaSeed < 0.00632) &&
                          (dPhiIn < 0.0394) &&
                          (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
                          (relIsoWithEffectiveArea < 0.0658 + 0.963/elePt) &&
                          (eInverseMinusPInverse < 0.0721) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTight = (sigmaIetaIeta < 0.0353) &&
                         (dEtaSeed < 0.00501) &&
                         (dPhiIn < 0.0236) &&
                         (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
                         (relIsoWithEffectiveArea < 0.0445 + 0.963/elePt) &&
                         (eInverseMinusPInverse < 0.0197) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isLooseNoIso = (sigmaIetaIeta < 0.0425) &&
                         (dEtaSeed < 0.00674) &&
                         (dPhiIn < 0.169) &&
                         (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
                         (eInverseMinusPInverse < 0.111) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);

               isMediumNoIso = (sigmaIetaIeta < 0.0387) &&
                          (dEtaSeed < 0.00632) &&
                          (dPhiIn < 0.0394) &&
                          (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
                          (eInverseMinusPInverse < 0.0721) &&
                          (mHits <= 1) &&
                          (isPassConVeto == true);

               isTightNoIso = (sigmaIetaIeta < 0.0353) &&
                         (dEtaSeed < 0.00501) &&
                         (dPhiIn < 0.0236) &&
                         (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
                         (eInverseMinusPInverse < 0.0197) &&
                         (mHits <= 1) &&
                         (isPassConVeto == true);
           } // end else (fabs(eleEta) > 1.479)

           recoElectronPt.push_back(iElectron->pt());
           recoElectronEta.push_back(iElectron->eta());
           recoElectronPhi.push_back(iElectron->phi());
           recoElectronEnergy.push_back(iElectron->energy());
           recoElectronPDGId.push_back(iElectron->pdgId());
           recoElectronIsolation.push_back(relIsoWithEffectiveArea);
           recoElectronIdLoose.push_back(isLoose);
           recoElectronIdMedium.push_back(isMedium);
           recoElectronIdTight.push_back(isTight);
           recoElectronIdLooseNoIso.push_back(isLooseNoIso);
           recoElectronIdMediumNoIso.push_back(isMediumNoIso);
           recoElectronIdTightNoIso.push_back(isTightNoIso);
           recoElectronEcalTrkEnergyPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyPostCorr"));
           recoElectronEcalTrkEnergyErrPostCorr.push_back(iElectron->userFloat("ecalTrkEnergyErrPostCorr"));
           recoElectronEnergyScaleValue.push_back(iElectron->userFloat("energyScaleValue"));
           recoElectronEnergyScaleUp.push_back(iElectron->userFloat("energyScaleUp"));
           recoElectronEnergyScaleDown.push_back(iElectron->userFloat("energyScaleDown"));
           recoElectronEnergySigmaValue.push_back(iElectron->userFloat("energySigmaValue"));
           recoElectronEnergySigmaUp.push_back(iElectron->userFloat("energySigmaUp"));
           recoElectronEnergySigmaDown.push_back(iElectron->userFloat("energySigmaDown"));

           recoElectronRefToMuon.push_back(0);
           recoElectronRefToTau.push_back(0);
           recoElectronRefToTauElectronCleaned.push_back(0);
       } // end for loop on electrons
   } // end if pElectron->size()>0
   
   // --- prepare tau vector ---
   if (pTau->size()>0)
   {
       for(edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           recoTauPt.push_back(iTau->pt());
           recoTauEta.push_back(iTau->eta());
           recoTauPhi.push_back(iTau->phi());
           recoTauEnergy.push_back(iTau->energy());
           recoTauPDGId.push_back(iTau->pdgId());
           recoTauDecayMode.push_back(iTau->decayMode());
           recoTauDecayModeFinding.push_back(iTau->tauID("decayModeFinding"));
           recoTauDecayModeFindingNewDMs.push_back(iTau->tauID("decayModeFindingNewDMs"));
           recoTauRefToMuon.push_back(0);
           recoTauRefToElectron.push_back(0);

           if (iTau->isTauIDAvailable("byIsolationMVArun2017v2DBoldDMwLTraw2017"))
           {
               recoTauIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
               recoTauIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           } // end if TauMVA ID available

           if (iTau->isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
           {
               recoTauDeepVSeraw.push_back(iTau->tauID("byDeepTau2017v2p1VSeraw"));
               recoTauDeepVSjetraw.push_back(iTau->tauID("byDeepTau2017v2p1VSjetraw"));
               recoTauDeepVSmuraw.push_back(iTau->tauID("byDeepTau2017v2p1VSmuraw"));

               recoTauDeepVSeLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSmu"));

               recoTauDeepVSeMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSe"));
               recoTauDeepVSjetMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSmu"));

               recoTauDeepVSeTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSe"));
               recoTauDeepVSjetTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSmu"));

               recoTauDeepVSeVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSjet"));
               recoTauDeepVSmuVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSmu"));

               recoTauDeepVSeVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSjet"));

               recoTauDeepVSeVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
               recoTauDeepVSjetVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
           } // end if DeepTau ID available

           if (iTau->isTauIDAvailable("againstMuonLoose3"))
           {
               recoTauAntiMuMVALoose.push_back(iTau->tauID("againstMuonLoose3"));
               recoTauAntiMuMVATight.push_back(iTau->tauID("againstMuonTight3"));
               recoTauAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw"));
               recoTauAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA6"));
               recoTauAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA6"));
               recoTauAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA6"));
               recoTauAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA6"));
               recoTauAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA6"));
           } // end if old anti-lep discs are available

           else{
               recoTauAntiMuMVALoose.push_back(iTau->tauID("againstMuonLooseSimple"));
               recoTauAntiMuMVATight.push_back(iTau->tauID("againstMuonTightSimple"));
               recoTauAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw2018"));
               recoTauAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA62018"));
               recoTauAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA62018"));
               recoTauAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA62018"));
               recoTauAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA62018"));
               recoTauAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA62018"));
           } // end if new anti-lep discs are available
       } // end for loop on taus
   } // end if pTau->size()>0

   // --- prepare tau vector (muon cleaned) ---
   if (pTauMuonCleaned->size()>0)
   {
       for(edm::View<pat::Tau>::const_iterator iTau=pTauMuonCleaned->begin(); iTau!=pTauMuonCleaned->end(); iTau++)
       {
           recoTauMuonCleanedPt.push_back(iTau->pt());
           recoTauMuonCleanedEta.push_back(iTau->eta());
           recoTauMuonCleanedPhi.push_back(iTau->phi());
           recoTauMuonCleanedEnergy.push_back(iTau->energy());
           recoTauMuonCleanedPDGId.push_back(iTau->pdgId());
           recoTauMuonCleanedDecayMode.push_back(iTau->decayMode());
           recoTauMuonCleanedDecayModeFinding.push_back(iTau->tauID("decayModeFinding"));
           recoTauMuonCleanedDecayModeFindingNewDMs.push_back(iTau->tauID("decayModeFindingNewDMs"));
           recoTauMuonCleanedRefToMuon.push_back(0);
           recoTauMuonCleanedRefToElectron.push_back(0);

           if (iTau->isTauIDAvailable("byIsolationMVArun2017v2DBoldDMwLTraw2017"))
           {
               recoTauMuonCleanedIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
               recoTauMuonCleanedIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauMuonCleanedIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           } // end if TauMVA ID available

           if (iTau->isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
           {
               recoTauMuonCleanedDeepVSeraw.push_back(iTau->tauID("byDeepTau2017v2p1VSeraw"));
               recoTauMuonCleanedDeepVSjetraw.push_back(iTau->tauID("byDeepTau2017v2p1VSjetraw"));
               recoTauMuonCleanedDeepVSmuraw.push_back(iTau->tauID("byDeepTau2017v2p1VSmuraw"));

               recoTauMuonCleanedDeepVSeLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSjet"));
               recoTauMuonCleanedDeepVSmuLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSmu"));

               recoTauMuonCleanedDeepVSeMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSjet"));
               recoTauMuonCleanedDeepVSmuMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSmu"));

               recoTauMuonCleanedDeepVSeTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSjet"));
               recoTauMuonCleanedDeepVSmuTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSmu"));

               recoTauMuonCleanedDeepVSeVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSjet"));
               recoTauMuonCleanedDeepVSmuVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSmu"));

               recoTauMuonCleanedDeepVSeVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSjet"));

               recoTauMuonCleanedDeepVSeVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSjet"));

               recoTauMuonCleanedDeepVSeVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSjet"));

               recoTauMuonCleanedDeepVSeVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
               recoTauMuonCleanedDeepVSjetVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
           } // end if DeepTau ID available

           if (iTau->isTauIDAvailable("againstMuonLoose3"))
           {
               recoTauMuonCleanedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLoose3"));
               recoTauMuonCleanedAntiMuMVATight.push_back(iTau->tauID("againstMuonTight3"));
               recoTauMuonCleanedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw"));
               recoTauMuonCleanedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA6"));
               recoTauMuonCleanedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA6"));
               recoTauMuonCleanedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA6"));
               recoTauMuonCleanedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA6"));
               recoTauMuonCleanedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA6"));
           } // end if old anti-lep discs are available

           else{
               recoTauMuonCleanedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLooseSimple"));
               recoTauMuonCleanedAntiMuMVATight.push_back(iTau->tauID("againstMuonTightSimple"));
               recoTauMuonCleanedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw2018"));
               recoTauMuonCleanedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA62018"));
               recoTauMuonCleanedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA62018"));
               recoTauMuonCleanedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA62018"));
               recoTauMuonCleanedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA62018"));
               recoTauMuonCleanedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA62018"));
           } // end if new anti-lep discs are available
       } // end for loop on taus
   } // end if pTauMuonCleaned->size()>0

   // --- prepare tau vector (electron cleaned) ---
   if (pTauElectronCleaned->size()>0)
   {
       for(edm::View<pat::Tau>::const_iterator iTau=pTauElectronCleaned->begin(); iTau!=pTauElectronCleaned->end(); iTau++)
       {
           recoTauElectronCleanedPt.push_back(iTau->pt());
           recoTauElectronCleanedEta.push_back(iTau->eta());
           recoTauElectronCleanedPhi.push_back(iTau->phi());
           recoTauElectronCleanedEnergy.push_back(iTau->energy());
           recoTauElectronCleanedPDGId.push_back(iTau->pdgId());
           recoTauElectronCleanedDecayMode.push_back(iTau->decayMode());
           recoTauElectronCleanedDecayModeFinding.push_back(iTau->tauID("decayModeFinding"));
           recoTauElectronCleanedDecayModeFindingNewDMs.push_back(iTau->tauID("decayModeFindingNewDMs"));
           recoTauElectronCleanedRefToMuon.push_back(0);
           recoTauElectronCleanedRefToElectron.push_back(0);

           if (iTau->isTauIDAvailable("byIsolationMVArun2017v2DBoldDMwLTraw2017"))
           {
               recoTauElectronCleanedIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
               recoTauElectronCleanedIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauElectronCleanedIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           } // end if TauMVA ID available

           if (iTau->isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
           {
               recoTauElectronCleanedDeepVSeraw.push_back(iTau->tauID("byDeepTau2017v2p1VSeraw"));
               recoTauElectronCleanedDeepVSjetraw.push_back(iTau->tauID("byDeepTau2017v2p1VSjetraw"));
               recoTauElectronCleanedDeepVSmuraw.push_back(iTau->tauID("byDeepTau2017v2p1VSmuraw"));

               recoTauElectronCleanedDeepVSeLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSjet"));
               recoTauElectronCleanedDeepVSmuLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSmu"));

               recoTauElectronCleanedDeepVSeMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSjet"));
               recoTauElectronCleanedDeepVSmuMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSmu"));

               recoTauElectronCleanedDeepVSeTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSjet"));
               recoTauElectronCleanedDeepVSmuTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSmu"));

               recoTauElectronCleanedDeepVSeVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSjet"));
               recoTauElectronCleanedDeepVSmuVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSmu"));

               recoTauElectronCleanedDeepVSeVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSjet"));

               recoTauElectronCleanedDeepVSeVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSjet"));

               recoTauElectronCleanedDeepVSeVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSjet"));

               recoTauElectronCleanedDeepVSeVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
               recoTauElectronCleanedDeepVSjetVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
           } // end if DeepTau ID available

           if (iTau->isTauIDAvailable("againstMuonLoose3"))
           {
               recoTauElectronCleanedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLoose3"));
               recoTauElectronCleanedAntiMuMVATight.push_back(iTau->tauID("againstMuonTight3"));
               recoTauElectronCleanedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw"));
               recoTauElectronCleanedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA6"));
               recoTauElectronCleanedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA6"));
               recoTauElectronCleanedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA6"));
               recoTauElectronCleanedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA6"));
               recoTauElectronCleanedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA6"));
           } // end if old anti-lep discs are available

           else{
               recoTauElectronCleanedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLooseSimple"));
               recoTauElectronCleanedAntiMuMVATight.push_back(iTau->tauID("againstMuonTightSimple"));
               recoTauElectronCleanedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw2018"));
               recoTauElectronCleanedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA62018"));
               recoTauElectronCleanedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA62018"));
               recoTauElectronCleanedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA62018"));
               recoTauElectronCleanedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA62018"));
               recoTauElectronCleanedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA62018"));
           } // end if new anti-lep discs are available
       } // end for loop on taus
   } // end if pTauElectronCleaned->size()>0

   // --- prepare tau vector (boosted) ---
   if (pTauBoosted->size()>0)
   {
       for(edm::View<pat::Tau>::const_iterator iTau=pTauBoosted->begin(); iTau!=pTauBoosted->end(); iTau++)
       {
           recoTauBoostedPt.push_back(iTau->pt());
           recoTauBoostedEta.push_back(iTau->eta());
           recoTauBoostedPhi.push_back(iTau->phi());
           recoTauBoostedEnergy.push_back(iTau->energy());
           recoTauBoostedPDGId.push_back(iTau->pdgId());
           recoTauBoostedDecayMode.push_back(iTau->decayMode());
           recoTauBoostedDecayModeFinding.push_back(iTau->tauID("decayModeFinding"));
           recoTauBoostedDecayModeFindingNewDMs.push_back(iTau->tauID("decayModeFindingNewDMs"));
           recoTauBoostedRefToMuon.push_back(0);
           recoTauBoostedRefToElectron.push_back(0);

           if (iTau->isTauIDAvailable("byIsolationMVArun2017v2DBoldDMwLTraw2017"))
           {
               recoTauBoostedIsoMVArawValue.push_back(iTau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
               recoTauBoostedIsoMVAVVLoose.push_back(iTau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVAVLoose.push_back(iTau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVALoose.push_back(iTau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVAMedium.push_back(iTau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVATight.push_back(iTau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVAVTight.push_back(iTau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
               recoTauBoostedIsoMVAVVTight.push_back(iTau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
           } // end if TauMVA ID available

           if (iTau->isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
           {
               recoTauBoostedDeepVSeraw.push_back(iTau->tauID("byDeepTau2017v2p1VSeraw"));
               recoTauBoostedDeepVSjetraw.push_back(iTau->tauID("byDeepTau2017v2p1VSjetraw"));
               recoTauBoostedDeepVSmuraw.push_back(iTau->tauID("byDeepTau2017v2p1VSmuraw"));

               recoTauBoostedDeepVSeLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSjet"));
               recoTauBoostedDeepVSmuLoose.push_back(iTau->tauID("byLooseDeepTau2017v2p1VSmu"));

               recoTauBoostedDeepVSeMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSjet"));
               recoTauBoostedDeepVSmuMedium.push_back(iTau->tauID("byMediumDeepTau2017v2p1VSmu"));

               recoTauBoostedDeepVSeTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSjet"));
               recoTauBoostedDeepVSmuTight.push_back(iTau->tauID("byTightDeepTau2017v2p1VSmu"));

               recoTauBoostedDeepVSeVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSjet"));
               recoTauBoostedDeepVSmuVLoose.push_back(iTau->tauID("byVLooseDeepTau2017v2p1VSmu"));

               recoTauBoostedDeepVSeVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetVTight.push_back(iTau->tauID("byVTightDeepTau2017v2p1VSjet"));

               recoTauBoostedDeepVSeVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetVVLoose.push_back(iTau->tauID("byVVLooseDeepTau2017v2p1VSjet"));

               recoTauBoostedDeepVSeVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetVVTight.push_back(iTau->tauID("byVVTightDeepTau2017v2p1VSjet"));

               recoTauBoostedDeepVSeVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
               recoTauBoostedDeepVSjetVVVLoose.push_back(iTau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
           } // end if DeepTau ID available

           if (iTau->isTauIDAvailable("againstMuonLoose3"))
           {
               recoTauBoostedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLoose3"));
               recoTauBoostedAntiMuMVATight.push_back(iTau->tauID("againstMuonTight3"));
               recoTauBoostedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw"));
               recoTauBoostedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA6"));
               recoTauBoostedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA6"));
               recoTauBoostedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA6"));
               recoTauBoostedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA6"));
               recoTauBoostedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA6"));
           } // end if old anti-lep discs are available

           else{
               recoTauBoostedAntiMuMVALoose.push_back(iTau->tauID("againstMuonLooseSimple"));
               recoTauBoostedAntiMuMVATight.push_back(iTau->tauID("againstMuonTightSimple"));
               recoTauBoostedAntiEleMVArawValue.push_back(iTau->tauID("againstElectronMVA6Raw2018"));
               recoTauBoostedAntiEleMVAVLoose.push_back(iTau->tauID("againstElectronVLooseMVA62018"));
               recoTauBoostedAntiEleMVALoose.push_back(iTau->tauID("againstElectronLooseMVA62018"));
               recoTauBoostedAntiEleMVAMedium.push_back(iTau->tauID("againstElectronMediumMVA62018"));
               recoTauBoostedAntiEleMVATight.push_back(iTau->tauID("againstElectronTightMVA62018"));
               recoTauBoostedAntiEleMVAVTight.push_back(iTau->tauID("againstElectronVTightMVA62018"));
           } // end if new anti-lep discs are available
       } // end for loop on taus
   } // end if pTauBoosted->size()>0

   // --- prepare for the common source particle reference records for muon/electron candidates
   if (pElectron->size()>0 && pMu->size()>0)
   {
       int refMuonValue = 1;
       int countMuon = 0;
       for (edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
       {
           bool findMatchedMu = false;

           int refElectronValue = refMuonValue;
           int countElectron = 0;
           for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
           {
               bool findMatchedEle = false;
               for (unsigned int indexMuon = 0; indexMuon < iMuon->numberOfSourceCandidatePtrs(); indexMuon++)
               {
                   for (unsigned int indexEle = 0; indexEle < iElectron->numberOfSourceCandidatePtrs(); indexEle++)
                   {
                       if (iElectron->sourceCandidatePtr(indexEle).key() == iMuon->sourceCandidatePtr(indexMuon).key())
                       {
                           findMatchedMu = true;
                           findMatchedEle = true;
                           break;
                       } // end if find the same source of electron and muon
                   } // end for loop on electron source particles
                   if (findMatchedEle) break;
               } // end for loop on muon source particles
               
               if (findMatchedEle)
               {
                   recoElectronRefToMuon.at(countElectron) = refElectronValue;
               } // end if findMatchedEle

               countElectron += 1;
           } // end for loop on electron candidates

           if (findMatchedMu)
           {
               recoMuonRefToElectron.at(countMuon) = refMuonValue;
               refMuonValue += 1;
           } // end if findMatchedMu

           countMuon += 1;
       } // end for loop on muon candidates
   } // end if pElectron->size()>0 && pMu->size()>0

   // --- prepare muon VS tau reference for overlapped mu/tau candidates
   if (pMu->size()>0 && pTau->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refMuonValue = refTauValue;
           int countMuon = 0;
           for (edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
           {
               bool findMatchedMu = false;
               for (unsigned int indexMu = 0; indexMu < iMuon->numberOfSourceCandidatePtrs(); indexMu++)
               {
                   for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
                   {
                       if (iMuon->sourceCandidatePtr(indexMu).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedMu = true;
                           findMatchedTau = true;
                           break;
                       } // end if the muon source and tau source have the same key
                   } // end for loop on tau source particles
                   if (findMatchedMu) break;
               } // end for loop on mu source particles

               if (findMatchedMu)
               {
                   recoMuonRefToTau.at(countMuon) = refMuonValue;
               } // end if findMatchedMu

               countMuon += 1;
           } // end for loop on muon candidates

           if (findMatchedTau)
           {
               recoTauRefToMuon.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pMu->size()>0 && pTau->size()>0

   // --- prepare muon VS tau (muon cleaned) reference for overlapped mu/tau candidates
   if (pMu->size()>0 && pTauMuonCleaned->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTauMuonCleaned->begin(); iTau!=pTauMuonCleaned->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refMuonValue = refTauValue;
           int countMuon = 0;
           for (edm::View<pat::Muon>::const_iterator iMuon=pMu->begin(); iMuon!=pMu->end(); iMuon++)
           {
               bool findMatchedMu = false;
               for (unsigned int indexMu = 0; indexMu < iMuon->numberOfSourceCandidatePtrs(); indexMu++)
               {
                   for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
                   {
                       if (iMuon->sourceCandidatePtr(indexMu).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedMu = true;
                           findMatchedTau = true;
                           break;
                       } // end if the muon source and tau source have the same key
                   } // end for loop on tau source particles
                   if (findMatchedMu) break;
               } // end for loop on mu source particles

               if (findMatchedMu)
               {
                   recoMuonRefToTauMuonCleaned.at(countMuon) = refMuonValue;
               } // end if findMatchedMu

               countMuon += 1;
           } // end for loop on muon candidates

           if (findMatchedTau)
           {
               recoTauMuonCleanedRefToMuon.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pMu->size()>0 && pTauMuonCleaned->size()>0

   // --- prepare electron VS tau reference for overlapped ele/tau candidates
   if (pElectron->size()>0 && pTau->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTau->begin(); iTau!=pTau->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refElectronValue = refTauValue;
           int countElectron = 0;
           for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
           {
               bool findMatchedEle = false;
               for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
               {
                   for (unsigned int indexEle = 0; indexEle < iElectron->numberOfSourceCandidatePtrs(); indexEle++)
                   {
                       if (iElectron->sourceCandidatePtr(indexEle).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedEle = true;
                           findMatchedTau = true;
                           break;
                       } // end if find the same source of electron and tau
                   } // end for loop on electron source particles
                   if (findMatchedEle) break;
               } // end for loop on tau source particles

               if (findMatchedEle)
               {
                   recoElectronRefToTau.at(countElectron) = refElectronValue;
               } // end if findMatchedEle

               countElectron += 1;
           } // end for loop on electron candidates

           if (findMatchedTau)
           {
               recoTauRefToElectron.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pElectron->size()>0 && pTau->size()>0

   // --- prepare electron VS tau (electron cleaned) reference for overlapped ele/tau candidates
   if (pElectron->size()>0 && pTauElectronCleaned->size()>0)
   {
       int refTauValue = 1;
       int countTau = 0;
       for (edm::View<pat::Tau>::const_iterator iTau=pTauElectronCleaned->begin(); iTau!=pTauElectronCleaned->end(); iTau++)
       {
           bool findMatchedTau = false;

           int refElectronValue = refTauValue;
           int countElectron = 0;
           for(edm::View<pat::Electron>::const_iterator iElectron=pElectron->begin(); iElectron!=pElectron->end(); iElectron++)
           {
               bool findMatchedEle = false;
               for (unsigned int indexTau = 0; indexTau < iTau->numberOfSourceCandidatePtrs(); indexTau++)
               {
                   for (unsigned int indexEle = 0; indexEle < iElectron->numberOfSourceCandidatePtrs(); indexEle++)
                   {
                       if (iElectron->sourceCandidatePtr(indexEle).key() == iTau->sourceCandidatePtr(indexTau).key())
                       {
                           findMatchedEle = true;
                           findMatchedTau = true;
                           break;
                       } // end if find the same source of electron and tau
                   } // end for loop on electron source particles
                   if (findMatchedEle) break;
               } // end for loop on tau source particles

               if (findMatchedEle)
               {
                   recoElectronRefToTauElectronCleaned.at(countElectron) = refElectronValue;
               } // end if findMatchedEle

               countElectron += 1;
           } // end for loop on electron candidates

           if (findMatchedTau)
           {
               recoTauElectronCleanedRefToElectron.at(countTau) = refTauValue;
               refTauValue += 1;
           } // end if findMatchedTau

           countTau += 1;
       } // end for loop on tau candidates
   } // end if pElectron->size()>0 && pTauElectronCleaned->size()>0

   // --- prepare jet vector ---
   if (pJet->size()>0)
   {
       for(edm::View<pat::Jet>::const_iterator iJet=pJet->begin(); iJet!=pJet->end(); iJet++)
       {
           if (iJet->pt() < 17 || fabs(iJet->eta()) > 2.5) continue;
           recoJetPt.push_back(iJet->pt());
           recoJetEta.push_back(iJet->eta());
           recoJetPhi.push_back(iJet->phi());
           recoJetEnergy.push_back(iJet->energy());
           // --- btag for jet ---
           // reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Jets
           // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Supported_Algorithms_and_Operati
           recoJetCSV.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
           recoJetDeepDiTauValue.push_back(iJet->userFloat("ditau2017v1"));
           recoJetDeepDiTauValueMD.push_back(iJet->userFloat("ditau2017MDv1"));
           recoJetIdLoose.push_back(iJet->userInt("idLoose"));
           recoJetIdTight.push_back(iJet->userInt("idTight"));
           recoJetIdTightLepVeto.push_back(iJet->userInt("idTightLepVeto"));
           recoJetIdPileUp.push_back(iJet->userInt("puID"));
       } // end for loop on jets
   } // end if pJet->size()>0

   // --- prepare MET vector ---
   if (pMet->size()>0)
   {
       for(edm::View<pat::MET>::const_iterator iMet=pMet->begin(); iMet!=pMet->end(); iMet++)
       {
           recoMET.push_back(iMet->pt());
           recoMETPhi.push_back(iMet->phi());
           recoMETPx.push_back(iMet->px());
           recoMETPy.push_back(iMet->py());
       } // end for loop on METs
   } // end if pMet->size()>0

   // --- fill the object tree ---
   objectTree->Fill();

   // ---- clear all the vectors for next event ----
   // --- reconstructed muons ---
   recoMuonPt.clear();
   recoMuonEta.clear();
   recoMuonPhi.clear();
   recoMuonEnergy.clear();
   recoMuonPDGId.clear();
   recoMuonIsolation.clear();
   recoMuonDXY.clear();
   recoMuonDZ.clear();
   recoMuonNTrackerLayers.clear();
   recoMuonTriggerFlag.clear();
   recoMuonRefToElectron.clear();
   recoMuonRefToTau.clear();
   recoMuonRefToTauMuonCleaned.clear();
   recoMuonIdLoose.clear();
   recoMuonIdMedium.clear();
   recoMuonIdTight.clear();

   // --- reconstructed electrons ---
   recoElectronPt.clear();
   recoElectronEta.clear();
   recoElectronPhi.clear();
   recoElectronEnergy.clear();
   recoElectronPDGId.clear();
   recoElectronIsolation.clear();
   recoElectronIdLoose.clear();
   recoElectronIdMedium.clear();
   recoElectronIdTight.clear();
   recoElectronIdLooseNoIso.clear();
   recoElectronIdMediumNoIso.clear();
   recoElectronIdTightNoIso.clear();
   recoElectronEcalTrkEnergyPostCorr.clear();
   recoElectronEcalTrkEnergyErrPostCorr.clear();
   recoElectronEnergyScaleValue.clear();
   recoElectronEnergyScaleUp.clear();
   recoElectronEnergyScaleDown.clear();
   recoElectronEnergySigmaValue.clear();
   recoElectronEnergySigmaUp.clear();
   recoElectronEnergySigmaDown.clear();
   recoElectronRefToMuon.clear();
   recoElectronRefToTau.clear();
   recoElectronRefToTauElectronCleaned.clear();

   // --- reconstructed taus ---
   recoTauPt.clear();
   recoTauEta.clear();
   recoTauPhi.clear();
   recoTauEnergy.clear();
   recoTauPDGId.clear();
   recoTauDecayMode.clear();
   recoTauDecayModeFinding.clear();
   recoTauDecayModeFindingNewDMs.clear();
   recoTauRefToMuon.clear();
   recoTauRefToElectron.clear();

   recoTauIsoMVArawValue.clear();
   recoTauIsoMVAVVLoose.clear();
   recoTauIsoMVAVLoose.clear();
   recoTauIsoMVALoose.clear();
   recoTauIsoMVAMedium.clear();
   recoTauIsoMVATight.clear();
   recoTauIsoMVAVTight.clear();
   recoTauIsoMVAVVTight.clear();
   
   recoTauAntiMuMVALoose.clear();
   recoTauAntiMuMVATight.clear();

   recoTauAntiEleMVArawValue.clear();
   recoTauAntiEleMVAVLoose.clear();
   recoTauAntiEleMVALoose.clear();
   recoTauAntiEleMVAMedium.clear();
   recoTauAntiEleMVATight.clear();
   recoTauAntiEleMVAVTight.clear();

   recoTauDeepVSeraw.clear();
   recoTauDeepVSjetraw.clear();
   recoTauDeepVSmuraw.clear();
 
   recoTauDeepVSeLoose.clear();
   recoTauDeepVSjetLoose.clear();
   recoTauDeepVSmuLoose.clear();

   recoTauDeepVSeMedium.clear();
   recoTauDeepVSjetMedium.clear();
   recoTauDeepVSmuMedium.clear();

   recoTauDeepVSeTight.clear();
   recoTauDeepVSjetTight.clear();
   recoTauDeepVSmuTight.clear();

   recoTauDeepVSeVLoose.clear();
   recoTauDeepVSjetVLoose.clear();
   recoTauDeepVSmuVLoose.clear();

   recoTauDeepVSeVTight.clear();
   recoTauDeepVSjetVTight.clear();

   recoTauDeepVSeVVLoose.clear();
   recoTauDeepVSjetVVLoose.clear();

   recoTauDeepVSeVVTight.clear();
   recoTauDeepVSjetVVTight.clear();

   recoTauDeepVSeVVVLoose.clear();
   recoTauDeepVSjetVVVLoose.clear();

   // --- reconstructed taus (muon cleaned) ---
   recoTauMuonCleanedPt.clear();
   recoTauMuonCleanedEta.clear();
   recoTauMuonCleanedPhi.clear();
   recoTauMuonCleanedEnergy.clear();
   recoTauMuonCleanedPDGId.clear();
   recoTauMuonCleanedDecayMode.clear();
   recoTauMuonCleanedDecayModeFinding.clear();
   recoTauMuonCleanedDecayModeFindingNewDMs.clear();
   recoTauMuonCleanedRefToMuon.clear();
   recoTauMuonCleanedRefToElectron.clear();

   recoTauMuonCleanedIsoMVArawValue.clear();
   recoTauMuonCleanedIsoMVAVVLoose.clear();
   recoTauMuonCleanedIsoMVAVLoose.clear();
   recoTauMuonCleanedIsoMVALoose.clear();
   recoTauMuonCleanedIsoMVAMedium.clear();
   recoTauMuonCleanedIsoMVATight.clear();
   recoTauMuonCleanedIsoMVAVTight.clear();
   recoTauMuonCleanedIsoMVAVVTight.clear();
   
   recoTauMuonCleanedAntiMuMVALoose.clear();
   recoTauMuonCleanedAntiMuMVATight.clear();

   recoTauMuonCleanedAntiEleMVArawValue.clear();
   recoTauMuonCleanedAntiEleMVAVLoose.clear();
   recoTauMuonCleanedAntiEleMVALoose.clear();
   recoTauMuonCleanedAntiEleMVAMedium.clear();
   recoTauMuonCleanedAntiEleMVATight.clear();
   recoTauMuonCleanedAntiEleMVAVTight.clear();

   recoTauMuonCleanedDeepVSeraw.clear();
   recoTauMuonCleanedDeepVSjetraw.clear();
   recoTauMuonCleanedDeepVSmuraw.clear();
 
   recoTauMuonCleanedDeepVSeLoose.clear();
   recoTauMuonCleanedDeepVSjetLoose.clear();
   recoTauMuonCleanedDeepVSmuLoose.clear();

   recoTauMuonCleanedDeepVSeMedium.clear();
   recoTauMuonCleanedDeepVSjetMedium.clear();
   recoTauMuonCleanedDeepVSmuMedium.clear();

   recoTauMuonCleanedDeepVSeTight.clear();
   recoTauMuonCleanedDeepVSjetTight.clear();
   recoTauMuonCleanedDeepVSmuTight.clear();

   recoTauMuonCleanedDeepVSeVLoose.clear();
   recoTauMuonCleanedDeepVSjetVLoose.clear();
   recoTauMuonCleanedDeepVSmuVLoose.clear();

   recoTauMuonCleanedDeepVSeVTight.clear();
   recoTauMuonCleanedDeepVSjetVTight.clear();

   recoTauMuonCleanedDeepVSeVVLoose.clear();
   recoTauMuonCleanedDeepVSjetVVLoose.clear();

   recoTauMuonCleanedDeepVSeVVTight.clear();
   recoTauMuonCleanedDeepVSjetVVTight.clear();

   recoTauMuonCleanedDeepVSeVVVLoose.clear();
   recoTauMuonCleanedDeepVSjetVVVLoose.clear();

   // --- reconstructed taus (electron cleaned) ---
   recoTauElectronCleanedPt.clear();
   recoTauElectronCleanedEta.clear();
   recoTauElectronCleanedPhi.clear();
   recoTauElectronCleanedEnergy.clear();
   recoTauElectronCleanedPDGId.clear();
   recoTauElectronCleanedDecayMode.clear();
   recoTauElectronCleanedDecayModeFinding.clear();
   recoTauElectronCleanedDecayModeFindingNewDMs.clear();
   recoTauElectronCleanedRefToMuon.clear();
   recoTauElectronCleanedRefToElectron.clear();

   recoTauElectronCleanedIsoMVArawValue.clear();
   recoTauElectronCleanedIsoMVAVVLoose.clear();
   recoTauElectronCleanedIsoMVAVLoose.clear();
   recoTauElectronCleanedIsoMVALoose.clear();
   recoTauElectronCleanedIsoMVAMedium.clear();
   recoTauElectronCleanedIsoMVATight.clear();
   recoTauElectronCleanedIsoMVAVTight.clear();
   recoTauElectronCleanedIsoMVAVVTight.clear();
   
   recoTauElectronCleanedAntiMuMVALoose.clear();
   recoTauElectronCleanedAntiMuMVATight.clear();

   recoTauElectronCleanedAntiEleMVArawValue.clear();
   recoTauElectronCleanedAntiEleMVAVLoose.clear();
   recoTauElectronCleanedAntiEleMVALoose.clear();
   recoTauElectronCleanedAntiEleMVAMedium.clear();
   recoTauElectronCleanedAntiEleMVATight.clear();
   recoTauElectronCleanedAntiEleMVAVTight.clear();

   recoTauElectronCleanedDeepVSeraw.clear();
   recoTauElectronCleanedDeepVSjetraw.clear();
   recoTauElectronCleanedDeepVSmuraw.clear();
 
   recoTauElectronCleanedDeepVSeLoose.clear();
   recoTauElectronCleanedDeepVSjetLoose.clear();
   recoTauElectronCleanedDeepVSmuLoose.clear();

   recoTauElectronCleanedDeepVSeMedium.clear();
   recoTauElectronCleanedDeepVSjetMedium.clear();
   recoTauElectronCleanedDeepVSmuMedium.clear();

   recoTauElectronCleanedDeepVSeTight.clear();
   recoTauElectronCleanedDeepVSjetTight.clear();
   recoTauElectronCleanedDeepVSmuTight.clear();

   recoTauElectronCleanedDeepVSeVLoose.clear();
   recoTauElectronCleanedDeepVSjetVLoose.clear();
   recoTauElectronCleanedDeepVSmuVLoose.clear();

   recoTauElectronCleanedDeepVSeVTight.clear();
   recoTauElectronCleanedDeepVSjetVTight.clear();

   recoTauElectronCleanedDeepVSeVVLoose.clear();
   recoTauElectronCleanedDeepVSjetVVLoose.clear();

   recoTauElectronCleanedDeepVSeVVTight.clear();
   recoTauElectronCleanedDeepVSjetVVTight.clear();

   recoTauElectronCleanedDeepVSeVVVLoose.clear();
   recoTauElectronCleanedDeepVSjetVVVLoose.clear();

   // --- reconstructed taus (boosted) ---
   recoTauBoostedPt.clear();
   recoTauBoostedEta.clear();
   recoTauBoostedPhi.clear();
   recoTauBoostedEnergy.clear();
   recoTauBoostedPDGId.clear();
   recoTauBoostedDecayMode.clear();
   recoTauBoostedDecayModeFinding.clear();
   recoTauBoostedDecayModeFindingNewDMs.clear();
   recoTauBoostedRefToMuon.clear();
   recoTauBoostedRefToElectron.clear();

   recoTauBoostedIsoMVArawValue.clear();
   recoTauBoostedIsoMVAVVLoose.clear();
   recoTauBoostedIsoMVAVLoose.clear();
   recoTauBoostedIsoMVALoose.clear();
   recoTauBoostedIsoMVAMedium.clear();
   recoTauBoostedIsoMVATight.clear();
   recoTauBoostedIsoMVAVTight.clear();
   recoTauBoostedIsoMVAVVTight.clear();
   
   recoTauBoostedAntiMuMVALoose.clear();
   recoTauBoostedAntiMuMVATight.clear();

   recoTauBoostedAntiEleMVArawValue.clear();
   recoTauBoostedAntiEleMVAVLoose.clear();
   recoTauBoostedAntiEleMVALoose.clear();
   recoTauBoostedAntiEleMVAMedium.clear();
   recoTauBoostedAntiEleMVATight.clear();
   recoTauBoostedAntiEleMVAVTight.clear();

   recoTauBoostedDeepVSeraw.clear();
   recoTauBoostedDeepVSjetraw.clear();
   recoTauBoostedDeepVSmuraw.clear();
 
   recoTauBoostedDeepVSeLoose.clear();
   recoTauBoostedDeepVSjetLoose.clear();
   recoTauBoostedDeepVSmuLoose.clear();

   recoTauBoostedDeepVSeMedium.clear();
   recoTauBoostedDeepVSjetMedium.clear();
   recoTauBoostedDeepVSmuMedium.clear();

   recoTauBoostedDeepVSeTight.clear();
   recoTauBoostedDeepVSjetTight.clear();
   recoTauBoostedDeepVSmuTight.clear();

   recoTauBoostedDeepVSeVLoose.clear();
   recoTauBoostedDeepVSjetVLoose.clear();
   recoTauBoostedDeepVSmuVLoose.clear();

   recoTauBoostedDeepVSeVTight.clear();
   recoTauBoostedDeepVSjetVTight.clear();

   recoTauBoostedDeepVSeVVLoose.clear();
   recoTauBoostedDeepVSjetVVLoose.clear();

   recoTauBoostedDeepVSeVVTight.clear();
   recoTauBoostedDeepVSjetVVTight.clear();

   recoTauBoostedDeepVSeVVVLoose.clear();
   recoTauBoostedDeepVSjetVVVLoose.clear();

   // --- reconstructed jets ---
   recoJetPt.clear();
   recoJetEta.clear();
   recoJetPhi.clear();
   recoJetEnergy.clear();
   recoJetCSV.clear();
   recoJetDeepDiTauValue.clear();
   recoJetDeepDiTauValueMD.clear();
   recoJetIdLoose.clear();
   recoJetIdTight.clear();
   recoJetIdTightLepVeto.clear();
   recoJetIdPileUp.clear();
   
   // --- reconstructed MET ---
   recoMET.clear();
   recoMETPhi.clear();
   recoMETPx.clear();
   recoMETPy.clear();

   // ---- gen muons ----
   genMuonPt.clear();
   genMuonEta.clear();
   genMuonPhi.clear();
   genMuonMass.clear();
   genMuonPDGId.clear();
   genMuonMotherPDGId.clear();

   // ---- gen electrons ----
   genElectronPt.clear();
   genElectronEta.clear();
   genElectronPhi.clear();
   genElectronMass.clear();
   genElectronPDGId.clear();
   genElectronMotherPDGId.clear();

   // ---- gen tau_mus ----
   genTauMuPt.clear();
   genTauMuEta.clear();
   genTauMuPhi.clear();
   genTauMuMass.clear();
   genTauMuPDGId.clear();
   genTauMuMotherPDGId.clear();
   genTauMuVisPt.clear();
   genTauMuVisMass.clear();

   // ---- gen tau_es ----
   genTauElePt.clear();
   genTauEleEta.clear();
   genTauElePhi.clear();
   genTauEleMass.clear();
   genTauElePDGId.clear();
   genTauEleMotherPDGId.clear();
   genTauEleVisPt.clear();
   genTauEleVisMass.clear();

   // ---- gen tau_hs ----
   genTauHadPt.clear();
   genTauHadEta.clear();
   genTauHadPhi.clear();
   genTauHadMass.clear();
   genTauHadPDGId.clear();
   genTauHadMotherPDGId.clear();
   genTauHadVisPt.clear();
   genTauHadVisMass.clear();
   genTauHadNPionZero.clear();
   genTauHadNChargedHadrons.clear();
}

// ------------ function for adding up all the visible daughter particles of tau_mu/tau_e ----------------
std::vector<const reco::Candidate*>
DiMuDiTauAnalyzer::findTauMuEleVisDaughters(const reco::Candidate* inputDaughter)
{
    std::vector<const reco::Candidate*> visParticles;
    if (inputDaughter->status() == 1)
    {
        if (fabs(inputDaughter->pdgId()) == 11 || fabs(inputDaughter->pdgId()) == 13 || fabs(inputDaughter->pdgId()) == 22)
        {
            visParticles.push_back(inputDaughter);
        } // end if final state particle is mu/ele/gamma
    } // end if input daughter is final state particle

    else{
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (grandDaughter->status() == 1)
            {
                if (fabs(grandDaughter->pdgId()) == 11 || fabs(grandDaughter->pdgId()) == 13 || fabs(grandDaughter->pdgId()) == 22)
                {
                    visParticles.push_back(grandDaughter);
                } // end if final state particle is mu/ele/gamma
            } // end if grand daughter is final state particle

            else{
                auto auxVisParticles = findTauMuEleVisDaughters(grandDaughter);
                visParticles.insert(visParticles.end(), auxVisParticles.begin(), auxVisParticles.end());
            } // end else grand daughter is(not) final state particle
        } // end for loop on grand daughters
    } // end else input daughter is(not) final state particle

    return visParticles;
}

// ------------ function for adding up all the visible daughter particles of tau_h ----------------
std::vector<const reco::Candidate*>
DiMuDiTauAnalyzer::findTauHadVisDaughters(const reco::Candidate* inputDaughter)
{
    std::vector<const reco::Candidate*> visParticles;
    if (inputDaughter->status() == 1)
    {
        if (fabs(inputDaughter->pdgId()) != 12 && fabs(inputDaughter->pdgId()) != 14 && fabs(inputDaughter->pdgId()) != 16)
        {
            visParticles.push_back(inputDaughter);
        } // end if final state particle is not neutrinos
    } // end if input daughter is final state particle

    else{
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (grandDaughter->status() == 1)
            {
                if (fabs(grandDaughter->pdgId()) != 12 && fabs(grandDaughter->pdgId()) != 14 && fabs(grandDaughter->pdgId()) != 16)
                {
                    visParticles.push_back(grandDaughter);
                } // end if final state particle is not neutrinos
            } // end if grand daughter is final state particle

            else{
                auto auxVisParticles = findTauHadVisDaughters(grandDaughter);
                visParticles.insert(visParticles.end(), auxVisParticles.begin(), auxVisParticles.end());
            } // end else grand daughter is(not) final state particle
        } // end for loop on grand daughters
    } // end else input daughter is(not) final state particle

    return visParticles;
}

// ------------ function for collecting all the pizeros from tau_h decay ----------------
int DiMuDiTauAnalyzer::findTauPiZeros(const reco::Candidate* inputDaughter)
{
    int numPiZero = 0;
    if (fabs(inputDaughter->pdgId()) == 111) numPiZero++;
    else if (fabs(inputDaughter->pdgId()) == 15)
    {
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            if (fabs(grandDaughter->pdgId()) == 111) numPiZero++;
            else if (fabs(grandDaughter->pdgId()) == 15)
            {
                numPiZero += findTauPiZeros(grandDaughter);
            } // end if the grand-daughter is still tau_h due to FSR
        } // end for loop on the grand-daughters
    } // end if the input daughter is still tau_h because of FSR

    return numPiZero;
}

// ------------ function for collecting all the charged hadrons from tau_h decay ----------------
int DiMuDiTauAnalyzer::findTauChargedHadrons(const reco::Candidate* inputDaughter)
{
    int numChargedHadrons = 0;
    bool chargedHadronsFromTau = fabs(inputDaughter->charge()) != 0 && fabs(inputDaughter->pdgId()) != 11 && fabs(inputDaughter->pdgId()) != 13 && fabs(inputDaughter->pdgId()) != 15; 

    if (chargedHadronsFromTau) numChargedHadrons++;
    else if (fabs(inputDaughter->pdgId()) == 15)
    {
        int nGrandDaughters = inputDaughter->numberOfDaughters();
        for (int iGrandDau = 0; iGrandDau < nGrandDaughters; iGrandDau++)
        {
            const reco::Candidate* grandDaughter = inputDaughter->daughter(iGrandDau);
            chargedHadronsFromTau = fabs(grandDaughter->charge()) != 0 && fabs(grandDaughter->pdgId()) != 11 && fabs(grandDaughter->pdgId()) != 13 && fabs(grandDaughter->pdgId()) != 15;
            if (chargedHadronsFromTau) numChargedHadrons++;
            else if (fabs(grandDaughter->pdgId()) == 15)
            {
                numChargedHadrons += findTauChargedHadrons(grandDaughter);
            } // end if the grand-daughter is still tau_h due to FSR
        } // end for loop on the grand-daughters
    } // end if the input daughter is still tau_h because of FSR

    return numChargedHadrons;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiMuDiTauAnalyzer::beginJob()
{
    // -- initialize the skimmed object vector tree --
    edm::Service<TFileService> fileService;
    objectTree = fileService->make<TTree>("objectTree","objectTree");

    objectTree->Branch("recoMuonPt", &recoMuonPt);
    objectTree->Branch("recoMuonEta", &recoMuonEta);
    objectTree->Branch("recoMuonPhi", &recoMuonPhi);
    objectTree->Branch("recoMuonEnergy", &recoMuonEnergy);
    objectTree->Branch("recoMuonPDGId", &recoMuonPDGId);
    objectTree->Branch("recoMuonIsolation", &recoMuonIsolation);
    objectTree->Branch("recoMuonDXY", &recoMuonDXY);
    objectTree->Branch("recoMuonDZ", &recoMuonDZ);
    objectTree->Branch("recoMuonNTrackerLayers", &recoMuonNTrackerLayers);
    objectTree->Branch("recoMuonTriggerFlag", &recoMuonTriggerFlag);
    objectTree->Branch("recoMuonRefToElectron", &recoMuonRefToElectron);
    objectTree->Branch("recoMuonRefToTau", &recoMuonRefToTau);
    objectTree->Branch("recoMuonRefToTauMuonCleaned", &recoMuonRefToTauMuonCleaned);
    objectTree->Branch("recoMuonIdLoose", &recoMuonIdLoose);
    objectTree->Branch("recoMuonIdMedium", &recoMuonIdMedium);
    objectTree->Branch("recoMuonIdTight", &recoMuonIdTight);

    objectTree->Branch("recoElectronPt", &recoElectronPt);
    objectTree->Branch("recoElectronEta", &recoElectronEta);
    objectTree->Branch("recoElectronPhi", &recoElectronPhi);
    objectTree->Branch("recoElectronEnergy", &recoElectronEnergy);
    objectTree->Branch("recoElectronPDGId", &recoElectronPDGId);
    objectTree->Branch("recoElectronIsolation", &recoElectronIsolation);
    objectTree->Branch("recoElectronIdLoose", &recoElectronIdLoose);
    objectTree->Branch("recoElectronIdMedium", &recoElectronIdMedium);
    objectTree->Branch("recoElectronIdTight", &recoElectronIdTight);
    objectTree->Branch("recoElectronIdLooseNoIso", &recoElectronIdLooseNoIso);
    objectTree->Branch("recoElectronIdMediumNoIso", &recoElectronIdMediumNoIso);
    objectTree->Branch("recoElectronIdTightNoIso", &recoElectronIdTightNoIso);
    objectTree->Branch("recoElectronEcalTrkEnergyPostCorr", &recoElectronEcalTrkEnergyPostCorr);
    objectTree->Branch("recoElectronEcalTrkEnergyErrPostCorr", &recoElectronEcalTrkEnergyErrPostCorr);
    objectTree->Branch("recoElectronEnergyScaleValue", &recoElectronEnergyScaleValue);
    objectTree->Branch("recoElectronEnergyScaleUp", &recoElectronEnergyScaleUp);
    objectTree->Branch("recoElectronEnergyScaleDown", &recoElectronEnergyScaleDown);
    objectTree->Branch("recoElectronEnergySigmaValue", &recoElectronEnergySigmaValue);
    objectTree->Branch("recoElectronEnergySigmaUp", &recoElectronEnergySigmaUp);
    objectTree->Branch("recoElectronEnergySigmaDown", &recoElectronEnergySigmaDown);
    objectTree->Branch("recoElectronRefToMuon", &recoElectronRefToMuon);
    objectTree->Branch("recoElectronRefToTau", &recoElectronRefToTau);
    objectTree->Branch("recoElectronRefToTauElectronCleaned", &recoElectronRefToTauElectronCleaned);

    objectTree->Branch("recoTauPt", &recoTauPt);
    objectTree->Branch("recoTauEta", &recoTauEta);
    objectTree->Branch("recoTauPhi", &recoTauPhi);
    objectTree->Branch("recoTauEnergy", &recoTauEnergy);
    objectTree->Branch("recoTauPDGId", &recoTauPDGId);
    objectTree->Branch("recoTauDecayMode", &recoTauDecayMode);
    objectTree->Branch("recoTauDecayModeFinding", &recoTauDecayModeFinding);
    objectTree->Branch("recoTauDecayModeFindingNewDMs", &recoTauDecayModeFindingNewDMs);
    objectTree->Branch("recoTauRefToMuon", &recoTauRefToMuon);
    objectTree->Branch("recoTauRefToElectron", &recoTauRefToElectron);

    objectTree->Branch("recoTauDeepVSeraw", &recoTauDeepVSeraw);
    objectTree->Branch("recoTauDeepVSjetraw", &recoTauDeepVSjetraw);
    objectTree->Branch("recoTauDeepVSmuraw", &recoTauDeepVSmuraw);

    objectTree->Branch("recoTauDeepVSeLoose", &recoTauDeepVSeLoose);
    objectTree->Branch("recoTauDeepVSjetLoose", &recoTauDeepVSjetLoose);
    objectTree->Branch("recoTauDeepVSmuLoose", &recoTauDeepVSmuLoose);

    objectTree->Branch("recoTauDeepVSeMedium", &recoTauDeepVSeMedium);
    objectTree->Branch("recoTauDeepVSjetMedium", &recoTauDeepVSjetMedium);
    objectTree->Branch("recoTauDeepVSmuMedium", &recoTauDeepVSmuMedium);

    objectTree->Branch("recoTauDeepVSeTight", &recoTauDeepVSeTight);
    objectTree->Branch("recoTauDeepVSjetTight", &recoTauDeepVSjetTight);
    objectTree->Branch("recoTauDeepVSmuTight", &recoTauDeepVSmuTight);

    objectTree->Branch("recoTauDeepVSeVLoose", &recoTauDeepVSeVLoose);
    objectTree->Branch("recoTauDeepVSjetVLoose", &recoTauDeepVSjetVLoose);
    objectTree->Branch("recoTauDeepVSmuVLoose", &recoTauDeepVSmuVLoose);

    objectTree->Branch("recoTauDeepVSeVTight", &recoTauDeepVSeVTight);
    objectTree->Branch("recoTauDeepVSjetVTight", &recoTauDeepVSjetVTight);

    objectTree->Branch("recoTauDeepVSeVVLoose", &recoTauDeepVSeVVLoose);
    objectTree->Branch("recoTauDeepVSjetVVLoose", &recoTauDeepVSjetVVLoose);

    objectTree->Branch("recoTauDeepVSeVVTight", &recoTauDeepVSeVVTight);
    objectTree->Branch("recoTauDeepVSjetVVTight", &recoTauDeepVSjetVVTight);

    objectTree->Branch("recoTauDeepVSeVVVLoose", &recoTauDeepVSeVVVLoose);
    objectTree->Branch("recoTauDeepVSjetVVVLoose", &recoTauDeepVSjetVVVLoose);

    objectTree->Branch("recoTauIsoMVArawValue", &recoTauIsoMVArawValue);
    objectTree->Branch("recoTauIsoMVAVVLoose", &recoTauIsoMVAVVLoose);
    objectTree->Branch("recoTauIsoMVAVLoose", &recoTauIsoMVAVLoose);
    objectTree->Branch("recoTauIsoMVALoose", &recoTauIsoMVALoose);
    objectTree->Branch("recoTauIsoMVAMedium", &recoTauIsoMVAMedium);
    objectTree->Branch("recoTauIsoMVATight", &recoTauIsoMVATight);
    objectTree->Branch("recoTauIsoMVAVTight", &recoTauIsoMVAVTight);
    objectTree->Branch("recoTauIsoMVAVVTight", &recoTauIsoMVAVVTight);

    objectTree->Branch("recoTauAntiMuMVALoose", &recoTauAntiMuMVALoose);
    objectTree->Branch("recoTauAntiMuMVATight", &recoTauAntiMuMVATight);
    
    objectTree->Branch("recoTauAntiEleMVArawValue", &recoTauAntiEleMVArawValue);
    objectTree->Branch("recoTauAntiEleMVAVLoose", &recoTauAntiEleMVAVLoose);
    objectTree->Branch("recoTauAntiEleMVALoose", &recoTauAntiEleMVALoose);
    objectTree->Branch("recoTauAntiEleMVAMedium", &recoTauAntiEleMVAMedium);
    objectTree->Branch("recoTauAntiEleMVATight", &recoTauAntiEleMVATight);
    objectTree->Branch("recoTauAntiEleMVAVTight", &recoTauAntiEleMVAVTight);

    objectTree->Branch("recoTauMuonCleanedPt", &recoTauMuonCleanedPt);
    objectTree->Branch("recoTauMuonCleanedEta", &recoTauMuonCleanedEta);
    objectTree->Branch("recoTauMuonCleanedPhi", &recoTauMuonCleanedPhi);
    objectTree->Branch("recoTauMuonCleanedEnergy", &recoTauMuonCleanedEnergy);
    objectTree->Branch("recoTauMuonCleanedPDGId", &recoTauMuonCleanedPDGId);
    objectTree->Branch("recoTauMuonCleanedDecayMode", &recoTauMuonCleanedDecayMode);
    objectTree->Branch("recoTauMuonCleanedDecayModeFinding", &recoTauMuonCleanedDecayModeFinding);
    objectTree->Branch("recoTauMuonCleanedDecayModeFindingNewDMs", &recoTauMuonCleanedDecayModeFindingNewDMs);
    objectTree->Branch("recoTauMuonCleanedRefToMuon", &recoTauMuonCleanedRefToMuon);
    objectTree->Branch("recoTauMuonCleanedRefToElectron", &recoTauMuonCleanedRefToElectron);

    objectTree->Branch("recoTauMuonCleanedDeepVSeraw", &recoTauMuonCleanedDeepVSeraw);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetraw", &recoTauMuonCleanedDeepVSjetraw);
    objectTree->Branch("recoTauMuonCleanedDeepVSmuraw", &recoTauMuonCleanedDeepVSmuraw);

    objectTree->Branch("recoTauMuonCleanedDeepVSeLoose", &recoTauMuonCleanedDeepVSeLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetLoose", &recoTauMuonCleanedDeepVSjetLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSmuLoose", &recoTauMuonCleanedDeepVSmuLoose);

    objectTree->Branch("recoTauMuonCleanedDeepVSeMedium", &recoTauMuonCleanedDeepVSeMedium);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetMedium", &recoTauMuonCleanedDeepVSjetMedium);
    objectTree->Branch("recoTauMuonCleanedDeepVSmuMedium", &recoTauMuonCleanedDeepVSmuMedium);

    objectTree->Branch("recoTauMuonCleanedDeepVSeTight", &recoTauMuonCleanedDeepVSeTight);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetTight", &recoTauMuonCleanedDeepVSjetTight);
    objectTree->Branch("recoTauMuonCleanedDeepVSmuTight", &recoTauMuonCleanedDeepVSmuTight);

    objectTree->Branch("recoTauMuonCleanedDeepVSeVLoose", &recoTauMuonCleanedDeepVSeVLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetVLoose", &recoTauMuonCleanedDeepVSjetVLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSmuVLoose", &recoTauMuonCleanedDeepVSmuVLoose);

    objectTree->Branch("recoTauMuonCleanedDeepVSeVTight", &recoTauMuonCleanedDeepVSeVTight);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetVTight", &recoTauMuonCleanedDeepVSjetVTight);

    objectTree->Branch("recoTauMuonCleanedDeepVSeVVLoose", &recoTauMuonCleanedDeepVSeVVLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetVVLoose", &recoTauMuonCleanedDeepVSjetVVLoose);

    objectTree->Branch("recoTauMuonCleanedDeepVSeVVTight", &recoTauMuonCleanedDeepVSeVVTight);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetVVTight", &recoTauMuonCleanedDeepVSjetVVTight);

    objectTree->Branch("recoTauMuonCleanedDeepVSeVVVLoose", &recoTauMuonCleanedDeepVSeVVVLoose);
    objectTree->Branch("recoTauMuonCleanedDeepVSjetVVVLoose", &recoTauMuonCleanedDeepVSjetVVVLoose);

    objectTree->Branch("recoTauMuonCleanedIsoMVArawValue", &recoTauMuonCleanedIsoMVArawValue);
    objectTree->Branch("recoTauMuonCleanedIsoMVAVVLoose", &recoTauMuonCleanedIsoMVAVVLoose);
    objectTree->Branch("recoTauMuonCleanedIsoMVAVLoose", &recoTauMuonCleanedIsoMVAVLoose);
    objectTree->Branch("recoTauMuonCleanedIsoMVALoose", &recoTauMuonCleanedIsoMVALoose);
    objectTree->Branch("recoTauMuonCleanedIsoMVAMedium", &recoTauMuonCleanedIsoMVAMedium);
    objectTree->Branch("recoTauMuonCleanedIsoMVATight", &recoTauMuonCleanedIsoMVATight);
    objectTree->Branch("recoTauMuonCleanedIsoMVAVTight", &recoTauMuonCleanedIsoMVAVTight);
    objectTree->Branch("recoTauMuonCleanedIsoMVAVVTight", &recoTauMuonCleanedIsoMVAVVTight);

    objectTree->Branch("recoTauMuonCleanedAntiMuMVALoose", &recoTauMuonCleanedAntiMuMVALoose);
    objectTree->Branch("recoTauMuonCleanedAntiMuMVATight", &recoTauMuonCleanedAntiMuMVATight);
    
    objectTree->Branch("recoTauMuonCleanedAntiEleMVArawValue", &recoTauMuonCleanedAntiEleMVArawValue);
    objectTree->Branch("recoTauMuonCleanedAntiEleMVAVLoose", &recoTauMuonCleanedAntiEleMVAVLoose);
    objectTree->Branch("recoTauMuonCleanedAntiEleMVALoose", &recoTauMuonCleanedAntiEleMVALoose);
    objectTree->Branch("recoTauMuonCleanedAntiEleMVAMedium", &recoTauMuonCleanedAntiEleMVAMedium);
    objectTree->Branch("recoTauMuonCleanedAntiEleMVATight", &recoTauMuonCleanedAntiEleMVATight);
    objectTree->Branch("recoTauMuonCleanedAntiEleMVAVTight", &recoTauMuonCleanedAntiEleMVAVTight);

    objectTree->Branch("recoTauElectronCleanedPt", &recoTauElectronCleanedPt);
    objectTree->Branch("recoTauElectronCleanedEta", &recoTauElectronCleanedEta);
    objectTree->Branch("recoTauElectronCleanedPhi", &recoTauElectronCleanedPhi);
    objectTree->Branch("recoTauElectronCleanedEnergy", &recoTauElectronCleanedEnergy);
    objectTree->Branch("recoTauElectronCleanedPDGId", &recoTauElectronCleanedPDGId);
    objectTree->Branch("recoTauElectronCleanedDecayMode", &recoTauElectronCleanedDecayMode);
    objectTree->Branch("recoTauElectronCleanedDecayModeFinding", &recoTauElectronCleanedDecayModeFinding);
    objectTree->Branch("recoTauElectronCleanedDecayModeFindingNewDMs", &recoTauElectronCleanedDecayModeFindingNewDMs);
    objectTree->Branch("recoTauElectronCleanedRefToMuon", &recoTauElectronCleanedRefToMuon);
    objectTree->Branch("recoTauElectronCleanedRefToElectron", &recoTauElectronCleanedRefToElectron);

    objectTree->Branch("recoTauElectronCleanedDeepVSeraw", &recoTauElectronCleanedDeepVSeraw);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetraw", &recoTauElectronCleanedDeepVSjetraw);
    objectTree->Branch("recoTauElectronCleanedDeepVSmuraw", &recoTauElectronCleanedDeepVSmuraw);

    objectTree->Branch("recoTauElectronCleanedDeepVSeLoose", &recoTauElectronCleanedDeepVSeLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetLoose", &recoTauElectronCleanedDeepVSjetLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSmuLoose", &recoTauElectronCleanedDeepVSmuLoose);

    objectTree->Branch("recoTauElectronCleanedDeepVSeMedium", &recoTauElectronCleanedDeepVSeMedium);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetMedium", &recoTauElectronCleanedDeepVSjetMedium);
    objectTree->Branch("recoTauElectronCleanedDeepVSmuMedium", &recoTauElectronCleanedDeepVSmuMedium);

    objectTree->Branch("recoTauElectronCleanedDeepVSeTight", &recoTauElectronCleanedDeepVSeTight);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetTight", &recoTauElectronCleanedDeepVSjetTight);
    objectTree->Branch("recoTauElectronCleanedDeepVSmuTight", &recoTauElectronCleanedDeepVSmuTight);

    objectTree->Branch("recoTauElectronCleanedDeepVSeVLoose", &recoTauElectronCleanedDeepVSeVLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetVLoose", &recoTauElectronCleanedDeepVSjetVLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSmuVLoose", &recoTauElectronCleanedDeepVSmuVLoose);

    objectTree->Branch("recoTauElectronCleanedDeepVSeVTight", &recoTauElectronCleanedDeepVSeVTight);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetVTight", &recoTauElectronCleanedDeepVSjetVTight);

    objectTree->Branch("recoTauElectronCleanedDeepVSeVVLoose", &recoTauElectronCleanedDeepVSeVVLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetVVLoose", &recoTauElectronCleanedDeepVSjetVVLoose);

    objectTree->Branch("recoTauElectronCleanedDeepVSeVVTight", &recoTauElectronCleanedDeepVSeVVTight);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetVVTight", &recoTauElectronCleanedDeepVSjetVVTight);

    objectTree->Branch("recoTauElectronCleanedDeepVSeVVVLoose", &recoTauElectronCleanedDeepVSeVVVLoose);
    objectTree->Branch("recoTauElectronCleanedDeepVSjetVVVLoose", &recoTauElectronCleanedDeepVSjetVVVLoose);

    objectTree->Branch("recoTauElectronCleanedIsoMVArawValue", &recoTauElectronCleanedIsoMVArawValue);
    objectTree->Branch("recoTauElectronCleanedIsoMVAVVLoose", &recoTauElectronCleanedIsoMVAVVLoose);
    objectTree->Branch("recoTauElectronCleanedIsoMVAVLoose", &recoTauElectronCleanedIsoMVAVLoose);
    objectTree->Branch("recoTauElectronCleanedIsoMVALoose", &recoTauElectronCleanedIsoMVALoose);
    objectTree->Branch("recoTauElectronCleanedIsoMVAMedium", &recoTauElectronCleanedIsoMVAMedium);
    objectTree->Branch("recoTauElectronCleanedIsoMVATight", &recoTauElectronCleanedIsoMVATight);
    objectTree->Branch("recoTauElectronCleanedIsoMVAVTight", &recoTauElectronCleanedIsoMVAVTight);
    objectTree->Branch("recoTauElectronCleanedIsoMVAVVTight", &recoTauElectronCleanedIsoMVAVVTight);

    objectTree->Branch("recoTauElectronCleanedAntiMuMVALoose", &recoTauElectronCleanedAntiMuMVALoose);
    objectTree->Branch("recoTauElectronCleanedAntiMuMVATight", &recoTauElectronCleanedAntiMuMVATight);
    
    objectTree->Branch("recoTauElectronCleanedAntiEleMVArawValue", &recoTauElectronCleanedAntiEleMVArawValue);
    objectTree->Branch("recoTauElectronCleanedAntiEleMVAVLoose", &recoTauElectronCleanedAntiEleMVAVLoose);
    objectTree->Branch("recoTauElectronCleanedAntiEleMVALoose", &recoTauElectronCleanedAntiEleMVALoose);
    objectTree->Branch("recoTauElectronCleanedAntiEleMVAMedium", &recoTauElectronCleanedAntiEleMVAMedium);
    objectTree->Branch("recoTauElectronCleanedAntiEleMVATight", &recoTauElectronCleanedAntiEleMVATight);
    objectTree->Branch("recoTauElectronCleanedAntiEleMVAVTight", &recoTauElectronCleanedAntiEleMVAVTight);

    objectTree->Branch("recoTauBoostedPt", &recoTauBoostedPt);
    objectTree->Branch("recoTauBoostedEta", &recoTauBoostedEta);
    objectTree->Branch("recoTauBoostedPhi", &recoTauBoostedPhi);
    objectTree->Branch("recoTauBoostedEnergy", &recoTauBoostedEnergy);
    objectTree->Branch("recoTauBoostedPDGId", &recoTauBoostedPDGId);
    objectTree->Branch("recoTauBoostedDecayMode", &recoTauBoostedDecayMode);
    objectTree->Branch("recoTauBoostedDecayModeFinding", &recoTauBoostedDecayModeFinding);
    objectTree->Branch("recoTauBoostedDecayModeFindingNewDMs", &recoTauBoostedDecayModeFindingNewDMs);
    objectTree->Branch("recoTauBoostedRefToMuon", &recoTauBoostedRefToMuon);
    objectTree->Branch("recoTauBoostedRefToElectron", &recoTauBoostedRefToElectron);

    objectTree->Branch("recoTauBoostedDeepVSeraw", &recoTauBoostedDeepVSeraw);
    objectTree->Branch("recoTauBoostedDeepVSjetraw", &recoTauBoostedDeepVSjetraw);
    objectTree->Branch("recoTauBoostedDeepVSmuraw", &recoTauBoostedDeepVSmuraw);

    objectTree->Branch("recoTauBoostedDeepVSeLoose", &recoTauBoostedDeepVSeLoose);
    objectTree->Branch("recoTauBoostedDeepVSjetLoose", &recoTauBoostedDeepVSjetLoose);
    objectTree->Branch("recoTauBoostedDeepVSmuLoose", &recoTauBoostedDeepVSmuLoose);

    objectTree->Branch("recoTauBoostedDeepVSeMedium", &recoTauBoostedDeepVSeMedium);
    objectTree->Branch("recoTauBoostedDeepVSjetMedium", &recoTauBoostedDeepVSjetMedium);
    objectTree->Branch("recoTauBoostedDeepVSmuMedium", &recoTauBoostedDeepVSmuMedium);

    objectTree->Branch("recoTauBoostedDeepVSeTight", &recoTauBoostedDeepVSeTight);
    objectTree->Branch("recoTauBoostedDeepVSjetTight", &recoTauBoostedDeepVSjetTight);
    objectTree->Branch("recoTauBoostedDeepVSmuTight", &recoTauBoostedDeepVSmuTight);

    objectTree->Branch("recoTauBoostedDeepVSeVLoose", &recoTauBoostedDeepVSeVLoose);
    objectTree->Branch("recoTauBoostedDeepVSjetVLoose", &recoTauBoostedDeepVSjetVLoose);
    objectTree->Branch("recoTauBoostedDeepVSmuVLoose", &recoTauBoostedDeepVSmuVLoose);

    objectTree->Branch("recoTauBoostedDeepVSeVTight", &recoTauBoostedDeepVSeVTight);
    objectTree->Branch("recoTauBoostedDeepVSjetVTight", &recoTauBoostedDeepVSjetVTight);

    objectTree->Branch("recoTauBoostedDeepVSeVVLoose", &recoTauBoostedDeepVSeVVLoose);
    objectTree->Branch("recoTauBoostedDeepVSjetVVLoose", &recoTauBoostedDeepVSjetVVLoose);

    objectTree->Branch("recoTauBoostedDeepVSeVVTight", &recoTauBoostedDeepVSeVVTight);
    objectTree->Branch("recoTauBoostedDeepVSjetVVTight", &recoTauBoostedDeepVSjetVVTight);

    objectTree->Branch("recoTauBoostedDeepVSeVVVLoose", &recoTauBoostedDeepVSeVVVLoose);
    objectTree->Branch("recoTauBoostedDeepVSjetVVVLoose", &recoTauBoostedDeepVSjetVVVLoose);

    objectTree->Branch("recoTauBoostedIsoMVArawValue", &recoTauBoostedIsoMVArawValue);
    objectTree->Branch("recoTauBoostedIsoMVAVVLoose", &recoTauBoostedIsoMVAVVLoose);
    objectTree->Branch("recoTauBoostedIsoMVAVLoose", &recoTauBoostedIsoMVAVLoose);
    objectTree->Branch("recoTauBoostedIsoMVALoose", &recoTauBoostedIsoMVALoose);
    objectTree->Branch("recoTauBoostedIsoMVAMedium", &recoTauBoostedIsoMVAMedium);
    objectTree->Branch("recoTauBoostedIsoMVATight", &recoTauBoostedIsoMVATight);
    objectTree->Branch("recoTauBoostedIsoMVAVTight", &recoTauBoostedIsoMVAVTight);
    objectTree->Branch("recoTauBoostedIsoMVAVVTight", &recoTauBoostedIsoMVAVVTight);

    objectTree->Branch("recoTauBoostedAntiMuMVALoose", &recoTauBoostedAntiMuMVALoose);
    objectTree->Branch("recoTauBoostedAntiMuMVATight", &recoTauBoostedAntiMuMVATight);
    
    objectTree->Branch("recoTauBoostedAntiEleMVArawValue", &recoTauBoostedAntiEleMVArawValue);
    objectTree->Branch("recoTauBoostedAntiEleMVAVLoose", &recoTauBoostedAntiEleMVAVLoose);
    objectTree->Branch("recoTauBoostedAntiEleMVALoose", &recoTauBoostedAntiEleMVALoose);
    objectTree->Branch("recoTauBoostedAntiEleMVAMedium", &recoTauBoostedAntiEleMVAMedium);
    objectTree->Branch("recoTauBoostedAntiEleMVATight", &recoTauBoostedAntiEleMVATight);
    objectTree->Branch("recoTauBoostedAntiEleMVAVTight", &recoTauBoostedAntiEleMVAVTight);

    objectTree->Branch("recoJetPt", &recoJetPt);
    objectTree->Branch("recoJetEta", &recoJetEta);
    objectTree->Branch("recoJetPhi", &recoJetPhi);
    objectTree->Branch("recoJetEnergy", &recoJetEnergy);
    objectTree->Branch("recoJetCSV", &recoJetCSV);
    objectTree->Branch("recoJetDeepDiTauValue", &recoJetDeepDiTauValue);
    objectTree->Branch("recoJetDeepDiTauValueMD", &recoJetDeepDiTauValueMD);
    objectTree->Branch("recoJetIdLoose", &recoJetIdLoose);
    objectTree->Branch("recoJetIdTight", &recoJetIdTight);
    objectTree->Branch("recoJetIdTightLepVeto", &recoJetIdTightLepVeto);
    objectTree->Branch("recoJetIdPileUp", &recoJetIdPileUp);
    
    objectTree->Branch("recoMET", &recoMET);
    objectTree->Branch("recoMETPhi", &recoMETPhi);
    objectTree->Branch("recoMETPx", &recoMETPx);
    objectTree->Branch("recoMETPy", &recoMETPy);

    objectTree->Branch("recoNPrimaryVertex", &recoNPrimaryVertex, "recoNPrimaryVertex/I");
    objectTree->Branch("eventID", &eventID, "eventID/I");

    if (isMC)
    {
        objectTree->Branch("genMuonPt", &genMuonPt);
        objectTree->Branch("genMuonEta", &genMuonEta);
        objectTree->Branch("genMuonPhi", &genMuonPhi);
        objectTree->Branch("genMuonMass", &genMuonMass);
        objectTree->Branch("genMuonPDGId", &genMuonPDGId);
        objectTree->Branch("genMuonMotherPDGId", &genMuonMotherPDGId);

        objectTree->Branch("genElectronPt", &genElectronPt);
        objectTree->Branch("genElectronEta", &genElectronEta);
        objectTree->Branch("genElectronPhi", &genElectronPhi);
        objectTree->Branch("genElectronMass", &genElectronMass);
        objectTree->Branch("genElectronPDGId", &genElectronPDGId);
        objectTree->Branch("genElectronMotherPDGId", &genElectronMotherPDGId);

        objectTree->Branch("genTauMuPt", &genTauMuPt);
        objectTree->Branch("genTauMuEta", &genTauMuEta);
        objectTree->Branch("genTauMuPhi", &genTauMuPhi);
        objectTree->Branch("genTauMuMass", &genTauMuMass);
        objectTree->Branch("genTauMuPDGId", &genTauMuPDGId);
        objectTree->Branch("genTauMuMotherPDGId", &genTauMuMotherPDGId);
        objectTree->Branch("genTauMuVisPt", &genTauMuVisPt);
        objectTree->Branch("genTauMuVisMass", &genTauMuVisMass);

        objectTree->Branch("genTauElePt", &genTauElePt);
        objectTree->Branch("genTauEleEta", &genTauEleEta);
        objectTree->Branch("genTauElePhi", &genTauElePhi);
        objectTree->Branch("genTauEleMass", &genTauEleMass);
        objectTree->Branch("genTauElePDGId", &genTauElePDGId);
        objectTree->Branch("genTauEleMotherPDGId", &genTauEleMotherPDGId);
        objectTree->Branch("genTauEleVisPt", &genTauEleVisPt);
        objectTree->Branch("genTauEleVisMass", &genTauEleVisMass);

        objectTree->Branch("genTauHadPt", &genTauHadPt);
        objectTree->Branch("genTauHadEta", &genTauHadEta);
        objectTree->Branch("genTauHadPhi", &genTauHadPhi);
        objectTree->Branch("genTauHadMass", &genTauHadMass);
        objectTree->Branch("genTauHadPDGId", &genTauHadPDGId);
        objectTree->Branch("genTauHadMotherPDGId", &genTauHadMotherPDGId);
        objectTree->Branch("genTauHadVisPt", &genTauHadVisPt);
        objectTree->Branch("genTauHadVisMass", &genTauHadVisMass);
        objectTree->Branch("genTauHadNPionZero", &genTauHadNPionZero);
        objectTree->Branch("genTauHadNChargedHadrons", &genTauHadNChargedHadrons);

        objectTree->Branch("recoNPU", &recoNPU, "recoNPU/I");
        objectTree->Branch("trueNInteraction", &trueNInteraction, "trueNInteraction/I");
        objectTree->Branch("genEventWeight", &genEventWeight, "genEventWeight/F");
    } // end if isMC == true
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuDiTauAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuDiTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuDiTauAnalyzer);
