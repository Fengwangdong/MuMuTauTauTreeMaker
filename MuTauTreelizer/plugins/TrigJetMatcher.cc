// -*- C++ -*-
//
// Package:    testED/TrigJetMatcher
// Class:      TrigJetMatcher
// 
/**\class TrigJetMatcher TrigJetMatcher.cc testED/TrigJetMatcher/plugins/TrigJetMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fengwangdong Zhang
//         Created:  Mon, 15 Apr 2019 13:37:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

class TrigJetMatcher : public edm::stream::EDFilter<> {
   public:
      explicit TrigJetMatcher(const edm::ParameterSet&);
      ~TrigJetMatcher();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Jet>> jetsTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
      std::vector<std::string> trigNames_;
      double dRCut_;
      double jetPtCut_;
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
TrigJetMatcher::TrigJetMatcher(const edm::ParameterSet& iConfig):
    jetsTag_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetsTag"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("triggerObjects")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Jet>>();
   trigNames_ = iConfig.getParameter<std::vector<std::string>>("trigNames");
   dRCut_ = iConfig.getParameter<double>("dRCut");
   jetPtCut_ = iConfig.getParameter<double>("jetPtCut");
}


TrigJetMatcher::~TrigJetMatcher()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrigJetMatcher::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<std::vector<pat::Jet>> jetColl = std::make_unique<std::vector<pat::Jet>>();
   edm::Handle<edm::View<pat::Jet>> pJets;
   edm::Handle<edm::TriggerResults> pTriggerBits;
   edm::Handle<std::vector<pat::TriggerObjectStandAlone> > pTriggerObjects;
   
   iEvent.getByToken(jetsTag_, pJets);
   iEvent.getByToken(triggerBits_, pTriggerBits);
   iEvent.getByToken(triggerObjects_, pTriggerObjects);

   const edm::TriggerNames &names = iEvent.triggerNames(*pTriggerBits);
   std::vector<unsigned int> corrTrigSpot;

   // --- find the user customized trigger path names in the existing trigger path list in TriggerResults of HLT ---
   for (unsigned int i = 0, n = pTriggerBits->size(); i<n; ++i)
   {
       for (std::string iName : trigNames_){
           if (names.triggerName(i).find(iName) != std::string::npos)
           {
               corrTrigSpot.push_back(i);
           } // end if find matched trigger name
       } // end for loop on provided trigger names
   } // end for loop on Trigger Results (HLT)

   // --- prepare for the reco-jet object for trigger matching --- 
   pat::Jet recoLeadingJet;
   pat::TriggerObjectStandAlone trigObj;
   bool checkPassEvent = false;

   for (unsigned int num : corrTrigSpot) // Loop over the trigger paths of TriggerResults cluster that are compatible with customized trigger paths
   {
       const std::string& name = names.triggerName(num);
       bool checkObjMatch = false; // This variable will be used to check if an object(jet) passes at least one customized trigger path

       for (uint iobj = 0; iobj < pTriggerObjects->size(); iobj++)
       {
           trigObj = pTriggerObjects->at(iobj);
           trigObj.unpackPathNames(names); // get the list of trigger paths of each triggered object

           if (trigObj.hasPathName(name, true) && !checkObjMatch) // if there is no yet reco-jet candidates matched with the triggered object
           {
               for(edm::View<pat::Jet>::const_iterator iJet=pJets->begin(); iJet!=pJets->end(); ++iJet) // loop over all the reco-jets
               {
                   double dRCurr = deltaR(*iJet, trigObj); // use dR to match the reco-jet and the triggered object
                   if (iJet->pt() > jetPtCut_ && dRCurr < dRCut_) // match reco-jet with trigger-jet
                   {
                       recoLeadingJet = *iJet;
                       checkObjMatch = true;
                       break;
                   } // end if reco-jet matched with trigger-jet
               } // end loop of reconstructed jets
           } // end if there is a reco-jet candidate matched with the triggered object

           if (checkObjMatch)
           {
               checkPassEvent = true;
               break;
           } // end if checkObjMatch == true && jet trigger
       } // end loop of all the triggered objects

       if (checkPassEvent) break;
   } // end loop of trigger paths in TriggerResults cluster that are compatible with customized trigger paths

   // --- resort reco-jet (the one matched with trigger-jet will be the leading jet) --- 
   if (checkPassEvent)
   {
       jetColl->push_back(recoLeadingJet);

       for(edm::View<pat::Jet>::const_iterator iJet=pJets->begin(); iJet!=pJets->end(); ++iJet)
       {
           if (deltaR(*iJet, recoLeadingJet) < 0.0001 && (iJet->pt() - recoLeadingJet.pt()) < 0.00001) continue; 
           jetColl->push_back(*iJet);
       } // end for loop on reco-jets

       iEvent.put(std::move(jetColl));
   } // end if checkPassEvent == true

   return (checkPassEvent);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TrigJetMatcher::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TrigJetMatcher::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrigJetMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TrigJetMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrigJetMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrigJetMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrigJetMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrigJetMatcher);
