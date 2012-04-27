#ifndef fsrValidation_h
#define fsrValidation_h

// system files
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>
//#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// root includes
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace edm;
using namespace reco;
using namespace std;

bool debugJet=false; //If true it will activate the cout verbosity


/////
// class declaration
/////

class fsrValidation : public edm::EDAnalyzer {
   public:
      explicit fsrValidation(const edm::ParameterSet&);
      ~fsrValidation();

   private:
      virtual void beginJob() ;
      virtual void beginRun(edm::Run const &, const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      double distR(TLorentzVector ,TLorentzVector);
      template <class T> bool searchHistory(T const&);
      //bool fsrValidation::searchHistory(const reco::GenParticle* const &);

      // ----------member data ---------------------------

      //Retrieved from the .py
      edm::InputTag genParticleCollection_;
      edm::InputTag goodEPairTag;
      bool usingMC;
      bool usingPF;
      bool useSecondChoice;

      
     
// constant variable definition ========== 
      double edgeEB;
      double edgeEE;
      double edgeTrk;
      double deltaConeGen;
      double deltaConeGenLoose;
      TLorentzVector e1_test, e2_test;
      std::vector<TLorentzVector> gamma1_gen, gamma2_gen;

      // histo initialization
      TH1F * h_ptReco;
      TH1F * h_ptGen;
      TH1F * h_ptGen3;
      TH1F * h_ptGen13;
      TH1F * h_sumPtGamma;
      TH1F * h_ptRecoG1;
      TH1F * h_ptRecoG2;
      TH1F * h_ptRecoG3;
      TH1F * h_ptRecoG4;
      TH1F * h_ptRecoPtGenG1;
      TH1F * h_ptRecoPtGenG2;
      TH1F * h_ptRecoPtGenG3;
      TH1F * h_ptRecoPtGenG4;
      TH2F * h_ptRecoPtGenVsPtGenG1;
      TH2F * h_ptRecoPtGenVsPtGenG2;
      TH2F * h_ptRecoPtGenVsPtGenG3;
      TH2F * h_ptRecoPtGenVsPtGenG4;
      TH1F * h_sumPtGammaDr05;
      TH1F * h_sumPtGammaDr06;
      TH1F * h_sumPtGammaDr07;
      TH1F * h_sumPtGammaDr08;
      TH1F * h_sumPtGammaDr09;
      TH1F * h_sumPtGammaDr10;
      TH1F * h_sumPtGammaDr11;
     
};

fsrValidation::fsrValidation(const edm::ParameterSet& conf)

{
// variable definition =========================
   edgeEB     = 1.479;
   edgeEE     = 3.0;
   edgeTrk    = 2.4;
   edgeEE     = edgeTrk;

  genParticleCollection_ = conf.getUntrackedParameter<edm::InputTag>("genParticleCollection",edm::InputTag("genParticles"));
  goodEPairTag        = conf.getParameter<edm::InputTag>("goodEPair");
  usingPF             = conf.getUntrackedParameter<bool>("usingPF",false);
  usingMC             = conf.getUntrackedParameter<bool>("usingMC",false);
  useSecondChoice     = conf.getUntrackedParameter<bool>("useSecondChoice",false);

  //now do what ever initialization is needed
  edm::Service<TFileService> fs; 

  //cut from cfg 
  deltaConeGen        = conf.getParameter<double>("deltaRConeGen");
  deltaConeGenLoose   = conf.getParameter<double>("deltaRConeGenloose");
  
  h_ptReco = fs->make<TH1F>("h_ptReco","ptReco",200,0,200);
  h_ptGen = fs->make<TH1F>("h_ptGen","ptGen",200,0,200);
  h_ptGen3 = fs->make<TH1F>("h_ptGen3","ptGen3",200,0,200);
  h_ptGen13 = fs->make<TH1F>("h_ptGen13","ptGen1diff3",1000,0,50);
  h_sumPtGamma = fs->make<TH1F>("h_sumPtGamma","sumPtGamma",1000,0,50);
  h_ptRecoG1 = fs->make<TH1F>("h_ptRecoG1","ptRecoG1",200,0,200);
  h_ptRecoG2 = fs->make<TH1F>("h_ptRecoG2","ptRecoG2",200,0,200);
  h_ptRecoG3 = fs->make<TH1F>("h_ptRecoG3","ptRecoG3",200,0,200);
  h_ptRecoG4 = fs->make<TH1F>("h_ptRecoG4","ptRecoG4",200,0,200);
  h_ptRecoPtGenG1 = fs->make<TH1F>("h_ptRecoPtGenG1","ptRecoPtGenG1",200,-10,10);
  h_ptRecoPtGenG2 = fs->make<TH1F>("h_ptRecoPtGenG2","ptRecoPtGenG2",200,-10,10);
  h_ptRecoPtGenG3 = fs->make<TH1F>("h_ptRecoPtGenG3","ptRecoPtGenG3",200,-10,10);
  h_ptRecoPtGenG4 = fs->make<TH1F>("h_ptRecoPtGenG4","ptRecoPtGenG4",200,-10,10);
  h_ptRecoPtGenVsPtGenG1 = fs->make<TH2F>("h_ptRecoPtGenVsPtGenG1","ptRecoPtGenVsPtGenG1",200,0,200,200,-10,10);
  h_ptRecoPtGenVsPtGenG2 = fs->make<TH2F>("h_ptRecoPtGenVsPtGenG2","ptRecoPtGenVsPtGenG2",200,0,200,200,-10,10);
  h_ptRecoPtGenVsPtGenG3 = fs->make<TH2F>("h_ptRecoPtGenVsPtGenG3","ptRecoPtGenVsPtGenG3",200,0,200,200,-10,10);
  h_ptRecoPtGenVsPtGenG4 = fs->make<TH2F>("h_ptRecoPtGenVsPtGenG4","ptRecoPtGenVsPtGenG4",200,0,200,200,-10,10);
  h_sumPtGammaDr05 = fs->make<TH1F>("h_sumPtGammaDr05","sumPtGammaDr05",1000,0,50);
  h_sumPtGammaDr06 = fs->make<TH1F>("h_sumPtGammaDr06","sumPtGammaDr06",1000,0,50);
  h_sumPtGammaDr07 = fs->make<TH1F>("h_sumPtGammaDr07","sumPtGammaDr07",1000,0,50);
  h_sumPtGammaDr08 = fs->make<TH1F>("h_sumPtGammaDr08","sumPtGammaDr08",1000,0,50);
  h_sumPtGammaDr09 = fs->make<TH1F>("h_sumPtGammaDr09","sumPtGammaDr09",1000,0,50);
  h_sumPtGammaDr10 = fs->make<TH1F>("h_sumPtGammaDr10","sumPtGammaDr10",1000,0,50);
  h_sumPtGammaDr11 = fs->make<TH1F>("h_sumPtGammaDr11","sumPtGammaDr11",1000,0,50);

}


fsrValidation::~fsrValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


template <class T>
bool fsrValidation::searchHistory(T const& itSearch)
//bool fsrValidation::searchHistory(const reco::GenParticle* const& itSearch)
{
   bool pass = false;
   if (fabs(itSearch->pdgId()) ==11 && itSearch->status()==3 && itSearch->mother()->pdgId()==23) {
      e1_test.SetPtEtaPhiM(itSearch->pt(),itSearch->eta(),
                           itSearch->phi(),itSearch->mass());

      pass=true;
      return pass;
   } else if (itSearch->status()==3 ){
      pass=false;
      return pass;
   } else {
      pass=searchHistory(itSearch->mother());
   }
   return pass;
}

#endif
