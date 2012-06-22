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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"


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

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

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
      edm::InputTag particleCollection_;
      std::vector<edm::InputTag>  isoValInputTags_;
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
      TH1F * h_cosT_EG;
      TH1F * h_cosT_EG1;
      TH1F * h_cosT_EG2;
      TH1F * h_cosT_EG3;
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
      TH1F * h_sumPtGammaDr13;
      TH1F * h_sumPtGammaDr15;
      TH1F * h_sumPtGammaDr17;
      TH1F * h_sumPtGammaDr20;
      TH1F * h_sumPtGammaDr25;
      TH1F * h_sumPtGammaDr30;
      TH1F * h_sumPtAllGammaDr05;
      TH1F * h_sumPtAllGammaDr06;
      TH1F * h_sumPtAllGammaDr07;
      TH1F * h_sumPtAllGammaDr08;
      TH1F * h_sumPtAllGammaDr09;
      TH1F * h_sumPtAllGammaDr10;
      TH1F * h_sumPtAllGammaDr11;
      TH1F * h_sumPtAllGammaDr13;
      TH1F * h_sumPtAllGammaDr15;
      TH1F * h_sumPtAllGammaDr17;
      TH1F * h_sumPtAllGammaDr20;
      TH1F * h_sumPtAllGammaDr25;
      TH1F * h_sumPtAllGammaDr30;
      
      TH2F * h_angleEnergy;
      TH1F * h_iso;
      TH2F * h_iso_fsr_Dr005;
      TH2F * h_iso_fsr_Dr01;
      TH2F * h_iso_fsr_Dr008;
      TH2F * h_ptRecoPtGenG1_Iso;
      TH2F * h_ptRecoPtGenG2_Iso;
      TH2F * h_ptRecoPtGenG3_Iso;

      TH1F * h_isoV05;
      TH1F * h_isoV08;
      TH1F * h_isoV10;
      TH2F * h_iso_fsrPt;
      TH2F * h_iso_fsrPt_Dr005;
      TH2F * h_iso_fsrPt_Dr01;
      TH2F * h_iso_fsrPt_Dr008;
      TH2F * h_isoIsov_fsr_Dr005;
      TH2F * h_isoIsov_fsr_Dr01;
      TH2F * h_isoIsov_fsr_Dr008;
      TH2F * h_isoIsov_ePt_Dr005;
      TH2F * h_isoIsov_ePt_Dr01;
      TH2F * h_isoIsov_ePt_Dr008;

      TH1F * h_nOfClustersEB;
      TH1F * h_deltaPhiCluEB;
      TH1F * h_deltaEtaCluEB;
      TH1F * h_deltaRCluEB;
      TH2F * h_dPhidEtaCluEB;
      TH1F * h_deltaPhiCluMaxEB;
      TH1F * h_deltaEtaCluMaxEB;
      TH1F * h_deltaRCluMaxEB;
      TH2F * h_dPhidEtaCluMaxEB;

      TH1F * h_nOfClustersEE;
      TH1F * h_deltaXCluEE;
      TH1F * h_deltaYCluEE;
      TH1F * h_deltaPhiCluEE;
      TH1F * h_deltaEtaCluEE;
      TH1F * h_deltaRCluEE;
      TH2F * h_dPhidEtaCluEE;
      TH1F * h_deltaXCluMaxEE;
      TH1F * h_deltaYCluMaxEE;
      TH1F * h_deltaPhiCluMaxEE;
      TH1F * h_deltaEtaCluMaxEE;
      TH1F * h_deltaRCluMaxEE;
      TH2F * h_dPhidEtaCluMaxEE;
      
      TH2F * h_ptCluVsPtRecoEta7;
      TH2F * h_ptCluVsFsrEta7;
      TH2F * h_etaCluVsNCluEta7;
      TH1F * h_isTrackSeed;
      TH1F * h_isGsfRef;

};

fsrValidation::fsrValidation(const edm::ParameterSet& conf)

{
// variable definition =========================
   edgeEB     = 1.479;
   edgeEE     = 3.0;
   edgeTrk    = 2.4;
   edgeEE     = edgeTrk;

  genParticleCollection_ = conf.getUntrackedParameter<edm::InputTag>("genParticleCollection",edm::InputTag("genParticles"));
  goodEPairTag           = conf.getParameter<edm::InputTag>("goodEPair");
  particleCollection_    = conf.getParameter<edm::InputTag>("particleCollection");
  usingPF                = conf.getUntrackedParameter<bool>("usingPF",false);
  usingMC                = conf.getUntrackedParameter<bool>("usingMC",false);
  useSecondChoice        = conf.getUntrackedParameter<bool>("useSecondChoice",false);
  isoValInputTags_       = conf.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

  //now do what ever initialization is needed
  edm::Service<TFileService> fs; 

  //cut from cfg 
  deltaConeGen        = conf.getParameter<double>("deltaRConeGen");
  deltaConeGenLoose   = conf.getParameter<double>("deltaRConeGenloose");
 
  h_cosT_EG = fs->make<TH1F>("h_cosT_EG","h_cosT_EG", 10, -1, 1);
  h_cosT_EG1 = fs->make<TH1F>("h_cosT_EG1","h_cosT_EG1",  10, -1, 1);
  h_cosT_EG2 = fs->make<TH1F>("h_cosT_EG2","h_cosT_EG2", 10, -1, 1);
  h_cosT_EG3 = fs->make<TH1F>("h_cosT_EG3","h_cosT_EG3", 10, -1, 1);
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
  h_angleEnergy = fs->make<TH2F>("h_angleEnergy", "h_angleEnergy", 10, 0.5, 1.5, 200, -2, 200);
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
  h_sumPtGammaDr13 = fs->make<TH1F>("h_sumPtGammaDr13","sumPtGammaDr13",1000,0,50);
  h_sumPtGammaDr15 = fs->make<TH1F>("h_sumPtGammaDr15","sumPtGammaDr15",1000,0,50);
  h_sumPtGammaDr17 = fs->make<TH1F>("h_sumPtGammaDr17","sumPtGammaDr17",1000,0,50);
  h_sumPtGammaDr20 = fs->make<TH1F>("h_sumPtGammaDr20","sumPtGammaDr20",1000,0,50);
  h_sumPtGammaDr25 = fs->make<TH1F>("h_sumPtGammaDr25","sumPtGammaDr25",1000,0,50);
  h_sumPtGammaDr30 = fs->make<TH1F>("h_sumPtGammaDr30","sumPtGammaDr30",1000,0,50);
  h_sumPtAllGammaDr05 = fs->make<TH1F>("h_sumPtAllGammaDr05","sumPtAllGammaDr05",1000,0,50);
  h_sumPtAllGammaDr06 = fs->make<TH1F>("h_sumPtAllGammaDr06","sumPtAllGammaDr06",1000,0,50);
  h_sumPtAllGammaDr07 = fs->make<TH1F>("h_sumPtAllGammaDr07","sumPtAllGammaDr07",1000,0,50);
  h_sumPtAllGammaDr08 = fs->make<TH1F>("h_sumPtAllGammaDr08","sumPtAllGammaDr08",1000,0,50);
  h_sumPtAllGammaDr09 = fs->make<TH1F>("h_sumPtAllGammaDr09","sumPtAllGammaDr09",1000,0,50);
  h_sumPtAllGammaDr10 = fs->make<TH1F>("h_sumPtAllGammaDr10","sumPtAllGammaDr10",1000,0,50);
  h_sumPtAllGammaDr11 = fs->make<TH1F>("h_sumPtAllGammaDr11","sumPtAllGammaDr11",1000,0,50);
  h_sumPtAllGammaDr13 = fs->make<TH1F>("h_sumPtAllGammaDr13","sumPtAllGammaDr13",1000,0,50);
  h_sumPtAllGammaDr15 = fs->make<TH1F>("h_sumPtAllGammaDr15","sumPtAllGammaDr15",1000,0,50);
  h_sumPtAllGammaDr17 = fs->make<TH1F>("h_sumPtAllGammaDr17","sumPtAllGammaDr17",1000,0,50);
  h_sumPtAllGammaDr20 = fs->make<TH1F>("h_sumPtAllGammaDr20","sumPtAllGammaDr20",1000,0,50);
  h_sumPtAllGammaDr25 = fs->make<TH1F>("h_sumPtAllGammaDr25","sumPtAllGammaDr25",1000,0,50);
  h_sumPtAllGammaDr30 = fs->make<TH1F>("h_sumPtAllGammaDr30","sumPtAllGammaDr30",1000,0,50);
  
  h_iso= fs->make<TH1F>("h_iso","h_iso", 100,0,0.3);
  h_iso_fsr_Dr005= fs->make<TH2F>("h_iso_fsr_Dr005","h_iso_fsr_Dr005",100,0,50, 100, 0, 0.3);
  h_iso_fsr_Dr01= fs->make<TH2F>("h_iso_fsr_Dr01","h_iso_fsr_Dr01",100,0, 50, 100, 0, 0.3);
  h_iso_fsr_Dr008= fs->make<TH2F>("h_iso_fsr_Dr008","h_iso_fsr_Dr008",100,0,50, 100, 0, 0.3);
  h_ptRecoPtGenG1_Iso = fs->make<TH2F>("h_ptRecoPtGenG1_Iso","h_ptRecoPtGenG1_Iso",200,-10,10, 100, 0, 0.3);
  h_ptRecoPtGenG2_Iso = fs->make<TH2F>("h_ptRecoPtGenG2_Iso","h_ptRecoPtGenG2_Iso",200,-10,10, 100, 0, 0.3);
  h_ptRecoPtGenG3_Iso = fs->make<TH2F>("h_ptRecoPtGenG3_Iso","h_ptRecoPtGenG3_Iso",200,-10,10, 100, 0, 0.3);

  h_isoV05= fs->make<TH1F>("h_isoV05","h_isoV05", 100,0,0.3);
  h_isoV08= fs->make<TH1F>("h_isoV08","h_isoV08", 100,0,0.3);
  h_isoV10= fs->make<TH1F>("h_isoV10","h_isoV10", 100,0,0.3);
  h_iso_fsrPt= fs->make<TH2F>("h_iso_fsrPt","h_iso_fsrPt",100,0, 1, 100, 0, 0.3);
  h_iso_fsrPt_Dr005= fs->make<TH2F>("h_iso_fsrPt_Dr005","h_iso_fsrPt_Dr005",100,0, 1, 100, 0, 0.3);
  h_iso_fsrPt_Dr01= fs->make<TH2F>("h_iso_fsrPt_Dr01","h_iso_fsrPt_Dr01",100,0, 1, 100, 0, 0.3);
  h_iso_fsrPt_Dr008= fs->make<TH2F>("h_iso_fsrPt_Dr008","h_iso_fsrPt_Dr008",100,0,1, 100, 0, 0.3);
  h_isoIsov_fsr_Dr005= fs->make<TH2F>("h_isoIsov_fsr_Dr005","h_isoIsov_fsrPt_Dr005",100,0, 1, 200, 0, 0.1);
  h_isoIsov_fsr_Dr01= fs->make<TH2F>("h_isoIsov_fsr_Dr01","h_isoIsov_fsrPt_Dr01",100,0, 1, 200, 0, 0.1);
  h_isoIsov_fsr_Dr008= fs->make<TH2F>("h_isoIsov_fsr_Dr008","h_isoIsov_fsrPt_Dr008",100,0,1, 200, 0, 0.1);
  h_isoIsov_ePt_Dr005= fs->make<TH2F>("h_isoIsov_ePt_Dr005","h_isoIsov_ePt_Dr005",150,0, 150, 200, 0, 0.1);
  h_isoIsov_ePt_Dr01= fs->make<TH2F>("h_isoIsov_ePt_Dr01","h_isoIsov_ePt_Dr01",150,0, 150, 200, 0, 0.1);
  h_isoIsov_ePt_Dr008= fs->make<TH2F>("h_isoIsov_ePt_Dr008","h_isoIsov_ePt_Dr008",150,0,150, 200, 0, 0.1);
  
  h_nOfClustersEB = fs->make<TH1F>("h_nOfClustersEB","nOfClustersEB", 20,0,20);
  h_deltaPhiCluEB = fs->make<TH1F>("h_deltaPhiCluEB","deltaPhiCluEB", 100,0,0.4);
  h_deltaEtaCluEB = fs->make<TH1F>("h_deltaEtaCluEB","deltaEtaCluEB", 200,0,0.2);
  h_deltaRCluEB   = fs->make<TH1F>("h_deltaRCluEB","deltaRCluEB", 100,0,0.4);
  h_dPhidEtaCluEB = fs->make<TH2F>("h_dPhidEtaCluEB","dPhidEtaCluEB",100,0,0.4,100,0,0.1);
  h_deltaPhiCluMaxEB = fs->make<TH1F>("h_deltaPhiCluMaxEB","deltaPhiCluMaxEB", 100,0,0.4);
  h_deltaEtaCluMaxEB = fs->make<TH1F>("h_deltaEtaCluMaxEB","deltaEtaCluMaxEB", 200,0,0.2);
  h_deltaRCluMaxEB   = fs->make<TH1F>("h_deltaRCluMaxEB","deltaRCluMaxEB", 100,0,0.4);
  h_dPhidEtaCluMaxEB = fs->make<TH2F>("h_dPhidEtaCluMaxEB","dPhidEtaCluMaxEB",100,0,0.4,100,0,0.1);
  
  h_nOfClustersEE = fs->make<TH1F>("h_nOfClustersEE","nOfClustersEE", 20,0,20);
  h_deltaXCluEE = fs->make<TH1F>("h_deltaXCluEE","deltaXCluEE", 200,0,100);
  h_deltaYCluEE = fs->make<TH1F>("h_deltaYCluEE","deltaYCluEE", 200,0,100);
  h_deltaPhiCluEE = fs->make<TH1F>("h_deltaPhiCluEE","deltaPhiCluEE", 100,0,0.4);
  h_deltaEtaCluEE = fs->make<TH1F>("h_deltaEtaCluEE","deltaEtaCluEE", 200,0,0.2);
  h_deltaRCluEE   = fs->make<TH1F>("h_deltaRCluEE","deltaRCluEE", 100,0,0.4);
  h_dPhidEtaCluEE = fs->make<TH2F>("h_dPhidEtaCluEE","dPhidEtaCluEE",100,0,0.4,100,0,0.1);
  h_deltaXCluMaxEE = fs->make<TH1F>("h_deltaXCluMaxEE","deltaXCluMaxEE", 200,0,100);
  h_deltaYCluMaxEE = fs->make<TH1F>("h_deltaYCluMaxEE","deltaYCluMaxEE", 200,0,100);
  h_deltaPhiCluMaxEE = fs->make<TH1F>("h_deltaPhiCluMaxEE","deltaPhiCluMaxEE", 100,0,0.4);
  h_deltaEtaCluMaxEE = fs->make<TH1F>("h_deltaEtaCluMaxEE","deltaEtaCluMaxEE", 200,0,0.2);
  h_deltaRCluMaxEE   = fs->make<TH1F>("h_deltaRCluMaxEE","deltaRCluMaxEE", 100,0,0.4);
  h_dPhidEtaCluMaxEE = fs->make<TH2F>("h_dPhidEtaCluMaxEE","dPhidEtaCluMaxEE",100,0,0.4,100,0,0.1);

  h_ptCluVsPtRecoEta7 = fs->make<TH2F>("h_ptCluVsPtRecoEta7","ptCluVsPtReco", 200,0,10, 200,0,50);
  h_ptCluVsFsrEta7 = fs->make<TH2F>("h_ptCluVsFsrEta7","ptCluVsFsr", 200,0,10, 100,0,2);
  h_etaCluVsNCluEta7 = fs->make<TH2F>("h_etaCluVsNCluEta7","etaCluVsNClu", 200,-2.5,2.5, 15, 0, 15);
  h_isGsfRef = fs->make<TH1F>("h_isGsfRef","isGsfRef", 2, 0,2);
  h_isTrackSeed = fs->make<TH1F>("h_isTrackSeed","isTrackSeed", 6, 0,6);
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
