#ifndef FSRAnalyzer_h_
#define FSRAnalyzer_h_


// system include files
#include <memory>
//#include <typeinfo> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"



// --- Order GenParticles by Pt: 
bool orderPartByPt(const HepMC::GenParticle* p1, const HepMC::GenParticle* p2){
  return (p1->momentum().perp() > p2->momentum().perp());
}


//
// class declaration
//

class FSRAnalyzer : public edm::EDAnalyzer {
 public:
  explicit FSRAnalyzer(const edm::ParameterSet&);
  ~FSRAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);



 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  TVector3 extrapolateToECAL(TVector3 const &momentum, TVector3 const &position, int const &subdet);
  bool is_in_ellipse(HepMC::GenParticle const &ph, HepMC::GenParticle const &ele, double const &deltaR);



  // ----------member data ---------------------------

  edm::InputTag goodEPairTag;
  edm::InputTag genParticleCollection_;
  edm::InputTag hepMCProduct_;

  bool isMC;

  enum SubDetector { BARREL, ENDCAP };

  // --- histograms
  
  TH1F * h_nele_noFSR;
  TH1F * h_nele_FSR;

  TH1F * h_Mee_noFSR;
  TH1F * h_pT_Z_noFSR;
  TH1F * h_Y_Z_noFSR;
  TH1F * h_pT_e1_noFSR;
  TH1F * h_pT_e2_noFSR;

  TH1F * h_Mee_FSR;
  TH1F * h_pT_Z_FSR;
  TH1F * h_Y_Z_FSR;
  TH1F * h_pT_e1_FSR;
  TH1F * h_pT_e2_FSR;

  TH1F * h_ph_mom;
  TH1F * h_ph_grandma;

  TH1F * h_FSR_n;
  TH1F * h_FSR_et;
  TH1F * h_FSR_sumet;
  TH1F * h_FSR_e;
  TH1F * h_FSR_dr;
  TH1F * h_FSR_angle;
  TH1F * h_FSR_dr_w;
  TH1F * h_FSR_angle_w;
  TH2F * h_FSR_deta_dphi;
  TH2F * h_FSR_deta_dphi_w;
  TH2F * h_FSR_dr_et;
  TH2F * h_FSR_dr_e;

  TH1F * h_FSR_frac[100];  	
  TProfile * h_FSR_frac_p;  	
  TH1F * h_FSR_totet[100];  	
  TProfile * h_FSR_totet_p;  	
 
  TH1F * h_pt_ph[100];
  TH1F * h_pt_ph_fsr[100];  
  TH1F * h_pt_ph_e[100];    
  TH1F * h_pt_ph_pi0[100];  
  TH1F * h_pt_ph_other[100];
  TH1F * h_pt_ph_PU[100];

  TH2F * h_dr_et_ph[100];
  TH2F * h_dr_et_ph_fsr[100];
  TH2F * h_dr_et_ph_e[100];
  TH2F * h_dr_et_ph_pi0[100];
  TH2F * h_dr_et_ph_other[100];
  TH2F * h_dr_et_ph_PU[100];

  TH2F * h_dr_e_ph[100];
  TH2F * h_dr_e_ph_fsr[100];
  TH2F * h_dr_e_ph_e[100];
  TH2F * h_dr_e_ph_pi0[100];
  TH2F * h_dr_e_ph_other[100];

  TH1F * h_dz_ph[100];
  TH1F * h_dz_ph_fsr[100];  
  TH1F * h_dz_ph_other[100];  
  TH1F * h_dz_ph_PU[100];

  TH1F * h_id_ph_other[100];

  TH1F * h_frac_fsr[100];
  TH1F * h_frac_other[100];

  TProfile * h_frac_fsr_p;
  TProfile * h_frac_other_p;
  TProfile * h_frac_PU_p;

  TH1F * h_totet_fsr[100];
  TH1F * h_totet_other[100];

  TProfile * h_totet_fsr_p;
  TProfile * h_totet_other_p;
  TProfile * h_totet_PU_p;


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
FSRAnalyzer::FSRAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

  goodEPairTag           = iConfig.getParameter<edm::InputTag>("goodEPair");
  genParticleCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("genParticleCollection",edm::InputTag("genParticles"));
  hepMCProduct_          = iConfig.getUntrackedParameter<edm::InputTag>("generator",edm::InputTag("generator"));


  isMC = iConfig.getUntrackedParameter<bool>("usingMC",false);


  edm::Service<TFileService> fs; 

  // --- histograms
  
  h_nele_noFSR  = fs->make<TH1F>("h_nele_noFSR","# of electrons (w/o FSR); # of electrons",10,0.,10.);
  h_nele_FSR    = fs->make<TH1F>("h_nele_FSR","# of electrons (w/ FSR); # of electrons",10,0.,10.);

  h_Mee_noFSR   = fs->make<TH1F>("h_Mee_noFSR","ee invariant mass (w/o FSR); M_{ee} [GeV/c^{2}]",200,0.,200.);
  h_pT_Z_noFSR  = fs->make<TH1F>("h_pT_Z_noFSR","Z P_{T} (w/o FSR); P_{T} [GeV/c]",200,0.,200.);
  h_Y_Z_noFSR   = fs->make<TH1F>("h_Y_Z_noFSR","Z rapidity (w/o FSR); Rapidity",200,-2.5,2.5);
  h_pT_e1_noFSR = fs->make<TH1F>("h_pT_e1_noFSR","Leading electron P_{T} (w/o FSR); P_{T} [GeV/c]",200,0.,200.);
  h_pT_e2_noFSR = fs->make<TH1F>("h_pT_e2_noFSR","Trailing electron P_{T} (w/o FSR); P_{T} [GeV/c]",200,0.,200.);

  h_Mee_FSR     = fs->make<TH1F>("h_Mee_FSR","ee invariant mass (w/  FSR); M_{ee} [GeV/c^{2}]",200,0.,200.);
  h_pT_Z_FSR    = fs->make<TH1F>("h_pT_Z_FSR","Z P_{T} (w/  FSR); P_{T} [GeV/c]",200,0.,200.);
  h_Y_Z_FSR     = fs->make<TH1F>("h_Y_Z_FSR","Z rapidity (w/  FSR); Rapidity",200,-2.5,2.5);
  h_pT_e1_FSR   = fs->make<TH1F>("h_pT_e1_FSR","Leading electron P_{T} (w/  FSR); P_{T} [GeV/c]",200,0.,200.);
  h_pT_e2_FSR   = fs->make<TH1F>("h_pT_e2_FSR","Trailing electron P_{T} (w/  FSR); P_{T} [GeV/c]",200,0.,200.);

  h_ph_mom      = fs->make<TH1F>("h_ph_mom",  "photon mother; PDG code", 6000,0.,6000.);
  h_ph_grandma  = fs->make<TH1F>("h_ph_grandma", "photon grandmother; PDG code", 6000,0.,6000.);

  h_FSR_n       = fs->make<TH1F>("h_FSR_n", "# of FSR photons", 10,0.,10.);
  h_FSR_et      = fs->make<TH1F>("h_FSR_et", "FSR photons E_{T};E_{T} [GeV]",250,0.,50.);
  h_FSR_sumet   = fs->make<TH1F>("h_FSR_sumet", "FSR photons #Sigma E_{T};E_{T} [GeV]",250,0.,50.);
  h_FSR_e       = fs->make<TH1F>("h_FSR_e", "FSR photons E;E [GeV]",250,0.,50.);
  h_FSR_dr      = fs->make<TH1F>("h_FSR_dr", "#DeltaR between the FSR photon and electron;#DeltaR",600,0.,6.);
  h_FSR_angle   = fs->make<TH1F>("h_FSR_angle", "Angle between the FSR photon and electron;#alpha [rad]",200,0.,TMath::Pi());
  h_FSR_dr_w    = fs->make<TH1F>("h_FSR_dr_w", "#DeltaR between the FSR photon and electron;#DeltaR",600,0.,6.);
  h_FSR_angle_w = fs->make<TH1F>("h_FSR_angle_w", "Angle between the FSR photon and electron;#alpha [rad]",200,0.,TMath::Pi());
  h_FSR_deta_dphi = fs->make<TH2F>("h_FSR_deta_dphi", "#Delta#eta vs #Delta#phi;#Delta#phi;#Delta#eta",1000,-TMath::Pi(),TMath::Pi(),1000,-3.,3.);
  h_FSR_deta_dphi_w = fs->make<TH2F>("h_FSR_deta_dphi_w", "#Delta#eta vs #Delta#phi;#Delta#phi;#Delta#eta",1000,-TMath::Pi(),TMath::Pi(),1000,-3.,3.);
  h_FSR_dr_et   = fs->make<TH2F>("h_FSR_dr_et", "#DeltaR vs E_{T};E_{T} [GeV];#DeltaR",1000,0.,10.,1000,0.,5.);
  h_FSR_dr_e    = fs->make<TH2F>("h_FSR_dr_e ", "#DeltaR vs E;E [GeV];#DeltaR",1000,0.,10.,1000,0.,5.);


  for (int ih=0; ih<100; ++ih) {

    double dR_cut  = 0.01 * (1. + ih);

    
    // photon Et spectra
    TString hname  = Form("h_pt_ph_%d",ih);
    TString htitle = Form("photons in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph[ih]    = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_pt_ph_fsr_%d",ih);
    htitle = Form("FSR photons in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph_fsr[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_pt_ph_e_%d",ih);
    htitle = Form("photons from other e's in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph_e[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_pt_ph_pi0_%d",ih);
    htitle = Form("photons from #pi^{0}'s in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph_pi0[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_pt_ph_other_%d",ih);
    htitle = Form("photons from other particles in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph_other[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_pt_ph_PU_%d",ih);
    htitle = Form("photons from PU in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_pt_ph_PU[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 100,0.,20.);

    hname  = Form("h_id_ph_other_%d",ih);
    htitle = Form("photon mother in #DeltaR<%4.2f; PDG code",dR_cut);
    h_id_ph_other[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(), 6000,0.,6000.);


    // photon DeltaR vs Et
    hname  = Form("h_dr_et_ph_%d",ih);
    htitle = Form("photons in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_et_ph[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_et_ph_fsr_%d",ih);
    htitle = Form("FSR photons in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_et_ph_fsr[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_et_ph_e%d",ih);
    htitle = Form("photons from e's in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_et_ph_e[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_et_ph_pi0_%d",ih);
    htitle = Form("photons from #pi^{0}'s in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_et_ph_pi0[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_et_ph_other_%d",ih);
    htitle = Form("photons from other particles in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_et_ph_other[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);


    // photon DeltaR vs E
    hname  = Form("h_dr_e_ph_%d",ih);
    htitle = Form("photons in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_e_ph[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_e_ph_fsr_%d",ih);
    htitle = Form("FSR photons in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_e_ph_fsr[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_e_ph_e%d",ih);
    htitle = Form("photons from e's in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_e_ph_e[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_e_ph_pi0_%d",ih);
    htitle = Form("photons from #pi^{0}'s in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_e_ph_pi0[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);

    hname  = Form("h_dr_e_ph_other_%d",ih);
    htitle = Form("photons from other particles in #DeltaR<%4.2f;E_{T} [GeV];#DeltaR",dR_cut);
    h_dr_e_ph_other[ih] = fs->make<TH2F>(hname.Data(), htitle.Data(),1000,0.,10.,1000,0.,1.);


    // Delta Z
    hname  = Form("h_dz_ph_%d",ih);
    htitle = Form("z_{#gamma}-z_{e} in #DeltaR<%4.2f;#Deltaz [mm]",dR_cut);
    h_dz_ph[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),4000,-200.,200.);

    hname  = Form("h_dz_ph_fsr_%d",ih);
    htitle = Form("z_{#gamma}-z_{e} in #DeltaR<%4.2f (FSR);#Deltaz [mm]",dR_cut);
    h_dz_ph_fsr[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),4000,-200.,200.);

    hname  = Form("h_dz_ph_other_%d",ih);
    htitle = Form("z_{#gamma}-z_{e} in #DeltaR<%4.2f (other #gamma);#Deltaz [mm]",dR_cut);
    h_dz_ph_other[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),4000,-200.,200.);

    hname  = Form("h_dz_ph_PU_%d",ih);
    htitle = Form("z_{#gamma}-z_{e} in #DeltaR<%4.2f (from PU);#Deltaz [mm]",dR_cut);
    h_dz_ph_PU[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),4000,-200.,200.);


    // Et fraction
    hname  = Form("h_FSR_frac_%d",ih);
    htitle = Form("fraction of FSR energy in #DeltaR<%4.2f;fraction",dR_cut);
    h_FSR_frac[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),220,0.,1.1);
    h_FSR_frac[ih]->Sumw2();

    hname  = Form("h_frac_fsr_%d",ih);
    htitle = Form("fraction of FSR #gamma energy in #DeltaR<%4.2f;fraction",dR_cut);
    h_frac_fsr[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),220,0.,1.1);
    h_frac_fsr[ih]->Sumw2();

    hname  = Form("h_frac_other_%d",ih);
    htitle = Form("fraction of other #gamma's energy in #DeltaR<%4.2f;fraction",dR_cut);
    h_frac_other[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),220,0.,1.1);
    h_frac_other[ih]->Sumw2();

    // Et tot
    hname  = Form("h_FSR_totet_%d",ih);
    htitle = Form("tot E_{T} of FSR energy in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_FSR_totet[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),100,0.,50.);
    h_FSR_totet[ih]->Sumw2();

    hname  = Form("h_totet_fsr_%d",ih);
    htitle = Form("tot E_{T} of FSR #gamma energy in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_totet_fsr[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),100,0.,50.);
    h_totet_fsr[ih]->Sumw2();

    hname  = Form("h_totet_other_%d",ih);
    htitle = Form("tot E_{T} of other #gamma's energy in #DeltaR<%4.2f;E_{T} [GeV]",dR_cut);
    h_totet_other[ih] = fs->make<TH1F>(hname.Data(), htitle.Data(),100,0.,50.);
    h_totet_other[ih]->Sumw2();

  }

    
  h_FSR_frac_p = fs->make<TProfile>("h_FSR_frac_p", "fraction of FSR energy in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_FSR_frac_p->Sumw2();

  h_frac_fsr_p = fs->make<TProfile>("h_frac_fsr_p","fraction of FSR #gamma energy in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_frac_fsr_p->Sumw2();

  h_frac_other_p = fs->make<TProfile>("h_frac_other_p","fraction of other #gamma energy in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_frac_other_p->Sumw2();
   
  h_frac_PU_p = fs->make<TProfile>("h_frac_PU_p","fraction of PU #gamma energy in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_frac_PU_p->Sumw2();

  h_FSR_totet_p = fs->make<TProfile>("h_FSR_totet_p", "tot E_{T} of FSR energy in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_FSR_totet_p->Sumw2();

  h_totet_fsr_p = fs->make<TProfile>("h_totet_fsr_p","tot E_{T} of FSR #gamma's in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_totet_fsr_p->Sumw2();

  h_totet_other_p = fs->make<TProfile>("h_totet_other_p","tot E_{T} of other #gamma's in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_totet_other_p->Sumw2();
   
  h_totet_PU_p = fs->make<TProfile>("h_totet_PU_p","tot E_{T} of PU #gamma's in a #DeltaR cone;#DeltaR", 100,0.,1.,0.,1.1);
  h_totet_PU_p->Sumw2();
   
  
}


FSRAnalyzer::~FSRAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


#endif

