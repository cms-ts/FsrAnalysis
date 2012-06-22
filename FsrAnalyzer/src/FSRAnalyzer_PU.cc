#include "Fsr/FsrAnalyzer/interface/FSRAnalyzer_PU.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// member functions
//

// ------------ method called for each event  ------------
void 
FSRAnalyzer_PU::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  if (isMC) {

    edm::Handle<CrossingFrame<edm::HepMCProduct> > cf_hepmc;
    bool gotHepMC = iEvent.getByLabel("mix","generator",cf_hepmc);
    if (!gotHepMC) std::cout<<" Could not read HepMCProducts!!!!"<<std::endl;
    else {


      // --- Global variables:

      const unsigned int n_dr_bins = 100;

      vector <HepMC::GenParticle*> electrons_beforeFSR;
      vector <HepMC::GenParticle*> electrons_afterFSR;

      TLorentzVector p_ele0_noFSR(0.,0.,0.,0.), 
	             p_ele1_noFSR(0.,0.,0.,0.), 
	             p_ele0_FSR(0.,0.,0.,0.),
	             p_ele1_FSR(0.,0.,0.,0.);

      double sumEt_tot_fsr0 =  0.;
      double sumEt_tot_fsr1 =  0.;
      double sumEt_fsr_ele0[n_dr_bins]   = {0.};
      double sumEt_other_ele0[n_dr_bins] = {0.};
      double sumEt_PU_ele0[n_dr_bins]    = {0.};
      double sumEt_fsr_ele1[n_dr_bins]   = {0.};
      double sumEt_other_ele1[n_dr_bins] = {0.};
      double sumEt_PU_ele1[n_dr_bins]    = {0.};

      map <HepMC::GenParticle*, vector <HepMC::GenParticle*> > FSR_photons;


      // --- Access the MixingCollection of the HepMCProducts and find the Z vertex:

      MixCollection<HepMCProduct>* cfhepmcprod = new MixCollection<HepMCProduct>(cf_hepmc.product());

      int Zvtx_index = -999;
      vector<unsigned int> PUvtx_index;
      
      unsigned int ihepmc = 0;
      for (MixCollection<edm::HepMCProduct>::iterator cfihepmc=cfhepmcprod->begin(); cfihepmc!=cfhepmcprod->end(); ++cfihepmc){

	if ( cfihepmc.getTrigger() == 1 )
	  Zvtx_index = ihepmc;
	else
	  PUvtx_index.push_back(ihepmc);

	ihepmc++;

      }




      // =====================================================================================================
      //  Z vertex


      if ( Zvtx_index == -999 ) {
 	cout << " FSRAnalyzer_PU::analyze: no Z vertex found " << endl;
	return;
      }
      else{

	const HepMC::GenEvent *genEvent = cfhepmcprod->getObject(Zvtx_index).GetEvent();


	// --- Find the Z daughters (before and after the FSR):

	for ( HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin(); p != genEvent->particles_end(); ++p ) {
	  
	  if ( fabs((*p)->pdg_id()) !=11 ) continue;

	  TLorentzVector p_gen( (*p)->momentum().px(),
				(*p)->momentum().py(),
				(*p)->momentum().pz(),
				(*p)->momentum().e() );

	  HepMC::GenParticle* mother = (*((*p)->production_vertex()->particles_in_const_begin()));


	  // --- electrons before FSR:
	  if ( (*p)->status() == 3 && fabs(mother->pdg_id()) == 23 ) {
	
	    //if ( p_gen.Pt()>20. && fabs(p_gen.Eta())<2.4 )
	    electrons_beforeFSR.push_back(*p);
	
	  }

	  // --- electrons after FSR:
	  if ( (*p)->status() == 1 ) {
	
	    if ( p_gen.Pt()>20. && fabs(p_gen.Eta())<2.4 )
	      electrons_afterFSR.push_back(*p);

	  }


	} // for p


	// --- Event selection:

	// Keep only Z-->ee events
	if ( electrons_beforeFSR.size() != 2 ) return;
	if ( electrons_beforeFSR[0]->pdg_id()*electrons_beforeFSR[1]->pdg_id() > 0 ) return;

	h_nele_noFSR->Fill(electrons_beforeFSR.size());
	h_nele_FSR->Fill(electrons_afterFSR.size());
    
 	if ( electrons_afterFSR.size() < 2 ) return;

	// Order the electrons by Pt
	std::sort (electrons_afterFSR.begin(), electrons_afterFSR.end(), orderPartByPt);
	// This is to keep the correct correspondence between afterFSR and beforeFSR
	if ( electrons_afterFSR[0]->pdg_id() != electrons_beforeFSR[0]->pdg_id() ) 
	  std::swap(electrons_beforeFSR[0],electrons_beforeFSR[1]);

	// Check for opposite charges
	if ( electrons_afterFSR[0]->pdg_id()*electrons_afterFSR[1]->pdg_id() > 0 ) return;
    

	p_ele0_noFSR.SetPxPyPzE( electrons_beforeFSR[0]->momentum().px(),
				 electrons_beforeFSR[0]->momentum().py(),
				 electrons_beforeFSR[0]->momentum().pz(),
				 electrons_beforeFSR[0]->momentum().e() );
	p_ele1_noFSR.SetPxPyPzE( electrons_beforeFSR[1]->momentum().px(),
				 electrons_beforeFSR[1]->momentum().py(),
				 electrons_beforeFSR[1]->momentum().pz(),
				 electrons_beforeFSR[1]->momentum().e() );

	p_ele0_FSR.SetPxPyPzE( electrons_afterFSR[0]->momentum().px(),
			       electrons_afterFSR[0]->momentum().py(),
			       electrons_afterFSR[0]->momentum().pz(),
			       electrons_afterFSR[0]->momentum().e() );
	p_ele1_FSR.SetPxPyPzE( electrons_afterFSR[1]->momentum().px(),
			       electrons_afterFSR[1]->momentum().py(),
			       electrons_afterFSR[1]->momentum().pz(),
			       electrons_afterFSR[1]->momentum().e() );



	// Apply the Z mass window cut
	const double m_ll_FSR = (p_ele0_FSR+p_ele1_FSR).M();
	if ( m_ll_FSR<71. || m_ll_FSR>111. ) return;



	// --- Loop over the final state photons:
    
	for ( HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin(); p != genEvent->particles_end(); ++p ) {

	  if ( (*p)->status() != 1 || fabs((*p)->pdg_id()) != 22 ) continue;
      

	  TLorentzVector p_ph( (*p)->momentum().px(),
			       (*p)->momentum().py(),
			       (*p)->momentum().pz(),
			       (*p)->momentum().e() );

      
	  // Get the photon mother and grandmother
	  HepMC::GenParticle* mother      = (*((*p)->production_vertex()->particles_in_const_begin()));
	  HepMC::GenParticle* grandmother = NULL;
	  if (mother && mother->production_vertex()) 
	    grandmother = (*(mother->production_vertex()->particles_in_const_begin()));

      
	  if ( fabs(mother->pdg_id()) == 11 ) {

	    bool from_ele0 = false;
	    bool from_ele1 = false;

	    while ( grandmother->pdg_id()==mother->pdg_id() ){
	  
	      from_ele0 = ( grandmother == electrons_beforeFSR[0] ); 
	      from_ele1 = ( grandmother == electrons_beforeFSR[1] ); 

	      grandmother = (*(grandmother->production_vertex()->particles_in_const_begin()));
	      if (!grandmother) break;

	    }


	    if ( from_ele0 && fabs(grandmother->pdg_id())==23 ){
	  
	      sumEt_tot_fsr0 += p_ph.Pt();
	  
	      h_FSR_et->Fill(p_ph.Pt());
	      h_FSR_e->Fill(p_ph.E());
	      h_FSR_dr->Fill(p_ph.DeltaR(p_ele0_FSR));
	      h_FSR_angle->Fill(p_ph.Angle(p_ele0_FSR.Vect()));
	      h_FSR_dr_w->Fill(p_ph.DeltaR(p_ele0_FSR),p_ph.Pt());
	      h_FSR_angle_w->Fill(p_ph.Angle(p_ele0_FSR.Vect()),p_ph.Pt());

	      h_FSR_deta_dphi->Fill(p_ph.DeltaPhi(p_ele0_FSR), p_ph.Eta()-p_ele0_FSR.Eta());
	      h_FSR_deta_dphi_w->Fill(p_ph.DeltaPhi(p_ele0_FSR), p_ph.Eta()-p_ele0_FSR.Eta(),p_ph.Pt());
	      h_FSR_dr_et->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
	      h_FSR_dr_e->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
	      
	      FSR_photons[electrons_afterFSR[0]].push_back(*p);
	  

	    }

	    if ( from_ele1 && fabs(grandmother->pdg_id())==23 ){

	      sumEt_tot_fsr1 += p_ph.Pt();

	      h_FSR_et->Fill(p_ph.Pt());
	      h_FSR_e->Fill(p_ph.E());
	      h_FSR_dr->Fill(p_ph.DeltaR(p_ele1_FSR));
	      h_FSR_angle->Fill(p_ph.Angle(p_ele1_FSR.Vect()));
	      h_FSR_dr_w->Fill(p_ph.DeltaR(p_ele1_FSR),p_ph.Pt());
	      h_FSR_angle_w->Fill(p_ph.Angle(p_ele1_FSR.Vect()),p_ph.Pt());
	  
	      h_FSR_deta_dphi->Fill(p_ph.DeltaPhi(p_ele1_FSR),p_ph.Eta()-p_ele1_FSR.Eta());
	      h_FSR_deta_dphi_w->Fill(p_ph.DeltaPhi(p_ele1_FSR), p_ph.Eta()-p_ele1_FSR.Eta(),p_ph.Pt());
	      h_FSR_dr_et->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
	      h_FSR_dr_e->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));

	      FSR_photons[electrons_afterFSR[1]].push_back(*p);

	    }

	    h_ph_mom->Fill(fabs(mother->pdg_id()));
	    if ( grandmother )
	      h_ph_grandma->Fill(fabs(grandmother->pdg_id()));

	  }
	    
            
	  double dR0 = p_ph.DeltaR(p_ele0_FSR);
	  double dR1 = p_ph.DeltaR(p_ele1_FSR);

	  double dz0 = (*p)->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();
	  if ( grandmother && grandmother->production_vertex() )
	    dz0 = grandmother->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();
	  else if ( mother && mother->production_vertex() )
	    dz0 = mother->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();

	  double dz1 = (*p)->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();
	  if ( grandmother && grandmother->production_vertex() )
	    dz1 = grandmother->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();
	  else if ( mother && mother->production_vertex() )
	    dz1 = mother->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();


	  for ( unsigned int i=0; i<n_dr_bins; ++i ) {

	    double dR_cut = 0.01 * (1 + i);


	    if ( dR0 < dR_cut ) {

	      h_pt_ph[i]->Fill(p_ph.Pt());
	      h_dz_ph[i]->Fill(dz0);
	      h_dr_et_ph[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
	      h_dr_e_ph[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
	
	      if ( fabs(mother->pdg_id()) == 11 ){

		if ( grandmother && fabs(grandmother->pdg_id()) == 23 ){
		  sumEt_fsr_ele0[i] += p_ph.Pt();
		  h_pt_ph_fsr[i]->Fill(p_ph.Pt());  
		  h_dz_ph_fsr[i]->Fill(dz0);  
		  h_dr_et_ph_fsr[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
		  h_dr_e_ph_fsr[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
		}
		else {
		  h_pt_ph_e[i]->Fill(p_ph.Pt());  
		  h_dr_et_ph_e[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
		  h_dr_e_ph_e[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
		}
	    
	      }
	      else {
	    
		sumEt_other_ele0[i] += p_ph.Pt();
		h_pt_ph_other[i]->Fill(p_ph.Pt());
		h_dz_ph_other[i]->Fill(dz0);
		h_dr_et_ph_other[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
		h_dr_e_ph_other[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));

		h_id_ph_other[i]->Fill(fabs(mother->pdg_id()));

	      }
	  
	    }
	    else if ( dR1 < dR_cut ) {

	      h_pt_ph[i]->Fill(p_ph.Pt());
	      h_dz_ph[i]->Fill(dz1);
	      h_dr_et_ph[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
	      h_dr_e_ph[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));
	
	      if ( fabs(mother->pdg_id()) == 11) {
		
		if( grandmother && fabs(grandmother->pdg_id()) == 23 ){
	      
		  sumEt_fsr_ele1[i] += p_ph.Pt();
		  h_pt_ph_fsr[i]->Fill(p_ph.Pt());  
		  h_dz_ph_fsr[i]->Fill(dz1);  
		  h_dr_et_ph_fsr[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
		  h_dr_e_ph_fsr[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));

		}
		else {
		  h_pt_ph_e[i]->Fill(p_ph.Pt());  
		  h_dr_et_ph_e[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
		  h_dr_e_ph_e[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));
		}
	      }
	      else {
	    
		sumEt_other_ele1[i] += p_ph.Pt();
		h_pt_ph_other[i]->Fill(p_ph.Pt());
		h_dr_et_ph_other[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
		h_dr_e_ph_other[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));
		h_id_ph_other[i]->Fill(fabs(mother->pdg_id()));

		h_dz_ph_other[i]->Fill(dz1);
	      }
	  
	    }

	  } // for i

	} // for p
 
      
      }


      // =====================================================================================================
      //  PU vertices


      for (unsigned int ivtx=0; ivtx<PUvtx_index.size(); ++ivtx) {

	const HepMC::GenEvent * genEvent = cfhepmcprod->getObject(PUvtx_index[ivtx]).GetEvent();

	for ( HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin(); p != genEvent->particles_end(); ++p ) {

	  if ( (*p)->status() != 1 || fabs((*p)->pdg_id()) != 22 ) continue;
      

	  // Get the photon mother and grandmother
	  HepMC::GenParticle* mother      = (*((*p)->production_vertex()->particles_in_const_begin()));
	  HepMC::GenParticle* grandmother = NULL;
	  if (mother && mother->production_vertex()) 
	    grandmother = (*(mother->production_vertex()->particles_in_const_begin()));


	  TLorentzVector p_ph( (*p)->momentum().px(),
			       (*p)->momentum().py(),
			       (*p)->momentum().pz(),
			       (*p)->momentum().e() );



	  double dR0 = p_ph.DeltaR(p_ele0_FSR);
	  double dR1 = p_ph.DeltaR(p_ele1_FSR);
	  
	  double dz0 = (*p)->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();
	  if ( grandmother && grandmother->production_vertex() )
	    dz0 = grandmother->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();
	  else if ( mother && mother->production_vertex() )
	    dz0 = mother->production_vertex()->point3d().z() - electrons_beforeFSR[0]->production_vertex()->point3d().z();

	  double dz1 = (*p)->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();
	  if ( grandmother && grandmother->production_vertex() )
	    dz1 = grandmother->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();
	  else if ( mother && mother->production_vertex() )
	    dz1 = mother->production_vertex()->point3d().z() - electrons_beforeFSR[1]->production_vertex()->point3d().z();

	  
	  for ( unsigned int i=0; i<n_dr_bins; ++i ) {

	    double dR_cut = 0.01 * (1 + i);

	    if ( dR0 < dR_cut ) {
	
	      h_pt_ph[i]->Fill(p_ph.Pt());
	      h_dz_ph[i]->Fill(dz0);
	      h_dr_et_ph[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
	      h_dr_e_ph[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
	
	      sumEt_PU_ele0[i] += p_ph.Pt();
	      h_pt_ph_PU[i]->Fill(p_ph.Pt()); 
	      h_dz_ph_PU[i]->Fill(dz0);
	      h_dr_et_ph_PU[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele0_FSR));
	      h_dr_e_ph_PU[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele0_FSR));
	
	    }
	    else if ( dR1 < dR_cut ) {

	      h_pt_ph[i]->Fill(p_ph.Pt());
	      h_dz_ph[i]->Fill(dz1);
	      h_dr_et_ph[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
	      h_dr_e_ph[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));
	
	      sumEt_PU_ele1[i] += p_ph.Pt();
	      h_pt_ph_PU[i]->Fill(p_ph.Pt()); 
	      h_dz_ph_PU[i]->Fill(dz1);
	      h_dr_et_ph_PU[i]->Fill(p_ph.Pt(),p_ph.DeltaR(p_ele1_FSR));
	      h_dr_e_ph_PU[i]->Fill(p_ph.E(),p_ph.DeltaR(p_ele1_FSR));
	
	    }
	  
	  } // for i


	} // for p


      } // for ivtx


      // =====================================================================================================
      //  Fill the histograms

      h_nvtx->Fill(cfhepmcprod->size());

      double m_ll_noFSR  = (p_ele0_noFSR+p_ele1_noFSR).M();
      double pt_ll_noFSR = (p_ele0_noFSR+p_ele1_noFSR).Pt();
      double y_ll_noFSR  = (p_ele0_noFSR+p_ele1_noFSR).Rapidity();
      double m_ll_FSR    = (p_ele0_FSR+p_ele1_FSR).M();
      double pt_ll_FSR   = (p_ele0_FSR+p_ele1_FSR).Pt();
      double y_ll_FSR    = (p_ele0_FSR+p_ele1_FSR).Rapidity();

      h_Mee_noFSR  ->Fill(m_ll_noFSR);
      h_pT_Z_noFSR ->Fill(pt_ll_noFSR);
      h_Y_Z_noFSR  ->Fill(y_ll_noFSR);
      h_pT_e1_noFSR->Fill(p_ele0_noFSR.Pt());
      h_pT_e2_noFSR->Fill(p_ele1_noFSR.Pt());

      h_Mee_FSR  ->Fill(m_ll_FSR);
      h_pT_Z_FSR ->Fill(pt_ll_FSR);
      h_Y_Z_FSR  ->Fill(y_ll_FSR);
      h_pT_e1_FSR->Fill(p_ele0_FSR.Pt());
      h_pT_e2_FSR->Fill(p_ele1_FSR.Pt());


      h_FSR_n->Fill(FSR_photons[electrons_afterFSR[0]].size());
      h_FSR_n->Fill(FSR_photons[electrons_afterFSR[1]].size());

      h_FSR_sumet->Fill(sumEt_tot_fsr0);
      h_FSR_sumet->Fill(sumEt_tot_fsr1);
    
      for ( unsigned int i=0; i<n_dr_bins; ++i ) {

	double dR_cut = 0.01 * (1 + i) - 0.005;

	double et_tot0 = sumEt_fsr_ele0[i] + sumEt_other_ele0[i] + sumEt_PU_ele0[i];
	double et_tot1 = sumEt_fsr_ele1[i] + sumEt_other_ele1[i] + sumEt_PU_ele1[i];


	if ( sumEt_tot_fsr0 > 0. ){
	  h_FSR_frac[i]->Fill(sumEt_fsr_ele0[i]/sumEt_tot_fsr0);
	  h_FSR_frac_p->Fill(dR_cut, sumEt_fsr_ele0[i]/sumEt_tot_fsr0);
	}
	h_FSR_totet[i]->Fill(sumEt_fsr_ele0[i]);
	h_FSR_totet_p->Fill(dR_cut, sumEt_fsr_ele0[i]);

	if ( et_tot0 > 0. ){
	  h_frac_fsr[i]->Fill(sumEt_fsr_ele0[i]/et_tot0);
	  h_frac_other[i]->Fill(sumEt_other_ele0[i]/et_tot0);

	  h_frac_fsr_p->Fill(dR_cut, sumEt_fsr_ele0[i]/et_tot0);
	  h_frac_other_p->Fill(dR_cut, sumEt_other_ele0[i]/et_tot0);
	  h_frac_PU_p->Fill(dR_cut, sumEt_PU_ele0[i]/et_tot0);
	  
	  h_totet_fsr[i]->Fill(sumEt_fsr_ele0[i]);
	  h_totet_other[i]->Fill(sumEt_other_ele0[i]);

	  h_totet_fsr_p->Fill(dR_cut, sumEt_fsr_ele0[i]);
	  h_totet_other_p->Fill(dR_cut, sumEt_other_ele0[i]);
	  h_totet_PU_p->Fill(dR_cut, sumEt_PU_ele0[i]);

	}
      
	if ( sumEt_tot_fsr1 > 0. ){
	  h_FSR_frac[i]->Fill(sumEt_fsr_ele1[i]/sumEt_tot_fsr1);
	  h_FSR_frac_p->Fill(dR_cut, sumEt_fsr_ele1[i]/sumEt_tot_fsr1);
	}
	h_FSR_totet[i]->Fill(sumEt_fsr_ele1[i]);
	h_FSR_totet_p->Fill(dR_cut, sumEt_fsr_ele1[i]);
      
	if ( et_tot1 > 0. ){
	  h_frac_fsr[i]->Fill(sumEt_fsr_ele1[i]/et_tot1);
	  h_frac_other[i]->Fill(sumEt_other_ele1[i]/et_tot1);
      
	  h_frac_fsr_p->Fill(dR_cut, sumEt_fsr_ele1[i]/et_tot1);
	  h_frac_other_p->Fill(dR_cut, sumEt_other_ele1[i]/et_tot1);
	  h_frac_PU_p->Fill(dR_cut, sumEt_PU_ele1[i]/et_tot1);
        
	  h_totet_fsr[i]->Fill(sumEt_fsr_ele1[i]);
	  h_totet_other[i]->Fill(sumEt_other_ele1[i]);
      
	  h_totet_fsr_p->Fill(dR_cut, sumEt_fsr_ele1[i]);
	  h_totet_other_p->Fill(dR_cut, sumEt_other_ele1[i]);
	  h_totet_PU_p->Fill(dR_cut, sumEt_PU_ele1[i]);
	  
	}

      } // for i
      

      // --- Cleanup the heap:
      if(cfhepmcprod) delete cfhepmcprod;
      
    } // if (!gotHepMC)


  } // if (isMC)


}


// ------------ method called once each job just before starting event loop  ------------
void 
FSRAnalyzer_PU::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FSRAnalyzer_PU::endJob() 
{

}  


// -------------------------------------------------------------------------
TVector3 
FSRAnalyzer_PU::extrapolateToECAL(TVector3 const &mom, TVector3 const &pos, int const &subdet)
{
  
  TVector3 intersection_point(-999.,-999.,-999.);

  // ECAL radius and half-length:
  const double R = 1290.;
  const double L = 2680.;

  // particle's direction:
  double a = mom.X()/mom.Mag();
  double b = mom.Y()/mom.Mag();
  double c = mom.Z()/mom.Mag();

  // particle's origin point:
  double x0 = pos.X(); 
  double y0 = pos.Y(); 
  double z0 = pos.Z(); 

  if ( subdet == BARREL ) {

    double t1 = (-sqrt(a*a*R*R+b*b*R*R-b*b*x0*x0-a*a*y0*y0+2.*a*b*x0*y0) - a*x0 - b*y0)/(a*a+b*b);
    double t2 = ( sqrt(a*a*R*R+b*b*R*R-b*b*x0*x0-a*a*y0*y0+2.*a*b*x0*y0) - a*x0 - b*y0)/(a*a+b*b);

    if ( (z0+c*t1)*c > 0. )
      intersection_point.SetXYZ(x0+a*t1, y0+b*t1, z0+c*t1);
    else
      intersection_point.SetXYZ(x0+a*t2, y0+b*t2, z0+c*t2);

  }
  else if ( subdet == ENDCAP ) {

    if ( mom.Eta() > 0. ) {
      
      double t = (L-z0)/c;
      intersection_point.SetXYZ(x0+a*t, y0+b*t, L);

    }
    else {
    
      double t = (-L-z0)/c;
      intersection_point.SetXYZ(x0+a*t, y0+b*t, -L);

    }
  
  }   

  return intersection_point;
  
}


// -------------------------------------------------------------------------
bool 
FSRAnalyzer_PU::is_in_ellipse(HepMC::GenParticle const &ph, HepMC::GenParticle const &ele, double const &deltaR)
{

  TVector3 ph_mom( ph.momentum().px(), ph.momentum().py(), ph.momentum().pz() );
  TVector3 ph_pos( ph.production_vertex()->position().x(), 
		   ph.production_vertex()->position().y(),
		   ph.production_vertex()->position().z() );

  TVector3 ele_mom( ele.momentum().px(), ele.momentum().py(), ele.momentum().pz() );
  TVector3 ele_pos( ele.production_vertex()->position().x(), 
		    ele.production_vertex()->position().y(),
		    ele.production_vertex()->position().z() );

  const int ph_det  = ( fabs(ph_mom.Eta()) <1.479 ? BARREL : ENDCAP );
  const int ele_det = ( fabs(ele_mom.Eta())<1.479 ? BARREL : ENDCAP );

  const double ele_eta = ele_mom.Eta();


  // Get the electron and photon intersection points on the ECAL surface
  TVector3 ph_intersection  = extrapolateToECAL( ph_mom, ph_pos, ph_det );
  TVector3 ele_intersection = extrapolateToECAL( ele_mom, ele_pos, ele_det );


  // if the electron points to the ECAL barrel:
  if ( fabs(ele_eta) <= 1.479 ) {
    
    double angle = TMath::PiOver2()-ele_intersection.Phi();
    ele_intersection.RotateZ(angle);
    ph_intersection.RotateZ(angle);
    ph_mom.RotateZ(angle);
    ph_pos.RotateZ(angle);
    ele_mom.RotateZ(angle);
    ele_pos.RotateZ(angle);

    if ( ele_eta > 0. )
      ele_mom.SetMagThetaPhi(ele_mom.Mag(),ele_mom.Theta()+deltaR,ele_mom.Phi());
    else
      ele_mom.SetMagThetaPhi(ele_mom.Mag(),ele_mom.Theta()-deltaR,ele_mom.Phi());
    
    TVector3 ellipse_perigee = extrapolateToECAL( ele_mom, ele_pos, ele_det );


    if ( ele_eta > 0. )
      ele_mom.SetMagThetaPhi(ele_mom.Mag(),ele_mom.Theta()-2.*deltaR,ele_mom.Phi());
    else
      ele_mom.SetMagThetaPhi(ele_mom.Mag(),ele_mom.Theta()+2.*deltaR,ele_mom.Phi());
    
    TVector3 ellipse_apogee = extrapolateToECAL( ele_mom, ele_pos, ele_det );


    const double b = ele_intersection.Perp()*TMath::Tan(deltaR);
    const double a = 0.5*fabs(ellipse_apogee.Z()-ellipse_perigee.Z());

    const double x0 = ele_intersection.X();
    const double z0 = 0.5*(ellipse_apogee.Z()+ellipse_perigee.Z());
    
    if ( ( (ph_intersection.X()-x0)*(ph_intersection.X()-x0)/(b*b) +   
    	   (ph_intersection.Z()-z0)*(ph_intersection.Z()-z0)/(a*a) ) < 1. ) 
      return true;
    
  }
  else {

    const double R  = ele_intersection.Z()*TMath::Tan(deltaR);
    const double x0 = ele_intersection.X();
    const double y0 = ele_intersection.Y();

    if ( ( (ph_intersection.X()-x0)*(ph_intersection.X()-x0) +   
	   (ph_intersection.Y()-y0)*(ph_intersection.Y()-y0) ) < R*R ) 
      return true;
    
  }

  return false;

}

// ------------ method called when starting to processes a run  ------------
void 
FSRAnalyzer_PU::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FSRAnalyzer_PU::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FSRAnalyzer_PU::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FSRAnalyzer_PU::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FSRAnalyzer_PU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FSRAnalyzer_PU);
