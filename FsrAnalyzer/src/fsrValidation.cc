// -*- C++ -*-
//

#include "Fsr/FsrAnalyzer/interface/fsrValidation.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Common/interface/RefVector.h"

//
// member functions
//

class GreaterPt{
public:
  bool operator()( const math::XYZTLorentzVector& a, const math::XYZTLorentzVector& b) {
    return a.Pt() > b.Pt();
  }
};

double
fsrValidation::distR(TLorentzVector a ,TLorentzVector b){

   double deltaPhi = fabs(a.Phi()-b.Phi());
   if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
   double delta = sqrt( deltaPhi*deltaPhi  + ((a.Eta()-b.Eta())*(a.Eta()-b.Eta())));
   return delta;
}
// ------------ method called for each event  ------------
void
fsrValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool Debug=false;
  double cat1=0.05;
  double cat2=2.0;
  double cat3=0;
 
   ///////////////////////
   ///// Z Analysis
   ///////////////////////

    edm::Handle<reco::GsfElectronCollection > goodEPair;
   edm::Handle<reco::PFCandidateCollection> goodPfEPair;
   int checkGoodEPairSize=0;
   if (!usingPF){
      iEvent.getByLabel (goodEPairTag, goodEPair);
      checkGoodEPairSize=goodEPair->size();}
   else {
      iEvent.getByLabel (goodEPairTag, goodPfEPair);
      checkGoodEPairSize=goodPfEPair->size();}
   
   edm::Handle<reco::GenParticleCollection> genPart;
   if (usingMC){
      iEvent.getByLabel (genParticleCollection_,genPart);
   }


   if (checkGoodEPairSize==2)
 
   {   
      //double dist=0.,maxDist=0.05;
      int eGenCount=0;
      //int g1GenCount=0, g2GenCount=0;
      
           
      //TLorentz vector of the two Z boson electrons, at GEN level
     TLorentzVector e1_gen,e2_gen;
     TLorentzVector e1_genS3,e2_genS3;
     gamma1_gen.clear();
     gamma2_gen.clear();

      //TLorentz vector of the two Z boson electrons, at RECO level
      TLorentzVector e1_reco, e2_reco;
      
      if (!usingPF){     
	 reco::GsfElectronCollection::const_iterator it=goodEPair->begin();     
	 e1_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 it++;
	 e2_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
      }
      else {     
	 reco::PFCandidateCollection::const_iterator it=goodPfEPair->begin();     
	 e1_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 it++;
	 e2_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
      }     
      if (usingMC){
	 reco::GenParticleCollection::const_iterator it_e1_gen;
	 reco::GenParticleCollection::const_iterator it_e2_gen;
	 for(reco::GenParticleCollection::const_iterator itgen=genPart->begin();itgen!=genPart->end();itgen++)
	 {
	    bool secondChoice = false;
	    bool isRightEle = false;
	    if ( fabs(itgen->pdgId())==11 && itgen->status()==1 ){
	       isRightEle = searchHistory(&(*itgen));
	       if (isRightEle){
		  TLorentzVector e_gen;
		  e_gen.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
		  double deltaR1= distR(e_gen,e1_reco);
		  double deltaR2= distR(e_gen,e2_reco);	
		  if (deltaR1 < deltaConeGen ){
		     e1_gen = e_gen;
		     it_e1_gen= itgen;
		     e1_genS3 = e1_test;
		     eGenCount++;
		  } 
		  else if (deltaR2 < deltaConeGen ){
		     e2_gen = e_gen;
		     it_e2_gen= itgen;
		     e2_genS3 = e1_test;
		     eGenCount++;
		  }
		  else if (useSecondChoice){
		     if (deltaR1 < deltaR2){
			e1_gen = e_gen;
			it_e1_gen= itgen;
			e1_genS3 = e1_test;
			eGenCount++;
			secondChoice = true;
		     }
		     else if (deltaR1 > deltaR2){
			e2_gen = e_gen;
			it_e2_gen= itgen;
			e2_genS3 = e1_test;
			eGenCount++;
			secondChoice = true;
		     }
		  }
		  if (secondChoice) cout << "Manovra di emergenza - distanza da 1 = "<< deltaR1<<" - distanza da 2 = "<< deltaR2<<endl;
	       }
	    }
	 }
	 //if (secondChoice) cout << "Ho dovuto usare la manovra di emergenza"<<endl;
	 if (eGenCount>2) cout << "Ho trovato piu' di due electtroni a livello GEN nel run "<< iEvent.id().event()<< endl;
	 if (eGenCount<2) cout << "Ho trovato meno di due electtroni a livello GEN nel run "<< iEvent.id().event()<< endl;

	 //if (it_e1_gen->mother()->status()==3) gamma1_gen.SetPtEtaPhiM(0,0,0,0);
	 //if (it_e2_gen->mother()->status()==3) gamma2_gen.SetPtEtaPhiM(0,0,0,0);
	 
	 if (eGenCount>=2 && (it_e1_gen->mother()->status()!=3 || it_e2_gen->mother()->status()!=3)) {
	    
	    for(reco::GenParticleCollection::const_iterator itgen=genPart->begin();itgen!=genPart->end();itgen++){
		bool isRightPhoton=false;
	       if ( fabs(itgen->pdgId())==22 && itgen->status()==1){
		  isRightPhoton=searchHistory(&(*itgen));
		  if (isRightPhoton){
		    TLorentzVector a;
		    a.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());		    
		    double deltaR1= distR(e1_test,e1_genS3);
		    double deltaR2= distR(e1_test,e2_genS3);
		    if (deltaR1<deltaR2) {
		       gamma1_gen.push_back(a);
		       if (deltaR1>0.05) cout << "Qui la distanza e' troppo grande - elettrone 1"<<endl;}
		    else if (deltaR1>deltaR2){
		       gamma2_gen.push_back(a);
		       if (deltaR2>0.05) cout << "Qui la distanza e' troppo grande - elettrone 2"<<endl;}
		  }
	       }
	    }
	 }

	 // Ho trovato Ereco, Egen, EgenS3, e i Ggen, ora posso fare i plot  ===========================
	 
	 h_ptReco->Fill(e1_reco.Pt());
	 h_ptReco->Fill(e2_reco.Pt());
	 h_ptGen->Fill(e1_gen.Pt());
	 h_ptGen->Fill(e2_gen.Pt());
	 h_ptGen3->Fill(e1_genS3.Pt());
	 h_ptGen3->Fill(e2_genS3.Pt());
	 h_ptGen13->Fill(e1_genS3.Pt() - e1_gen.Pt());
	 h_ptGen13->Fill(e2_genS3.Pt() - e2_gen.Pt());
	 double sumPtG1 =0, 
	    sumPt05G1=0, sumPt06G1=0, sumPt07G1=0, sumPt08G1=0, sumPt09G1=0, sumPt10G1=0, sumPt11G1=0;
	 for (vector<TLorentzVector>::const_iterator itG = gamma1_gen.begin(); itG != gamma1_gen.end(); itG++){
	    sumPtG1=sumPtG1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.05) sumPt05G1=sumPt05G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.06) sumPt06G1=sumPt06G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.07) sumPt07G1=sumPt07G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.08) sumPt08G1=sumPt08G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.09) sumPt09G1=sumPt09G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.10) sumPt10G1=sumPt10G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.11) sumPt11G1=sumPt11G1 + itG->Pt();
	 }
	 double sumPtG2 =0, 
	    sumPt05G2=0, sumPt06G2=0, sumPt07G2=0, sumPt08G2=0, sumPt09G2=0, sumPt10G2=0, sumPt11G2=0;
	 for (vector<TLorentzVector>::const_iterator itG = gamma2_gen.begin(); itG != gamma2_gen.end(); itG++){
	    sumPtG2=sumPtG2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.05) sumPt05G2=sumPt05G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.06) sumPt06G2=sumPt06G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.07) sumPt07G2=sumPt07G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.08) sumPt08G2=sumPt08G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.09) sumPt09G2=sumPt09G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.10) sumPt10G2=sumPt10G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.11) sumPt11G2=sumPt11G2 + itG->Pt();
	 }
	 h_sumPtGamma->Fill(sumPtG1);
	 h_sumPtGamma->Fill(sumPtG2);
	 h_sumPtGammaDr05->Fill(sumPt05G1);
	 h_sumPtGammaDr05->Fill(sumPt05G2);
	 h_sumPtGammaDr06->Fill(sumPt06G1);
	 h_sumPtGammaDr06->Fill(sumPt06G2);
	 h_sumPtGammaDr07->Fill(sumPt07G1);
	 h_sumPtGammaDr07->Fill(sumPt07G2);
	 h_sumPtGammaDr08->Fill(sumPt08G1);
	 h_sumPtGammaDr08->Fill(sumPt08G2);
	 h_sumPtGammaDr09->Fill(sumPt09G1);
	 h_sumPtGammaDr09->Fill(sumPt09G2);
	 h_sumPtGammaDr10->Fill(sumPt10G1);
	 h_sumPtGammaDr10->Fill(sumPt10G2);
	 h_sumPtGammaDr11->Fill(sumPt11G1);
	 h_sumPtGammaDr11->Fill(sumPt11G2);

	 if (sumPtG1<cat1){	    
	    h_ptRecoG1->Fill(e1_reco.Pt());
	    h_ptRecoPtGenG1->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenVsPtGenG1->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	 } else if (sumPtG1< cat2){
	    h_ptRecoG2->Fill(e1_reco.Pt());
	    h_ptRecoPtGenG2->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenVsPtGenG2->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	 }else {	    
	    h_ptRecoG3->Fill(e1_reco.Pt());
	    h_ptRecoPtGenG3->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenVsPtGenG3->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	 }

	 if (sumPtG2<cat1){	    
	    h_ptRecoG1->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG1->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG1->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	 } else if (sumPtG2< cat2){
	    h_ptRecoG2->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG2->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG2->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	 }else {	    
	    h_ptRecoG3->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG3->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG3->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	 }
      }

   } else if (Debug){std::cout << "WARNING: More than two electron selected"<< std::endl;}  

}


// ------------ method called once each job just before starting event loop  ------------
void 
fsrValidation::beginJob()
{

  cout<<endl;
  cout<<"##############################"<<endl;
  cout<<"#   Studies Parameters   #"<<endl;
  cout<<"##############################"<<endl;
  cout<<endl; 
  cout<<"DR isolation cone at gen Level="<<deltaConeGen<<endl;
  cout<<endl;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
fsrValidation::endJob() 
{
    
 //Histograms
 

}

// ------------ method called when starting to processes a run  ------------
void 
fsrValidation::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
fsrValidation::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
fsrValidation::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
fsrValidation::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(fsrValidation);
