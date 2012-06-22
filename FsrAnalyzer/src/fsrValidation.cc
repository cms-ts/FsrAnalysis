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
   
   ///////////////////////
   ///// Z Analysis
   ///////////////////////
   Handle<reco::PFCandidateCollection> particleCollection;
   iEvent.getByLabel(particleCollection_,particleCollection);

   // iso deposits
   IsoDepositVals isoVals(isoValInputTags_.size());
   for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
      iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
   }
   
   const IsoDepositVals * IsoVals = &isoVals;

   double lepIsoRho;  
   /////// Pileup density "rho" for lepton isolation subtraction /////
   //cout << "1- removing PU ..."<<endl;
   edm::Handle<double> rhoLepIso;
   const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
   iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
   if( *rhoLepIso == *rhoLepIso)  lepIsoRho = *rhoLepIso;
   else  lepIsoRho =  -999999.9;   
  
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
     
     double A_eff_PH, A_eff_NH;
      //TLorentz vector of the two Z boson electrons, at RECO level
      TLorentzVector e1_reco, e2_reco;
      double iso_e1, iso_e2;
      double iso_e1V05, iso_e2V05,iso_e1V08, iso_e2V08, iso_e1V10, iso_e2V10;

      reco::PFCandidateCollection::const_iterator e1recoit, e2recoit; 

      if (!usingPF){     
	 reco::GsfElectronCollection::const_iterator it=goodEPair->begin();     
	 e1_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 double isoPhotonUser=0;
	 for(reco::PFCandidateCollection::const_iterator itPart = particleCollection->begin(); 
	     itPart!= particleCollection->end(); ++itPart) {
	    if (itPart->pdgId()==22){
	       double deltaPhi = fabs (it->phi() - itPart->phi());
	       if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
	       double deltaR = sqrt ( pow(it->eta() - itPart->eta(),2) + deltaPhi*deltaPhi);
	       if (deltaR < 0.3 
		   //&& itPart->pt()>0.5 
		  ) isoPhotonUser += itPart->pt();
	    }
	 }
	 //iso_e1 = isoPhotonUser;
	 // reco::GsfElectronRef myElectronRef(goodEPair,0);
// 	 //double charged = (*IsoVals[0])[it->originalObjectRef()];//myElectronRef
// 	 //double photon  = (*IsoVals[1])[it->originalObjectRef()]; //myElectronRef
// 	 double photon  = (*(*IsoVals)[1])[myElectronRef]; //myElectronRef
// 	 //double neutral = (*IsoVals[2])[it->originalObjectRef()];//myElectronRef 
	 
	 // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
	 if(abs(it->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
	 else if(abs(it->eta())>1.0 && abs(it->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>1.479 && abs(it->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>2.0 && abs(it->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.2 && abs(it->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.3 && abs(it->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
	 else {A_eff_PH=0.11 ; A_eff_NH=0.021;} 
// 	 iso_e1 = max(photon - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e1 = max(isoPhotonUser - lepIsoRho*A_eff_PH  , 0.)/it->pt();

	 it++;
	 e2_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 isoPhotonUser=0;
	 for(reco::PFCandidateCollection::const_iterator itPart = particleCollection->begin(); 
	     itPart!= particleCollection->end(); ++itPart) {
	    if (itPart->pdgId()==22){
	       double deltaPhi = fabs (it->phi() - itPart->phi());
	       if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
	       double deltaR = sqrt ( pow(it->eta() - itPart->eta(),2) + deltaPhi*deltaPhi);
	       if (deltaR < 0.3 
		   //&& itPart->pt()>0.5 
		  ) isoPhotonUser += itPart->pt();
	    }
	 }
	 //iso_e2 = isoPhotonUser;
	 // reco::GsfElectronRef myElectronRef1(goodEPair,1);
// 	 //double charged = (*IsoVals[0])[it->originalObjectRef()];//myElectronRef
// 	 //double photon  = (*IsoVals[1])[it->originalObjectRef()]; //myElectronRef
// 	 photon  = (*(*IsoVals)[1])[myElectronRef1]; //myElectronRef
// 	 //double neutral = (*IsoVals[2])[it->originalObjectRef()];//myElectronRef 
	 
	 // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
	 if(abs(it->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
	 else if(abs(it->eta())>1.0 && abs(it->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>1.479 && abs(it->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>2.0 && abs(it->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.2 && abs(it->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.3 && abs(it->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
	 else {A_eff_PH=0.11 ; A_eff_NH=0.021;} 
// 	 iso_e2 = max(photon - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e2 = max(isoPhotonUser - lepIsoRho*A_eff_PH  , 0.)/it->pt();	 
      }
      else {     
	 reco::PFCandidateCollection::const_iterator it=goodPfEPair->begin(); 
	 e1recoit = it;
	 e1_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 double isoPhotonUser=0;
	 double isoPhotonUserV05=0, isoPhotonUserV08=0, isoPhotonUserV10=0;
	 for(reco::PFCandidateCollection::const_iterator itPart = particleCollection->begin(); 
	     itPart!= particleCollection->end(); ++itPart) {
	    if (itPart->pdgId()==22){
	       double deltaPhi = fabs (it->phi() - itPart->phi());
	       if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
	       double deltaR = sqrt ( pow(it->eta() - itPart->eta(),2) + deltaPhi*deltaPhi);
	       if (deltaR < 0.3 
		   //&& itPart->pt()>0.5 
		  ) {
		  isoPhotonUser += itPart->pt();
		  if (deltaR > 0.05) isoPhotonUserV05 += itPart->pt();
		  if (deltaR > 0.08) isoPhotonUserV08 += itPart->pt();
		  if (deltaR > 0.10) isoPhotonUserV10 += itPart->pt();
	       }
	    }
	 }
	 //iso_e1 = isoPhotonUser;
	 // cout << "control 1 "<<endl;
// 	 reco::PFCandidateRef myElectronRef(goodPfEPair,0);
// 	 cout << "control 2 "<<endl;
// 	 //double charged = (*IsoVals[0])[it->originalObjectRef()];//myElectronRef
// 	 //double photon  = (*IsoVals[1])[it->originalObjectRef()]; //myElectronRef
// 	 double photon  = (*(*IsoVals)[1])[myElectronRef]; //myElectronRef
// 	 cout << "control 3 "<<endl;
// 	 //double neutral = (*IsoVals[2])[it->originalObjectRef()];//myElectronRef 
	 
	 // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
	 if(abs(it->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
	 else if(abs(it->eta())>1.0 && abs(it->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>1.479 && abs(it->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>2.0 && abs(it->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.2 && abs(it->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.3 && abs(it->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
	 else {A_eff_PH=0.11 ; A_eff_NH=0.021;} 
// 	 iso_e1 = max(photon - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e1 = max(isoPhotonUser - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e1V05 = max(isoPhotonUserV05 - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e1V08 = max(isoPhotonUserV08 - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e1V10 = max(isoPhotonUserV10 - lepIsoRho*A_eff_PH  , 0.)/it->pt();

	 it++;
	 e2recoit= it;
	 e2_reco.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());
	 isoPhotonUser=0; isoPhotonUserV05=0; isoPhotonUserV08=0; isoPhotonUserV10=0;
	 for(reco::PFCandidateCollection::const_iterator itPart = particleCollection->begin(); 
	     itPart!= particleCollection->end(); ++itPart) {
	    if (itPart->pdgId()==22){
	       double deltaPhi = fabs (it->phi() - itPart->phi());
	       if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
	       double deltaR = sqrt ( pow(it->eta() - itPart->eta(),2) + deltaPhi*deltaPhi);
	       if (deltaR < 0.3 
		   //&& itPart->pt()>0.5 
		  ) {
		  isoPhotonUser += itPart->pt();
		  if (deltaR > 0.05) isoPhotonUserV05 += itPart->pt();
		  if (deltaR > 0.08) isoPhotonUserV08 += itPart->pt();
		  if (deltaR > 0.10) isoPhotonUserV10 += itPart->pt();
	       }
	    }
	 }
	 //iso_e2 = isoPhotonUser;	 	 
	 // reco::PFCandidateRef myElectronRef1(goodPfEPair,1);
// 	 //double charged = (*IsoVals[0])[it->originalObjectRef()];//myElectronRef
// 	 //double photon  = (*IsoVals[1])[it->originalObjectRef()]; //myElectronRef
// 	 photon  = (*(*IsoVals)[1])[myElectronRef1]; //myElectronRef
// 	 //double neutral = (*IsoVals[2])[it->originalObjectRef()];//myElectronRef 
	 
	 // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
	 if(abs(it->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
	 else if(abs(it->eta())>1.0 && abs(it->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>1.479 && abs(it->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
	 else if(abs(it->eta())>2.0 && abs(it->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.2 && abs(it->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
	 else if(abs(it->eta())>2.3 && abs(it->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
	 else {A_eff_PH=0.11 ; A_eff_NH=0.021;} 
// 	 iso_e2 = max(photon - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e2 = max(isoPhotonUser - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e2V05 = max(isoPhotonUserV05 - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e2V08 = max(isoPhotonUserV08 - lepIsoRho*A_eff_PH  , 0.)/it->pt();
 	 iso_e2V10 = max(isoPhotonUserV10 - lepIsoRho*A_eff_PH  , 0.)/it->pt();
      } 


      if (usingMC){
	 reco::GenParticleCollection::const_iterator it_e1_gen;
	 reco::GenParticleCollection::const_iterator it_e2_gen;
	 double sumPtAll05G1=0, sumPtAll06G1=0, sumPtAll07G1=0, sumPtAll08G1=0, sumPtAll09G1=0, sumPtAll10G1=0, sumPtAll11G1=0, 
	    sumPtAll13G1=0, sumPtAll15G1=0, sumPtAll17G1=0, sumPtAll20G1=0, sumPtAll25G1=0, sumPtAll30G1=0;
	 double sumPtAll05G2=0, sumPtAll06G2=0, sumPtAll07G2=0, sumPtAll08G2=0, sumPtAll09G2=0, sumPtAll10G2=0, sumPtAll11G2=0, 
	    sumPtAll13G2=0, sumPtAll15G2=0, sumPtAll17G2=0, sumPtAll20G2=0, sumPtAll25G2=0, sumPtAll30G2=0;
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
	    
	    if ( itgen->pdgId()==22 && itgen->status()==1 ){
	       TLorentzVector photonG;
	       photonG.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
	       double deltaR1= distR(photonG,e1_reco);
	       double deltaR2= distR(photonG,e2_reco);	        
	       if (deltaR1<0.05) sumPtAll05G1=sumPtAll05G1 + itgen->pt();
	       if (deltaR1<0.06) sumPtAll06G1=sumPtAll06G1 + itgen->pt();
	       if (deltaR1<0.07) sumPtAll07G1=sumPtAll07G1 + itgen->pt();
	       if (deltaR1<0.08) sumPtAll08G1=sumPtAll08G1 + itgen->pt();
	       if (deltaR1<0.09) sumPtAll09G1=sumPtAll09G1 + itgen->pt();
	       if (deltaR1<0.10) sumPtAll10G1=sumPtAll10G1 + itgen->pt();
	       if (deltaR1<0.11) sumPtAll11G1=sumPtAll11G1 + itgen->pt();
	       if (deltaR1<0.13) sumPtAll13G1=sumPtAll13G1 + itgen->pt();
	       if (deltaR1<0.15) sumPtAll15G1=sumPtAll15G1 + itgen->pt();
	       if (deltaR1<0.17) sumPtAll17G1=sumPtAll17G1 + itgen->pt();
	       if (deltaR1<0.20) sumPtAll20G1=sumPtAll20G1 + itgen->pt();
	       if (deltaR1<0.25) sumPtAll25G1=sumPtAll25G1 + itgen->pt();
	       if (deltaR1<0.30) sumPtAll30G1=sumPtAll30G1 + itgen->pt();

	       if (deltaR2<0.05) sumPtAll05G2=sumPtAll05G2 + itgen->pt();
	       if (deltaR2<0.06) sumPtAll06G2=sumPtAll06G2 + itgen->pt();
	       if (deltaR2<0.07) sumPtAll07G2=sumPtAll07G2 + itgen->pt();
	       if (deltaR2<0.08) sumPtAll08G2=sumPtAll08G2 + itgen->pt();
	       if (deltaR2<0.09) sumPtAll09G2=sumPtAll09G2 + itgen->pt();
	       if (deltaR2<0.10) sumPtAll10G2=sumPtAll10G2 + itgen->pt();
	       if (deltaR2<0.11) sumPtAll11G2=sumPtAll11G2 + itgen->pt();
	       if (deltaR2<0.13) sumPtAll13G2=sumPtAll13G2 + itgen->pt();
	       if (deltaR2<0.15) sumPtAll15G2=sumPtAll15G2 + itgen->pt();
	       if (deltaR2<0.17) sumPtAll17G2=sumPtAll17G2 + itgen->pt();
	       if (deltaR2<0.20) sumPtAll20G2=sumPtAll20G2 + itgen->pt();
	       if (deltaR2<0.25) sumPtAll25G2=sumPtAll25G2 + itgen->pt();
	       if (deltaR2<0.30) sumPtAll30G2=sumPtAll30G2 + itgen->pt();
	    }
	 }	 
	 h_sumPtAllGammaDr05->Fill(sumPtAll05G1);
	 h_sumPtAllGammaDr05->Fill(sumPtAll05G2);
	 h_sumPtAllGammaDr06->Fill(sumPtAll06G1);
	 h_sumPtAllGammaDr06->Fill(sumPtAll06G2);
	 h_sumPtAllGammaDr07->Fill(sumPtAll07G1);
	 h_sumPtAllGammaDr07->Fill(sumPtAll07G2);
	 h_sumPtAllGammaDr08->Fill(sumPtAll08G1);
	 h_sumPtAllGammaDr08->Fill(sumPtAll08G2);
	 h_sumPtAllGammaDr09->Fill(sumPtAll09G1);
	 h_sumPtAllGammaDr09->Fill(sumPtAll09G2);
	 h_sumPtAllGammaDr10->Fill(sumPtAll10G1);
	 h_sumPtAllGammaDr10->Fill(sumPtAll10G2);
	 h_sumPtAllGammaDr11->Fill(sumPtAll11G1);
	 h_sumPtAllGammaDr11->Fill(sumPtAll11G2);
	 h_sumPtAllGammaDr13->Fill(sumPtAll13G1);
	 h_sumPtAllGammaDr13->Fill(sumPtAll13G2);
	 h_sumPtAllGammaDr15->Fill(sumPtAll15G1);
	 h_sumPtAllGammaDr15->Fill(sumPtAll15G2);
	 h_sumPtAllGammaDr17->Fill(sumPtAll17G1);
	 h_sumPtAllGammaDr17->Fill(sumPtAll17G2);
	 h_sumPtAllGammaDr20->Fill(sumPtAll20G1);
	 h_sumPtAllGammaDr20->Fill(sumPtAll20G2);
	 h_sumPtAllGammaDr25->Fill(sumPtAll25G1);
	 h_sumPtAllGammaDr25->Fill(sumPtAll25G2);
	 h_sumPtAllGammaDr30->Fill(sumPtAll30G1);
	 h_sumPtAllGammaDr30->Fill(sumPtAll30G2);

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

	 double sumPtG1 =0; 
	 double cosT_g=0;	 
	 double cosT_e = -2 * TMath::ATan(-TMath::Exp((e1_gen.Eta())));
	 double sumPt05G1=0, sumPt06G1=0, sumPt07G1=0, sumPt08G1=0, sumPt09G1=0, sumPt10G1=0, sumPt11G1=0, 
	    sumPt13G1=0, sumPt15G1=0, sumPt17G1=0, sumPt20G1=0, sumPt25G1=0, sumPt30G1=0;
	 double iso_fsr1_Dr008=0,iso_fsr1_Dr005=0,iso_fsr1_Dr01=0;
	 vector<double> cosT_g1;
	 cosT_g1.clear();
	 for (vector<TLorentzVector>::const_iterator itG = gamma1_gen.begin(); itG != gamma1_gen.end(); itG++){
	     
	     sumPtG1=sumPtG1 + itG->Pt();
	     cosT_g = -2 * TMath::ATan(-TMath::Exp((itG->Eta())));
	     cosT_g1.push_back(TMath::Cos((cosT_e - cosT_g)));

	    if (distR((*itG),e1_reco)<0.05) sumPt05G1=sumPt05G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.06) sumPt06G1=sumPt06G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.07) sumPt07G1=sumPt07G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.08) sumPt08G1=sumPt08G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.09) sumPt09G1=sumPt09G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.10) sumPt10G1=sumPt10G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.11) sumPt11G1=sumPt11G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.13) sumPt13G1=sumPt13G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.15) sumPt15G1=sumPt15G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.17) sumPt17G1=sumPt17G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.20) sumPt20G1=sumPt20G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.25) sumPt25G1=sumPt25G1 + itG->Pt();
	    if (distR((*itG),e1_reco)<0.30) sumPt30G1=sumPt30G1 + itG->Pt();

	    if (distR((*itG),e1_reco)>0.08) iso_fsr1_Dr008 +=itG->Pt();
	    if (distR((*itG),e1_reco)>0.05) iso_fsr1_Dr005 +=itG->Pt();
	    if (distR((*itG),e1_reco)>0.1)  iso_fsr1_Dr01 +=itG->Pt();
	    

	 }
	 double sumPtG2 =0, 
	    sumPt05G2=0, sumPt06G2=0, sumPt07G2=0, sumPt08G2=0, sumPt09G2=0, sumPt10G2=0, sumPt11G2=0, 
	    sumPt13G2=0, sumPt15G2=0, sumPt17G2=0, sumPt20G2=0, sumPt25G2=0, sumPt30G2=0;
	 double iso_fsr2_Dr008=0,iso_fsr2_Dr005=0,iso_fsr2_Dr01=0;
	 cosT_e = -2 * TMath::ATan(-TMath::Exp((e2_gen.Eta())));
	 vector<double> cosT_g2;
	 cosT_g2.clear();
	 for (vector<TLorentzVector>::const_iterator itG = gamma2_gen.begin(); itG != gamma2_gen.end(); itG++){
	    sumPtG2=sumPtG2 + itG->Pt();
	    cosT_g = -2 * TMath::ATan(-TMath::Exp((itG->Eta())));
	    cosT_g2.push_back(TMath::Cos((cosT_e - cosT_g)));
	    if (distR((*itG),e2_reco)<0.05) sumPt05G2=sumPt05G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.06) sumPt06G2=sumPt06G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.07) sumPt07G2=sumPt07G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.08) sumPt08G2=sumPt08G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.09) sumPt09G2=sumPt09G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.10) sumPt10G2=sumPt10G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.11) sumPt11G2=sumPt11G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.13) sumPt13G2=sumPt13G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.15) sumPt15G2=sumPt15G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.17) sumPt17G2=sumPt17G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.20) sumPt20G2=sumPt20G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.25) sumPt25G2=sumPt25G2 + itG->Pt();
	    if (distR((*itG),e2_reco)<0.30) sumPt30G2=sumPt30G2 + itG->Pt();

	    if (distR((*itG),e2_reco)>0.08) iso_fsr2_Dr008 +=itG->Pt();
	    if (distR((*itG),e2_reco)>0.05) iso_fsr2_Dr005 +=itG->Pt();
	    if (distR((*itG),e2_reco)>0.1)  iso_fsr2_Dr01 +=itG->Pt();

	 }
	
	 h_iso -> Fill(iso_e1);
	 h_iso -> Fill(iso_e2);
	 h_sumPtGamma->Fill(sumPtG1);
	 h_sumPtGamma->Fill(sumPtG2);
	 
	 h_iso_fsr_Dr008 -> Fill(iso_fsr1_Dr008, iso_e1);
	 h_iso_fsr_Dr005 -> Fill(iso_fsr1_Dr005, iso_e1);
	 h_iso_fsr_Dr01 -> Fill(iso_fsr1_Dr01, iso_e1);
	 h_iso_fsr_Dr008 -> Fill(iso_fsr2_Dr008, iso_e2);
	 h_iso_fsr_Dr005 -> Fill(iso_fsr2_Dr005, iso_e2);
	 h_iso_fsr_Dr01 -> Fill(iso_fsr2_Dr01, iso_e2);

	 h_isoV05 -> Fill(iso_e1V05);
	 h_isoV05 -> Fill(iso_e2V05);
	 h_isoV08 -> Fill(iso_e1V08);
	 h_isoV08 -> Fill(iso_e2V08);
	 h_isoV10 -> Fill(iso_e1V10);
	 h_isoV10 -> Fill(iso_e2V10);
	 h_iso_fsrPt -> Fill(sumPtG1/e1_reco.Pt(), iso_e1);
	 h_iso_fsrPt_Dr008 -> Fill(sumPt08G1/e1_reco.Pt(), iso_e1);
	 h_iso_fsrPt_Dr005 -> Fill(sumPt05G1/e1_reco.Pt(), iso_e1);
	 h_iso_fsrPt_Dr01 -> Fill(sumPt10G1/e1_reco.Pt(), iso_e1);
	 h_iso_fsrPt -> Fill(sumPtG2/e2_reco.Pt(), iso_e2);
	 h_iso_fsrPt_Dr008 -> Fill(sumPt08G2/e2_reco.Pt(), iso_e2);
	 h_iso_fsrPt_Dr005 -> Fill(sumPt05G2/e2_reco.Pt(), iso_e2);
	 h_iso_fsrPt_Dr01 -> Fill(sumPt10G2/e2_reco.Pt(), iso_e2);

	 h_isoIsov_fsr_Dr008 -> Fill(sumPt08G1/e1_reco.Pt(), iso_e1-iso_e1V08);
	 h_isoIsov_fsr_Dr005 -> Fill(sumPt05G1/e1_reco.Pt(), iso_e1-iso_e1V05);
	 h_isoIsov_fsr_Dr01 -> Fill(sumPt10G1/e1_reco.Pt(), iso_e1-iso_e1V10);
	 h_isoIsov_ePt_Dr008 -> Fill(e1_reco.Pt(), iso_e1-iso_e1V08);
	 h_isoIsov_ePt_Dr005 -> Fill(e1_reco.Pt(), iso_e1-iso_e1V05);
	 h_isoIsov_ePt_Dr01 -> Fill(e1_reco.Pt(), iso_e1-iso_e1V10);
	 h_isoIsov_fsr_Dr008 -> Fill(sumPt08G2/e2_reco.Pt(), iso_e2-iso_e2V08);
	 h_isoIsov_fsr_Dr005 -> Fill(sumPt05G2/e2_reco.Pt(), iso_e2-iso_e2V05);
	 h_isoIsov_fsr_Dr01 -> Fill(sumPt10G2/e2_reco.Pt(), iso_e2-iso_e2V10);
	 h_isoIsov_ePt_Dr008 -> Fill(e2_reco.Pt(), iso_e2-iso_e2V08);
	 h_isoIsov_ePt_Dr005 -> Fill(e2_reco.Pt(), iso_e2-iso_e2V05);
	 h_isoIsov_ePt_Dr01 -> Fill(e2_reco.Pt(), iso_e2-iso_e2V10);
	 
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
	 h_sumPtGammaDr13->Fill(sumPt13G1);
	 h_sumPtGammaDr13->Fill(sumPt13G2);
	 h_sumPtGammaDr15->Fill(sumPt15G1);
	 h_sumPtGammaDr15->Fill(sumPt15G2);
	 h_sumPtGammaDr17->Fill(sumPt17G1);
	 h_sumPtGammaDr17->Fill(sumPt17G2);
	 h_sumPtGammaDr20->Fill(sumPt20G1);
	 h_sumPtGammaDr20->Fill(sumPt20G2);
	 h_sumPtGammaDr25->Fill(sumPt25G1);
	 h_sumPtGammaDr25->Fill(sumPt25G2);
	 h_sumPtGammaDr30->Fill(sumPt30G1);
	 h_sumPtGammaDr30->Fill(sumPt30G2);
	 
	 
	 for (vector<double>::const_iterator itG = cosT_g1.begin(); itG != cosT_g1.end(); itG++){
	 h_cosT_EG -> Fill(*itG);
	 h_angleEnergy -> Fill(*itG, sumPtG1);
	 if (sumPtG1<cat1) h_cosT_EG1->Fill(*itG);
	 else if (sumPtG1< cat2) h_cosT_EG2->Fill(*itG);
	 else h_cosT_EG3->Fill(*itG);
	 }
	 
	 for (vector<double>::const_iterator itG = cosT_g2.begin(); itG != cosT_g2.end(); itG++){
	 h_cosT_EG -> Fill(*itG);
	 h_angleEnergy -> Fill(*itG, sumPtG2);
	 if (sumPtG2<cat1) h_cosT_EG1->Fill(*itG);
	 else if (sumPtG2< cat2) h_cosT_EG2->Fill(*itG);
	 else h_cosT_EG3->Fill(*itG);
	 }


	 if (sumPtG1<cat1){
     	   		 
	    h_ptRecoG1->Fill(e1_reco.Pt()); 
	    h_ptRecoPtGenG1->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenG1_Iso->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1), iso_e1);
	    h_ptRecoPtGenVsPtGenG1->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	
	 } else if (sumPtG1< cat2){
     			 
	    h_ptRecoG2->Fill(e1_reco.Pt());
	    h_ptRecoPtGenG2->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenVsPtGenG2->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenG2_Iso->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1), iso_e1);

	 
	 }else {
	 
	    h_ptRecoG3->Fill(e1_reco.Pt());
	    h_ptRecoPtGenG3->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenVsPtGenG3->Fill(e1_gen.Pt()+sumPtG1, e1_reco.Pt()-(e1_gen.Pt()+sumPtG1));
	    h_ptRecoPtGenG3_Iso->Fill(e1_reco.Pt()-(e1_gen.Pt()+sumPtG1), iso_e1);

	 }

	 if (sumPtG2<cat1){	    
	    h_ptRecoG1->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG1->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG1->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenG1_Iso->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2), iso_e2);

	 } else if (sumPtG2< cat2){
	    h_ptRecoG2->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG2->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG2->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenG2_Iso->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2), iso_e2);

	 }else {	    
	    h_ptRecoG3->Fill(e2_reco.Pt());
	    h_ptRecoPtGenG3->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenVsPtGenG3->Fill(e2_gen.Pt()+sumPtG2, e2_reco.Pt()-(e2_gen.Pt()+sumPtG2));
	    h_ptRecoPtGenG3_Iso->Fill(e2_reco.Pt()-(e2_gen.Pt()+sumPtG2), iso_e2);

	 }
         
      // analysis of the distance between the main cluster and the other ones
      const reco::CaloCluster * mainCluster1 = e1recoit->superClusterRef()->seed().get();
      const reco::CaloCluster * mainCluster2 = e2recoit->superClusterRef()->seed().get();
      int nOfMain = 0, nOfClusters=0;
      double maxDistClu=-1, maxdEtaClu=-1, maxdPhiClu=-1, etaClu7=-999, maxdX=-1, maxdY=-1;
      // for the first electron
      for (reco::CaloCluster_iterator itClu = e1recoit->superClusterRef()->clustersBegin(); 
	   itClu!= e1recoit->superClusterRef()->clustersEnd(); itClu++){
	 nOfClusters++;
	 if ( (&(*itClu))->get() != mainCluster1){
	    
	    double dPhiClu = fabs ((*itClu)->phi() - mainCluster1->phi());
	    if (dPhiClu > acos(-1)) dPhiClu = 2*acos(-1) - dPhiClu;
	    double dEtaClu = fabs((*itClu)->eta() - mainCluster1->eta());
	    double dRClu = sqrt( dPhiClu*dPhiClu  + dEtaClu*dEtaClu );
	    double deltaX = fabs ((*itClu)->x() - mainCluster1->x());
	    double deltaY = fabs ((*itClu)->y() - mainCluster1->y());

	    if (fabs(mainCluster1->eta())<1.479) {
	       h_deltaPhiCluEB->Fill(dPhiClu);
	       h_deltaEtaCluEB->Fill(dEtaClu);
	       h_deltaRCluEB->Fill(dRClu);
	       h_dPhidEtaCluEB->Fill(dPhiClu,dEtaClu);
	    } else {
	       h_deltaXCluEE->Fill(deltaX);
	       h_deltaYCluEE->Fill(deltaY);
	       h_deltaPhiCluEE->Fill(dPhiClu);
	       h_deltaEtaCluEE->Fill(dEtaClu);
	       h_deltaRCluEE->Fill(dRClu);
	       h_dPhidEtaCluEE->Fill(dPhiClu,dEtaClu);
	    }

	    if (dEtaClu>0.07){
	       h_ptCluVsPtRecoEta7->Fill( (*itClu)->energy(),e1recoit->pt());
	       h_ptCluVsFsrEta7->Fill( (*itClu)->energy(),sumPtG1);
	       etaClu7= (*itClu)->eta();
	       if (e1recoit->gsfElectronRef().isNull()) {h_isGsfRef->Fill(0);}
	       else {
		  h_isGsfRef->Fill(1);
		  if (e1recoit->gsfElectronRef()->trackerDrivenSeed()==1) h_isTrackSeed->Fill(0);
		  if (e1recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(1);
		  if (e1recoit->gsfElectronRef()->trackerDrivenSeed()==1 && 
		      e1recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(2);
		  if (e1recoit->gsfElectronRef()->trackerDrivenSeed()==0 && 
		      e1recoit->gsfElectronRef()->ecalDrivenSeed()==0) h_isTrackSeed->Fill(3);
		  if (e1recoit->gsfElectronRef()->trackerDrivenSeed()==1 && 
		      e1recoit->gsfElectronRef()->ecalDrivenSeed()==0) h_isTrackSeed->Fill(4);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==0 && 
		      e2recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(5);
	       }
	    }
	    if (dRClu > maxDistClu){
	       maxDistClu = dRClu;
	       maxdEtaClu = dEtaClu;
	       maxdPhiClu = dPhiClu;
	       maxdX = deltaX;
	       maxdY = deltaY;
	    }	    
	 } else {nOfMain++;}
      }      
      if (etaClu7!=-999) h_etaCluVsNCluEta7->Fill(etaClu7,nOfClusters);
      if (fabs(mainCluster1->eta())<1.479) {
	 h_nOfClustersEB->Fill(nOfClusters);	 
	 if (nOfClusters>1){
	    h_deltaPhiCluMaxEB->Fill(maxdPhiClu);
	    h_deltaEtaCluMaxEB->Fill(maxdEtaClu);
	    h_deltaRCluMaxEB->Fill(maxDistClu);
	    h_dPhidEtaCluMaxEB->Fill(maxdPhiClu,maxdEtaClu);
	 }
      } else {
	 h_nOfClustersEE->Fill(nOfClusters);
	 if (nOfClusters>1){
	    h_deltaXCluMaxEE->Fill(maxdX);
	    h_deltaYCluMaxEE->Fill(maxdY);
	    h_deltaPhiCluMaxEE->Fill(maxdPhiClu);
	    h_deltaEtaCluMaxEE->Fill(maxdEtaClu);
	    h_deltaRCluMaxEE->Fill(maxDistClu);
	    h_dPhidEtaCluMaxEE->Fill(maxdPhiClu,maxdEtaClu);
	 }
      }
      //----------------------------------------------------
      // for the second electron
      nOfClusters=0;
      maxDistClu=-1;
      etaClu7=-999;
      for (reco::CaloCluster_iterator itClu = e2recoit->superClusterRef()->clustersBegin(); 
	   itClu!= e2recoit->superClusterRef()->clustersEnd(); itClu++){
	 nOfClusters++;
	 if ( (&(*itClu))->get() != mainCluster2){
	    
	    double dPhiClu = fabs ((*itClu)->phi() - mainCluster2->phi());
	    if (dPhiClu > acos(-1)) dPhiClu = 2*acos(-1) - dPhiClu;
	    double dEtaClu = fabs((*itClu)->eta() - mainCluster2->eta());
	    double dRClu = sqrt( dPhiClu*dPhiClu  + dEtaClu*dEtaClu );
	    double deltaX = fabs ((*itClu)->x() - mainCluster2->x());
	    double deltaY = fabs ((*itClu)->y() - mainCluster2->y());
	    if (fabs(mainCluster2->eta())<1.479) {
	       h_deltaPhiCluEB->Fill(dPhiClu);
	       h_deltaEtaCluEB->Fill(dEtaClu);
	       h_deltaRCluEB->Fill(dRClu);
	       h_dPhidEtaCluEB->Fill(dPhiClu,dEtaClu);
	    } else {
	       h_deltaXCluEE->Fill(deltaX);
	       h_deltaYCluEE->Fill(deltaY);
	       h_deltaPhiCluEE->Fill(dPhiClu);
	       h_deltaEtaCluEE->Fill(dEtaClu);
	       h_deltaRCluEE->Fill(dRClu);
	       h_dPhidEtaCluEE->Fill(dPhiClu,dEtaClu);
	    } 

	    if (dEtaClu>0.07){
	       h_ptCluVsPtRecoEta7->Fill( (*itClu)->energy(),e2recoit->pt());
	       h_ptCluVsFsrEta7->Fill( (*itClu)->energy(),sumPtG2);
	       etaClu7=(*itClu)->eta();
	       if (e2recoit->gsfElectronRef().isNull()) {h_isGsfRef->Fill(0);}
	       else {
		  h_isGsfRef->Fill(1);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==1) h_isTrackSeed->Fill(0);
		  if (e2recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(1);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==1 && 
		      e2recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(2);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==0 && 
		      e2recoit->gsfElectronRef()->ecalDrivenSeed()==0) h_isTrackSeed->Fill(3);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==1 && 
		      e2recoit->gsfElectronRef()->ecalDrivenSeed()==0) h_isTrackSeed->Fill(4);
		  if (e2recoit->gsfElectronRef()->trackerDrivenSeed()==0 && 
		      e2recoit->gsfElectronRef()->ecalDrivenSeed()==1) h_isTrackSeed->Fill(5);
	       }
	    }

	    if (dRClu > maxDistClu){
	       maxDistClu = dRClu;
	       maxdEtaClu = dEtaClu;
	       maxdPhiClu = dPhiClu;
	       maxdX = deltaX;
	       maxdY = deltaY;
	    }	    
	 } else {nOfMain++;}
      } 
      if (etaClu7!=-999) h_etaCluVsNCluEta7->Fill(etaClu7,nOfClusters);
      if (fabs(mainCluster2->eta())<1.479) {
	 h_nOfClustersEB->Fill(nOfClusters);
	 if (nOfClusters>1){
	    h_deltaPhiCluMaxEB->Fill(maxdPhiClu);
	    h_deltaEtaCluMaxEB->Fill(maxdEtaClu);
	    h_deltaRCluMaxEB->Fill(maxDistClu);
	    h_dPhidEtaCluMaxEB->Fill(maxdPhiClu,maxdEtaClu);
	 }    
      } else {
	 h_nOfClustersEE->Fill(nOfClusters);
	 if (nOfClusters>1){
	    h_deltaXCluMaxEE->Fill(maxdX);
	    h_deltaYCluMaxEE->Fill(maxdY);
	    h_deltaPhiCluMaxEE->Fill(maxdPhiClu);
	    h_deltaEtaCluMaxEE->Fill(maxdEtaClu);
	    h_deltaRCluMaxEE->Fill(maxDistClu);
	    h_dPhidEtaCluMaxEE->Fill(maxdPhiClu,maxdEtaClu);
	 }
      }
 
      //cout <<"il numero di main cluster e' = " << nOfMain << endl;
      //----------------------------------------------------
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
