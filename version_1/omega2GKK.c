#include <iostream>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
typedef std::vector<int> Vint;

void omega2GKK(){
	TChain *t1 = new TChain("TreeAna");
	t1->Add("omegaRecoil_09_v1.root");
	TClonesArray*    Gamma = new TClonesArray("TLorentzVector");
	TClonesArray*    Pip = new TClonesArray("TLorentzVector");
	TClonesArray*    Pim = new TClonesArray("TLorentzVector");
	TClonesArray*    Kp = new TClonesArray("TLorentzVector");
	TClonesArray*    Km = new TClonesArray("TLorentzVector"); 
	TClonesArray*    Electronp = new TClonesArray("TLorentzVector");
	TClonesArray*    Electronm = new TClonesArray("TLorentzVector");
	TClonesArray*    Muonp = new TClonesArray("TLorentzVector");
	TClonesArray*    Muonm = new TClonesArray("TLorentzVector");
	TClonesArray*    Protonp = new TClonesArray("TLorentzVector");
	TClonesArray*    Protonm = new TClonesArray("TLorentzVector"); 	
	Int_t nGamma,nGood,nElectronp,nElectronm,nMuonp,nMuonm,nPip,nPim,nKp,nKm,nProtonp,nProtonm;
	Int_t runid,evtid;
	t1->SetBranchAddress("runid",&runid);
	t1->SetBranchAddress("evtid",&evtid);
	t1->SetBranchAddress("Gamma",&Gamma);
	t1->SetBranchAddress("Pip",&Pip);
	t1->SetBranchAddress("Pim",&Pim);
	t1->SetBranchAddress("Km",&Km);
	t1->SetBranchAddress("Kp",&Kp);	
	t1->SetBranchAddress("Electronp",&Electronm);
	t1->SetBranchAddress("Electronm",&Electronm);
	t1->SetBranchAddress("Muonm",&Muonm);
	t1->SetBranchAddress("Muonp",&Muonp);
	t1->SetBranchAddress("nProtonp",&nProtonp);
	t1->SetBranchAddress("nProtonm",&nProtonm);
	t1->SetBranchAddress("nGamma", &nGamma);
	t1->SetBranchAddress("nGood", &nGood);
	t1->SetBranchAddress("nElectronp",&nElectronp);
	t1->SetBranchAddress("nElectronm",&nElectronm);
	t1->SetBranchAddress("nMuonp",&nMuonp);
	t1->SetBranchAddress("nMuonm",&nMuonm);
	t1->SetBranchAddress("nPip",&nPip);
	t1->SetBranchAddress("nPim",&nPim);
	t1->SetBranchAddress("nKp",&nKp);
	t1->SetBranchAddress("nKm",&nKm);
	t1->SetBranchAddress("nProtonp",&nProtonp);
	t1->SetBranchAddress("nProtonm",&nProtonm);

	Long64_t nevt = t1->GetEntries();
	int count=0;
	for(Long64_t k=0;k<nevt;k++){
		t1->GetEntry(k);
		if(nGamma!=4||nPim!=1||nPip!=1||nMuonm!=0||nMuonp!=0||nKm!=1||nKp!=1||nElectronm!=0||nElectronp!=0||nProtonm!=0|| nProtonp!=0) continue;
		TLorentzVector gamma(0,0,0,0),electronp(0,0,0,0),electronm(0,0,0,0),muonp(0,0,0,0),muonm(0,0,0,0),pip(0,0,0,0),pim(0,0,0,0),kp(0,0,0,0),km(0,0,0,0),protonp(0,0,0,0),protonm(0,0,0,0),X(0,0,0,0);
		for(int ig=0;ig<nGamma;ig++){
			gamma += *(TLorentzVector*)Gamma->At(ig);
		}
		if(gamma.M()>0) X=gamma;

		for(int ielectronp=0;ielectronp<nElectronp;ielectronp++){
			electronp += *(TLorentzVector*)Electronp->At(ielectronp);
		}
		if(electronp.M()>0) X += electronp;

		for(int ielectronm=0;ielectronm<nElectronm;ielectronm++){
			electronm += *(TLorentzVector*)Electronm->At(ielectronm);
		}
		if(electronm.M()>0) X += electronm;

		for(int imuonp=0;imuonp<nMuonp;imuonp++){
			muonp += *(TLorentzVector*)Muonp->At(imuonp);
		}
		if(muonp.M()>0) X += muonp;

		for(int imuonm=0;imuonm<nMuonm;imuonm++){
			muonm += *(TLorentzVector*)Muonm->At(imuonm);
		}
		if(muonm.M()>0) X += muonm;
			
		for(int ipip=0;ipip<nPip;ipip++){
			pip += *(TLorentzVector*)Pip->At(ipip);
		}
		if(pip.M()>0) X += pip;

		for(int ipim=0;ipim < nPim;ipim++){
			pim += *(TLorentzVector*)Pim->At(ipim);
		}
		if(pim.M()>0) X += pim;

		for(int ikp=0;ikp<nKp;ikp++){
			kp += *(TLorentzVector*)Kp->At(ikp);
		}
		if(kp.M()>0) X += kp;

		for(int ikm=0;ikm<nKm;ikm++){
			km += *(TLorentzVector*)Km->At(ikm);	
		}
		if(km.M()>0) X += km;
	
		for(int iprotonp=0;iprotonp<nProtonp;iprotonp++){
			protonp += *(TLorentzVector*)Protonp->At(iprotonp);
		}
		if(protonp.M()>0) X += protonp;

		for(int iprotonm=0;iprotonm<nProtonm;iprotonm++){
			protonm += *(TLorentzVector*)Protonm->At(iprotonm);
		}
		if(protonm.M()>0) X += protonm;	 

		if(X.M()>2.9&&X.M()<3.2){
		cout<<"  runid: "<<runid<<"           "<<"   evtid: "<<evtid<<endl;
		cout<<"  Gamma: "<<(*(TLorentzVector*)Gamma->At(0)).Px()<<"     "<<(*(TLorentzVector*)Gamma->At(0)).E()<<"      "<<endl;
		cout<<"      "<<(*(TLorentzVector*)Gamma->At(1)).Px()<<"     "<<(*(TLorentzVector*)Gamma->At(1)).E()<<"      "<<endl;;
		cout<<"      "<<(*(TLorentzVector*)Gamma->At(2)).Px()<<"     "<<(*(TLorentzVector*)Gamma->At(2)).E()<<"      "<<endl;
		cout<<"      "<<(*(TLorentzVector*)Gamma->At(3)).Px()<<"     "<<(*(TLorentzVector*)Gamma->At(3)).E()<<"      "<<endl;
		cout<<"   Kp: "<<(*(TLorentzVector*)Kp->At(0)).Px()<<"     "<<(*(TLorentzVector*)Kp->At(0)).E()<<"      "<<endl;
		cout<<"   Km: "<<(*(TLorentzVector*)Km->At(0)).Px()<<"     "<<(*(TLorentzVector*)Km->At(0)).E()<<"      "<<endl;
		cout<<"   Pip: "<<(*(TLorentzVector*)Pip->At(0)).Px()<<"     "<<(*(TLorentzVector*)Pip->At(0)).E()<<"      "<<endl;
		cout<<"   Pim: "<<(*(TLorentzVector*)Pim->At(0)).Px()<<"     "<<(*(TLorentzVector*)Pim->At(0)).E()<<"      "<<endl;
		count++;
		}
		if(count==33) break;
	}
}

