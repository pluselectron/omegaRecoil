#include <iostream>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
typedef std::vector<int> Vint;

TLorentzVector addallMomentun(int,TClonesArray* );

void omegarecoil_plot(){
	TStopwatch timer;
	TChain *t1 = new TChain("TreeAna");
	t1->Add("version_2/line01_11.root");
	t1->Add("version_2/line12_22.root");
	t1->Add("version_2/line12_22_1.root");
	t1->Add("version_2/line23_33.root");
	t1->Add("version_2/line23_33_1.root");
	t1->Add("version_2/line34_43.root");
	t1->Add("version_2/completeRoot2.root");
	t1->Add("version_2/completeRoot3.root");

	Int_t nGamma,nGood,nElectronp,nElectronm,nMuonp,nMuonm,nPip,nPim,nKp,nKm,nProtonp,nProtonm;
	Int_t igamma1_omega,igamma2_omega,ipip_omega,ipim_omega;
	Int_t evtid,runid;
	TClonesArray*    Gamma = new TClonesArray("TLorentzVector");
	TClonesArray*    Electronp = new TClonesArray("TLorentzVector");
	TClonesArray*    Electronm = new TClonesArray("TLorentzVector");
	TClonesArray*    Muonp = new TClonesArray("TLorentzVector");
	TClonesArray*    Muonm = new TClonesArray("TLorentzVector");
	TClonesArray*    Pip = new TClonesArray("TLorentzVector");
	TClonesArray*    Pim = new TClonesArray("TLorentzVector");
	TClonesArray*    Kp = new TClonesArray("TLorentzVector");
	TClonesArray*    Km = new TClonesArray("TLorentzVector");  	
	TClonesArray*    Protonp = new TClonesArray("TLorentzVector");
	TClonesArray*    Protonm = new TClonesArray("TLorentzVector");
	TClonesArray*	 Omega = new TClonesArray("TLorentzVector");

	t1->SetBranchAddress("Gamma",&Gamma);
	t1->SetBranchAddress("Electronp",&Electronm);
	t1->SetBranchAddress("Electronm",&Electronm);
	t1->SetBranchAddress("Muonm",&Muonm);
	t1->SetBranchAddress("Muonp",&Muonp);
	t1->SetBranchAddress("Pip",&Pip);
	t1->SetBranchAddress("Pim",&Pim);
	t1->SetBranchAddress("Km",&Km);
	t1->SetBranchAddress("Kp",&Kp);	
	t1->SetBranchAddress("Protonm",&Protonm);
	t1->SetBranchAddress("Protonp",&Protonm);
	t1->SetBranchAddress("Omega",&Omega);
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
	t1->SetBranchAddress("runid",&runid);
	t1->SetBranchAddress("evtid",&evtid);
	t1->SetBranchAddress("igamma1_omega",&igamma1_omega);
	t1->SetBranchAddress("igamma2_omega",&igamma2_omega);
	t1->SetBranchAddress("ipip_omega",&ipip_omega);
	t1->SetBranchAddress("ipim_omega",&ipim_omega);

	//Histogram
	TH1D *hrecoil_omega[40960];TH1D *hrecoil_rank[40960];
	for(int i=0;i<40960;i++){
		ostringstream myString;
		myString << "hrecoil_omega"<< i;
		hrecoil_omega[i] = new TH1D(myString.str().c_str(),"",(2.0-1.0)/0.005,1.0,2.0); 	  	
		hrecoil_omega[i]->GetXaxis()->SetTitle("recoil_omega");
		hrecoil_omega[i]->GetYaxis()->SetTitle("Events/5Mev");
		hrecoil_omega[i]->SetTitleSize(0.3);
	}

	//***********recoil_omega and X********
	TLorentzVector gamma1_omega(0,0,0,0),gamma2_omega(0,0,0,0),pip_omega(0,0,0,0),pim_omega(0,0,0,0);
	TLorentzVector omega(0,0,0,0),cms(0.011*3.097,0,0,3.097);
	Double_t recoil_omega;
	Int_t count[40960]={0};
	Long64_t i_hist;

	Long64_t nevt = t1->GetEntries();
//	for(Long64_t k=0;k<nevt;k++){
		t1->GetEntry(15);
		if(nGamma>10||nPim>4||nPip>4||nMuonm>2||nMuonp>2||nKm>2||nKp>2||nElectronm>2||nElectronp>2||nProtonm>2|| nProtonp>2) continue;
		i_hist=(nGamma-2)*4096 + (nPim-1)*1024 + (nPip-1)*256 + nMuonm*128 + nMuonp*64 + nKm*32 + nKp*16 + nElectronm*8 + nElectronp*4 + nProtonm*2 + nProtonp;		
		count[i_hist] +=1;

		ostringstream histTitle;
		histTitle<<"nG:"<<nGamma-2<<" - e+:"<<nElectronp<<" - e-:"<<nElectronm<<" - #mu+:"<<nMuonp<<" - #mu-:"<<nMuonm<<" - #pi+:"<<nPip-1<<" - #pi-:"<<nPim-1<<" - K+:"<<nKp<<" - K-:"<<nKm<<" - p+:"<<nProtonp<<" - p-:"<<nProtonm;
		if(count[i_hist]==100){
		hrecoil_omega[i_hist]->SetTitle(histTitle.str().c_str());
		//cout<<histTitle.str().c_str()<<endl;	nElectronp,Electronp
		}

		omega = (TLorentzVector*)Omega->At(0);
		if(omega.M()<0) continue;
cout<<" 111  "<<"   15  "<<endl;
/*		X = addallMomentun(nGamma,Gamma) + addallMomentun() + addallMomentun(nElectronm,Electronm) + addallMomentun(nMuonp,Muonp) + addallMomentun(nMuonm,Muonm) + 
		addallMomentun(nPip,Pip) + addallMomentun(nPim,Pim) + addallMomentun(nKp,Kp) + addallMomentun(nKm,Km) + addallMomentun(nProtonp,Protonp) + addallMomentun(nProtonm,Protonm);
*/
cout<<"  222  "<<endl;
//delete addallMomentun;
//delete Gamma;
//delete Electronp;delete Electronm;delete Muonp;delete Muonm;delete Pip;delete Pim; delete Kp;delete Km;delete Protonp;delete Protonm;
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

		
    	//if(X.M()>2.9&&X.M()<3.2){
    		recoil_omega = (cms-omega).M();
    		if(recoil_omega>0) hrecoil_omega[i_hist]->Fill(recoil_omega);
//cout<<" fill   "<<recoil_omega<<endl;
    	//}
cout<<" 333 "<<endl;				
//	}
	cout<<" FINISH FILL "<<endl;

	//Rank
	TH1D *h1;	
	TH1D *h2;
	for(int i1=0;i1<40960;i1++){
		for(int i2=40960-1;i2>i1;i2--){
			if(hrecoil_omega[i2]->GetEntries() >= hrecoil_omega[i2-1]->GetEntries()){
				h1=hrecoil_omega[i2-1];
				hrecoil_omega[i2-1]=hrecoil_omega[i2];
				hrecoil_omega[i2]=h1;
			}
		}
		//cout<<" i "<<i1<<" H_RECOIL_OMEGA ENTRIES "<<hrecoil_omega[i1]->GetEntries()<<endl;
	}
	cout<<" AFTER RANK BY ENTRIES "<<endl;

	//Draw
	TCanvas *crecoil_omega[10]; 
	int c1=0,c2=0;
	for(int i=0;i<10;i++){
		ostringstream myCanvas1,crecoil_save;
		myCanvas1<<"Crecoil_"<<i;
		crecoil_save<<"./recoil_omega/Crecoil_"<<i<<".pdf";
		crecoil_omega[i] = new TCanvas(myCanvas1.str().c_str(),myCanvas1.str().c_str(),1200,800);
		crecoil_omega[i]->Divide(2,2);
		for(int j=0;j<4;j++){ 
			crecoil_omega[i]->cd(j+1);
			if(hrecoil_omega[c1]->GetEntries()>0) hrecoil_omega[c1++]->Draw();
			else break;
		}
		crecoil_omega[i]->SaveAs(crecoil_save.str().c_str());
	}

	double cputime = timer.CpuTime();
	printf("RT=%7.3f s, Cpu=%7.3f s\n",timer.RealTime(),cputime);



}

/*TLorentzVector addallMomentun(int N,TClonesArray* Particle){
		TLorentzVector P(0,0,0,0),Q(0,0,0,0);
		for(int i=0;i<N;i++){
			P += *(TLorentzVector*)Particle->At(i);
//cout<<" i "<<P.M()<<endl;
		//Particle->BypassStreamer();;
		}
	if(P.M()>0) return P;
	else return		Q;
	
} 
*/








