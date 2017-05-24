#include <iostream>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
typedef std::vector<int> Vint;

void classify_plot(){
	TStopwatch timer;
	TChain *t1 = new TChain("TreeAna");
	t1->Add("omegaRecoil_09_v2.root");

	Int_t nGamma,nGood,nElectronp,nElectronm,nMuonp,nMuonm,nPip,nPim,nKp,nKm,nProtonp,nProtonm;
	Int_t igamma1_omega,igamma2_omega,ipip_omega,ipim_omega;
	Double_t m_omega;
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
	t1->SetBranchAddress("igamma1_omega",&igamma1_omega);
	t1->SetBranchAddress("igamma2_omega",&igamma2_omega);
	t1->SetBranchAddress("ipip_omega",&ipip_omega);
	t1->SetBranchAddress("ipim_omega",&ipim_omega);
	t1->SetBranchAddress("m_omega",&m_omega);

	//Histogram
	TH1D *hrecoil_omega[40960];TH1D *hrecoil_rank[40960];
	TH1D *hX[40960];TH1D *hX_Rank[40960];
	for(int i=0;i<40960;i++){
		ostringstream mystring, myString;
		myString << "hrecoil_omega"<< i;
		mystring <<"hX"<<i;
		hrecoil_omega[i] = new TH1D(myString.str().c_str(),"",(2.0-1.0)/0.005,1.0,2.0); 	  	
		hX[i] = new TH1D(mystring.str().c_str(),"",(2.0-1.0)/0.005,1.0,2.0);
		hrecoil_omega[i]->GetXaxis()->SetTitle("recoil_omega");
		hrecoil_omega[i]->GetYaxis()->SetTitle("Events/5Mev");
		hX[i]->GetXaxis()->SetTitle("X invariant  Mass");
		hX[i]->GetYaxis()->SetTitle("Events/5Mev");
	}

	//***********recoil_omega and X********
	TLorentzVector gamma1_omega(0,0,0,0),gamma2_omega(0,0,0,0),pip_omega(0,0,0,0),pim_omega(0,0,0,0);
	TLorentzVector omega(0,0,0,0),cms(0.011*3.097,0,0,3.097),Jpsi(0,0,0,0);
	Double_t recoil_omega,rest_invariMass;
	Int_t count[40960]={0};
	Long64_t i_hist;
	Long64_t nevt = t1->GetEntries();
	for(Long64_t k=0;k<nevt;k++){
		t1->GetEntry(k);
		if(nGamma>10||nPim>4||nPip>4||nMuonm>2||nMuonp>2||nKm>2||nKp>2||nElectronm>2||nElectronp>2||nProtonm>2|| nProtonp>2) continue;
		i_hist=(nGamma-2)*4096 + (nPim-1)*1024 + (nPip-1)*256 + nMuonm*128 + nMuonp*64 + nKm*32 + nKp*16 + nElectronm*8 + nElectronp*4 + nProtonm*2 + nProtonp;		
		count[i_hist] +=1;

		ostringstream histTitle;
		histTitle<<"nG:"<<nGamma-2<<"-nEp:"<<nElectronp<<"-Em:"<<nElectronm<<"-Mup:"<<nMuonp<<"-Mum:"<<nMuonm<<"-nPip:"<<nPip-1<<"-nPim:"<<nPim-1<<"-nKp:"<<nKp<<"-nKm:"<<nKm<<"-nPp:"<<nProtonp<<"-nProtonm:"<<nProtonm;
		if(count[i_hist]==100){
		hrecoil_omega[i_hist]->SetTitle(histTitle.str().c_str());
		hX[i_hist]->SetTitle(histTitle.str().c_str());
		//cout<<histTitle.str().c_str()<<endl;	
		}

		TLorentzVector gamma(0,0,0,0),electronp(0,0,0,0),electronm(0,0,0,0),muonp(0,0,0,0),muonm(0,0,0,0),pip(0,0,0,0),pim(0,0,0,0),kp(0,0,0,0),km(0,0,0,0),protonp(0,0,0,0),protonm(0,0,0,0),X(0,0,0,0);
		for(int ig=0;ig<nGamma;ig++){
			if(ig!=igamma1_omega&&ig!=igamma2_omega) gamma += *(TLorentzVector*)Gamma->At(ig);
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
			if(ipip!=ipip_omega) pip += *(TLorentzVector*)Pip->At(ipip);
		}
		if(pip.M()>0) X += pip;

		for(int ipim=0;ipim < nPim;ipim++){
			if(ipim!=ipim_omega) pim += *(TLorentzVector*)Pim->At(ipim);
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

		gamma1_omega = (TLorentzVector*)Gamma->At(igamma1_omega);
		gamma2_omega = (TLorentzVector*)Gamma->At(igamma2_omega);
		pip_omega = (TLorentzVector*)Pip->At(ipip_omega);
		pim_omega = (TLorentzVector*)Pim->At(ipim_omega);
		omega = gamma1_omega + gamma2_omega + pip_omega + pim_omega;
	
		if(omega.M()<0) continue;
		recoil_omega = (cms-omega).M();
		Jpsi = X + omega;
    	if(Jpsi.M()>2.9&&Jpsi.M()<3.2){
    		if(recoil_omega>0) hrecoil_omega[i_hist]->Fill(recoil_omega);
    		if(X.M()>0) hX[i_hist]->Fill(X.M());
    	}				
	}
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
			if(hX[i2]->GetEntries() >= hX[i2-1]->GetEntries()){
				h2 = hX[i2-1];
				hX[i2-1] = hX[i2];
				hX[i2] = h2;
			}
		}
		//cout<<" i "<<i1<<" H_RECOIL_OMEGA ENTRIES "<<hrecoil_omega[i1]->GetEntries()<<endl;
	}
	cout<<" AFTER RANK BY ENTRIES "<<endl;

	//Draw
	TCanvas *crecoil_omega[150]; TCanvas *cX[150];
	int c1=0,c2=0;
	for(int i=0;i<150;i++){
		ostringstream myCanvas1,myCanvas2,cX_save,crecoil_save;
		myCanvas1<<"Crecoil_"<<i;
		crecoil_save<<"./recoil_omega/Crecoil_"<<i<<".png";
		myCanvas2<<"CX_"<<i;
		cX_save<<"./Classify_plot/CX_"<<i<<".png";		
		crecoil_omega[i] = new TCanvas(myCanvas1.str().c_str(),myCanvas1.str().c_str(),1200,800);
		cX[i] = new TCanvas(myCanvas2.str().c_str(),myCanvas2.str().c_str(),1200,800);
		crecoil_omega[i]->Divide(2,2);
		cX[i]->Divide(2,2);
		for(int j=0;j<4;j++){ 
			crecoil_omega[i]->cd(j+1);
			if(hrecoil_omega[c1]->GetEntries()>0) hrecoil_omega[c1++]->Draw();
			cX[i]->cd(j+1);
			if(hX[c2]->GetEntries()>0) hX[c2++]->Draw();else break;
		}

		crecoil_omega[i]->SaveAs(crecoil_save.str().c_str());
		cX[i]->SaveAs(cX_save.str().c_str());
	}

	double cputime = timer.CpuTime();
	printf("RT=%7.3f s, Cpu=%7.3f s\n",timer.RealTime(),cputime);



}










