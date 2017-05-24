#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "omegaRecoilAlg/omegaRecoil.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TRandom.h"
const double mpi = 0.13957;
const double mpi0 = 0.1349766; 
const double momega = 0.78265;
const double mk = 0.493677;
const double meta = 0.547853;
const double me = 0.000510998910;
const double mmu = 0.105658367;
const double mp = 0.93827203; 
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double Ejpsi=3.097;
const double Epsip=3.686;
const double Econt=3.650;
const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<TLorentzVector> Vt4;
using namespace std;
/////////////////////////////////////////////////////////////////////////////
omegaRecoil::omegaRecoil(const std::string& name, ISvcLocator* pSvcLocator) :		//do not change this function
	Algorithm(name, pSvcLocator) {
		//Declare the properties  
		declareProperty("OutputFileName",  m_OutputFileName = "out.root");
		declareProperty("Vr0cut", m_vr0cut=2.0);
		declareProperty("Vz0cut", m_vz0cut=10.0);
		declareProperty("Ccoscut", m_ccoscut=0.93);
		declareProperty("IsoAngleCut", m_isoAngleCut=10.0);
		declareProperty("ReadSig", m_readSig = 1);
	}
/////////////////////////////////////////////////////////////////////////////
StatusCode omegaRecoil::initialize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
	TString s_OutputFileName(m_OutputFileName);
	s_OutputFileName.ReplaceAll("[\"","");
	s_OutputFileName.ReplaceAll("\"]","");
	saveFile = new TFile(s_OutputFileName, "recreate");
    //*********** data structure ******************
    saveFile = new TFile(s_OutputFileName, "recreate");
	Gamma = new TClonesArray("TLorentzVector"); 
	Pip = new TClonesArray("TLorentzVector");
	Pim = new TClonesArray("TLorentzVector"); 
	Kp = new TClonesArray("TLorentzVector");
	Km = new TClonesArray("TLorentzVector");
	Protonp = new TClonesArray("TLorentzVector");
	Protonm = new TClonesArray("TLorentzVector");
	Electronp = new TClonesArray("TLorentzVector");
	Electronm = new TClonesArray("TLorentzVector");
	Muonp = new TClonesArray("TLorentzVector");
	Muonm = new TClonesArray("TLorentzVector");
	//*********** tree structure ******************
	if(m_readSig==1){
		TreeAna = new TTree("TreeAna", "analysis");
		TreeAna->Branch("runid", &runid, "runid/I");
		TreeAna->Branch("evtid", &evtid, "evtid/I");
		TreeAna->Branch("Gamma", "TClonesArray", &Gamma, 256000, 0);
		TreeAna->Branch("Electronp","TClonesArray", &Electronp, 256000,0);
		TreeAna->Branch("Electronm","TClonesArray", &Electronm, 256000,0);
		TreeAna->Branch("Muonp","TClonesArray", &Muonp, 256000,0);
		TreeAna->Branch("Muonm","TClonesArray", &Muonm, 256000,0);  
		TreeAna->Branch("Pip","TClonesArray", &Pip, 256000,0);
		TreeAna->Branch("Pim","TClonesArray", &Pim, 256000,0);
		TreeAna->Branch("Kp","TClonesArray", &Kp, 256000,0);
		TreeAna->Branch("Km","TClonesArray", &Km, 256000,0);
		TreeAna->Branch("Protonp","TClonesArray", &Protonp, 256000,0);
		TreeAna->Branch("Protonm","TClonesArray", &Protonm, 256000,0);
		TreeAna->Branch("nGamma", &nGamma, "nGamma/I");
		TreeAna->Branch("nGood", &nGood, "nGood/I");
		TreeAna->Branch("nElectronp",&nElectronp,"nElectronp/I");
		TreeAna->Branch("nElectronm",&nElectronm,"nElectronm/I");
		TreeAna->Branch("nMuonp",&nMuonp,"nMuonp/I");
		TreeAna->Branch("nMuonm",&nMuonm,"nMuonm/I");
		TreeAna->Branch("nPip",&nPip,"nPip/I");
		TreeAna->Branch("nPim",&nPim,"nPim/I");
		TreeAna->Branch("nKp",&nKp,"nKp/I");
		TreeAna->Branch("nKm",&nKm,"nKm/I");
		TreeAna->Branch("nProtonp",&nProtonp,"nProtonp/I");
		TreeAna->Branch("nProtonm",&nProtonm,"nProtonm/I");
		TreeAna->Branch("m_omega", &m_omega, "m_omega/D");
		TreeAna->Branch("igamma1_omega",&igamma1_omega,"igamma1_omega/I");
		TreeAna->Branch("igamma2_omega",&igamma2_omega,"igamma2_omega/I");
		TreeAna->Branch("ipip_omega",&ipip_omega,"ipip_omega/I");
		TreeAna->Branch("ipim_omega",&ipim_omega,"ipim_omega/I");
		//Initialize TClonesArray
		Gamma->BypassStreamer();  
		Pip->BypassStreamer();
		Pim->BypassStreamer();
		Kp->BypassStreamer();
		Km->BypassStreamer();
		Protonp->BypassStreamer();
		Protonm->BypassStreamer();
		Electronp->BypassStreamer();
		Electronm->BypassStreamer();
		Muonp->BypassStreamer();
		Muonm->BypassStreamer();
	}//endof if
	//--------end of book--------
	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////
StatusCode omegaRecoil::execute() {
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");  
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	log << MSG::DEBUG <<"run, evtnum = "
		<< runNo << " , "
		<< event <<endreq;  
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::INFO << "get event tag OK" << endreq;
	log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;
	double ecms = 3.097;
	HepLorentzVector cms(0.011*ecms,0,0,ecms);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	if(event%1000 == 0){
		cout<<"Processing "<<event<<"th event..."<<endl;
	}
	//*****************Global Parameter*****************
	runid = runNo;
	evtid=event;
	m_omega=0;
	nElectronp=0;nElectronm=0;nMuonp=0;nMuonm=0;nPip=0;nPim=0;nKp=0;nKm=0;nProtonp=0;nProtonm=0;
	CLHEP::Hep3Vector xorigin(0,0,0);
	///*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*// 
	//*****************Primary Vertex*****************
	CLHEP::HepSymMatrix VtxErr(3,0);
	if(!Primary_Vertex(&xorigin,&VtxErr)) return StatusCode::SUCCESS; 
	//************************************************

	//***********Good Photon selection***************
	nGamma=0;
	Vint iGamma;
	iGamma.clear();
	Photon_Select(evtRecEvent,evtRecTrkCol,xorigin,&iGamma);
	nGamma = iGamma.size();
	log << MSG::DEBUG << "num Good Photon " << nGamma  << " , " <<evtRecEvent->totalNeutral()<<endreq;
	if(nGamma<2) return StatusCode::SUCCESS;
	//************************************************

	//**********Assign 4-momentum to each photon******
	Vp4 pGamma;
	pGamma.clear();
	Assign_Phomomen(evtRecTrkCol,nGamma,iGamma,xorigin,&pGamma);
	//************************************************

	//***********Good Charged Track Selection*********
	nGood=0;
	Vint iGood;
	Int_t nCharge=0;
	iGood.clear();
	Charged_Select(evtRecEvent,evtRecTrkCol,xorigin,&iGood,&nCharge);
	nGood = iGood.size();
	log << MSG::DEBUG << "num Good Charged Track "  << nGood << " , " <<"num Charged "<<nCharge<<endreq;
	if((nGood<2)||(nCharge!=0)) return StatusCode::SUCCESS;
	//************************************************

	//***********Particle Identify (PID)**************
	int idparticle[20]={-1}; 
	if(!Particle_Id(evtRecTrkCol,iGood,nGood,idparticle)) return StatusCode::SUCCESS;
	//************************************************

	//*****Assign 4-momentum to each Charged Track****
	Vp4 pAtrk;
	pAtrk.clear();
	Assign_Trkmomen(evtRecTrkCol,iGood,nGood,idparticle,&pAtrk);
	//************************************************

	//* Build all pi0 group by pi0-->gamma + gamma * 
	Vint gamma1_pi0, gamma2_pi0;
	double Rec_m = mpi0;
	gamma1_pi0.clear();
	gamma2_pi0.clear();
	Reconstruct_2gamma(evtRecTrkCol,iGamma,nGamma,pGamma,Rec_m,&gamma1_pi0,&gamma2_pi0);
	int nPi0=gamma1_pi0.size();
	if(nPi0<1) return StatusCode::SUCCESS;
	//***********************************
cout<<" NUMBER OF #pi0: "<<nPi0<<endl;
	//*********FIND Pi+ and Pi-***************
	Vp4 pPip,pPim;
    Vint iPip,iPim;
	pPip.clear();
	pPim.clear();
    iPip.clear();
    iPim.clear();
	for(int j=0;j<nGood;j++){
        if(idparticle[j]==4) {pPip.push_back(pAtrk[j]);iPip.push_back(iGood[j]);}
        if(idparticle[j]==5) {pPim.push_back(pAtrk[j]);iPim.push_back(iGood[j]);}
	} 
	int npip = pPip.size(); int npim = pPim.size();
	if(npip<1||npim<1) return StatusCode::SUCCESS;
cout<<" number of #pi +: "<<npip<<"   "<<"number of #pi -: "<<npim<<endl;
for(int i=0;i<npip;i++){cout<<"ID of PI+: "<<iPip[i]<<"  ";}
for(int i=0;i<npim;i++){cout<<"iD of PI-: "<<iPim[i]<<"  ";}
	//************PI+ PI- vertex fit***********************
    Vint viPip，viPim;
    viPip.clear();
    viPim.clear();
    Vtxfit(evtRecTrkCol,iPip,npip,iPim,npim,&viPip,&viPim);
    int n_vtxpippim = viPip.size();
    if(n_vtxpippim<1) return StatusCode::SUCCESS;
    //***********************************
cout<<"AFTER PI+ PI- VERTEX FIT THE NUMBER OF PI+PI-: "<<n_vtxpippim<<endl;
for(int i=0;i<n_vtxpippim;i++){cout<<" after vetex the ID of Pi+: "<<viPip[i]<<"  ";}
for(int i=0;i<n_vtxpippim;i++){cout<<" After vetex the ID of Pi-: "<<viPim[i]<<"  ";}
	//**Build the best Omega  by Omega-->pi0 + pi+ + pi- **
	Vint pip_omega,pim_omega,gamma1_omega,gamma2_omega;
	vector<double> deta_omega;    
	pip_omega.clear();
	pim_omega.clear();
	gamma1_omega.clear();
	gamma2_omega.clear();
	deta_omega.clear();
	for(int i0=0;i0<nPi0;i0++){
		for(int i1=0;i1<n_vtxpippim;i1++){
				HepLorentzVector omega = pGamma[gamma1_pi0[i0]]+pGamma[gamma2_pi0[i0]]+pPip[viPip[i1]]+pPim[viPim[i1]];
				if(fabs(omega.m()-momega)<0.030){
					pip_omega.push_back(viPip[i1]);
					pim_omega.push_back(viPim[i1]);
					gamma1_omega.push_back(gamma1_pi0[i0]);
					gamma2_omega.push_back(gamma2_pi0[i0]);
					deta_omega.push_back(pow(fabs(omega.m()-momega),2));
			}
		}
	}
cout<<" satisfy mass window of Omega "<<endl;
for(int i=0;i<pip_omega.size();i++){cout<<" satisfy Omega'S Pi+: "<<pip_omega[i]<<"  ";}
for(int i=0;i<pim_omega.size();i++){cout<<" Satisfy Omega'S Pi-: "<<pim_omega[i]<<"  ";}
	int nOmega = deta_omega.size();    
	if(nOmega<1) return StatusCode::SUCCESS;
	igamma1_omega=-1;igamma2_omega=-1;ipip_omega=-1;ipim_omega=-1;
	int min=-1;
	for(int i=0;i<nOmega;i++){
		if(deta_omega[i]>min){
			min=deta_omega[i]; 
			ipip_omega=pip_omega[i];
			ipim_omega=pim_omega[i];
			igamma1_omega=gamma1_omega[i];
			igamma2_omega=gamma2_omega[i];
		}        
	}
	//**************save in TClonsArray***************
	int n0=0,n1=0,n2=0,n3=0,n4=0,n5=0,n6=0,n7=0,n8=0,n9=0,n10=0; //e+e-u+u-pi+pi-K+K-P+P-Gamma
	int ntotal; ntotal=nGood+nGamma;
	for (int i = 0; i < nGood; ++i)
	{
		HepLorentzVector ptrk = pAtrk[i];
		TLorentzVector p_Particle;
		p_Particle.SetPxPyPzE(ptrk.px(),ptrk.py(),ptrk.pz(),ptrk.e());
		if(idparticle[i]==0) {new ((*Electronp)[n0]) TLorentzVector(p_Particle);n0++;}
		if(idparticle[i]==1) {new ((*Electronm)[n1]) TLorentzVector(p_Particle);n1++;}
		if(idparticle[i]==2) {new ((*Muonp)[n2]) TLorentzVector(p_Particle);n2++;}   
		if(idparticle[i]==3) {new ((*Muonm)[n3]) TLorentzVector(p_Particle);n3++;}   
		if(idparticle[i]==4) {new ((*Pip)[n4]) TLorentzVector(p_Particle);n4++;}
		if(idparticle[i]==5) {new ((*Pim)[n5]) TLorentzVector(p_Particle);n5++;} 
		if(idparticle[i]==6) {new ((*Kp)[n6]) TLorentzVector(p_Particle);n6++;}  
		if(idparticle[i]==7) {new ((*Km)[n7]) TLorentzVector(p_Particle);n7++;}  
		if(idparticle[i]==8) {new ((*Protonp)[n8]) TLorentzVector(p_Particle);n8++;} 
		if(idparticle[i]==9) {new ((*Protonm)[n9]) TLorentzVector(p_Particle);n9++;}
	}
	for (int i = 0; i < nGamma; ++i)
	{
		HepLorentzVector ptrk = pGamma[i];
		TLorentzVector p_Particle;
		p_Particle.SetPxPyPzE(ptrk.px(),ptrk.py(),ptrk.pz(),ptrk.e());
		new ((*Gamma)[n10]) TLorentzVector(p_Particle);n10++;
	}
	nElectronp=n0;nElectronm=n1;nMuonp=n2;nMuonm=n3;nPip=n4;nPim=n5;nKp=n6;nKm=n7;nProtonp=n8;nProtonm=n9;
	m_omega=(pGamma[igamma1_omega]+pGamma[igamma2_omega]+pPip[ipip_omega]+pPim[ipim_omega]).m();
	TreeAna->Fill();

	cout<<"nelectronp:"<<n0<<" "<<"nelectronm:"<<n1<<" "<<"nmounp:"<<n2<<" "<<"nmounm:"<<n3<<" "<<"npip:"<<n4<<" "
		<<"npim:"<<n5<<" "<<"nkp:"<<n6<<" "<<"nkm:"<<n7<<" "<<"nprotonp:"<<n8<<" "<<"nprotonm:"<<n9<<" "<<"nGamma:"<<n10<<endl;
	cout<<" M_pip: "<<(*((TLorentzVector*)Pip->At(0))).M()<<endl;
	cout<<" FINISH EXECUTE "<<"  "<<" M_omega: "<<m_omega<<" "<<" RUNID: "<<runid<<" " <<" EVTID: "<<evtid<<endl;
	return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode omegaRecoil::finalize() {
	cout << "In finalize()..." << endl;
	saveFile->cd();
	if(m_readSig==1) TreeAna->Write();
	saveFile->Close();
	return StatusCode::SUCCESS;
}
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/Member Function*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*//
//*****************Primary Vertex***********************************************************************
bool omegaRecoil::Primary_Vertex(Hep3Vector *xorigin,HepSymMatrix *VtxErr){
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex(); 
		double*  vv = vtxsvc->SigmaPrimaryVertex();  
		(*xorigin).setX(dbv[0]);
		(*xorigin).setY(dbv[1]);
		(*xorigin).setZ(dbv[2]);                
		(*VtxErr)[0][0] = vv[0]*vv[0];
		(*VtxErr)[1][1] = vv[1]*vv[1];
		(*VtxErr)[2][2] = vv[2]*vv[2];
		return true;
	}else return false; 
}
//***************** Good Photon selection***************************************************************
void omegaRecoil::Photon_Select(SmartDataPtr<EvtRecEvent> evtRecEvent,SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Hep3Vector xorigin,Vint *iGamma){
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++){
		if(i >= evtRecTrkCol->size()) break;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		// find the nearest charged track  
		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.; 
		for(int j = 0; j < evtRecEvent->totalCharged(); j++){
			if(j >= evtRecTrkCol->size()) break;
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition() - xorigin;
			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);
			thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			if(angd < dang){
				dang = angd;
				dthe = thed;
				dphi = phid;
			}
		}
		if(dang >= 200) continue;
		double eraw = emcTrk->energy();
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);
		// good photon cut will be set here     
		double the = emcpos.theta();
		double e_threshold = 10.0;
		if(fabs(cos(the)) < 0.8)   e_threshold = 0.025;
		else if((fabs(cos(the)) > 0.86) && (fabs(cos(the)) < 0.92)) e_threshold = 0.050;
		if(eraw < e_threshold) continue;
		if(fabs(dang) < m_isoAngleCut) continue;
		if(emcTrk->time()>14  || emcTrk->time()<0) continue;
		(*iGamma).push_back(i);        
	}
}
//*****************Assign 4-momentum to each photon*****************************************
void omegaRecoil::Assign_Phomomen(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Int_t nGamma,Vint iGamma,Hep3Vector xorigin,Vp4 *pGamma){
	for(int i = 0; i < nGamma; i++){
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGamma[i]; 
		RecEmcShower* emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z()); 
		Hep3Vector gammaDirection = emcpos - xorigin;
		double phi = gammaDirection.phi();
		double the = gammaDirection.theta();
		//if(fabs(gammaDirection.cosTheta())>0.93)continue;
		HepLorentzVector ptrk;
		ptrk.setPx(eraw*sin(the)*cos(phi));
		ptrk.setPy(eraw*sin(the)*sin(phi));
		ptrk.setPz(eraw*cos(the));
		ptrk.setE(eraw);
		(*pGamma).push_back(ptrk);
	}
}
//***************** Good Charged Track Selection*****************************************
void omegaRecoil::Charged_Select(SmartDataPtr<EvtRecEvent> evtRecEvent,SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Hep3Vector xorigin,Vint *iGood,Int_t *nCharge){
	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		if(i >= evtRecTrkCol->size()) break;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch=mdcTrk->p();
		double ptch=mdcTrk->pxy();
		double x0=mdcTrk->x();
		double y0=mdcTrk->y();
		double z0=mdcTrk->z();
		double phi0=mdcTrk->helix(1);
		double xv=xorigin.x();
		double yv=xorigin.y();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];
		double ccos=cos(mdcTrk->theta());
		if(fabs(ccos) > m_ccoscut) continue;
		if(fabs(Rvz0) >= m_vz0cut ) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
		if(mdcKalTrk->charge() == 0) continue;
		(*nCharge) +=mdcTrk->charge();
		(*iGood).push_back(i);
	}
}
// *********************************PID*****************************************
bool omegaRecoil::Particle_Id(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Vint iGood,Int_t nGood,int idparticle[]){
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		pid->init();
		pid->setMethod(pid->methodProbability());//particle species given by a probability density function(P(x;p,H))
		pid->setChiMinCut(4);//for dedx sys,tof1 tof2 sys,emc sys,muc sys etc all ChiMinCut ??
		pid->setRecTrack(*itTrk);//for all sys RecTrack ??
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useEmc() | pid->useMuc()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyMuon() | pid->onlyProton() | pid->onlyElectron());    // seperater Pion/Kaon/Muon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;
		Double_t pro_particle[5];
		pro_particle[0]  = pid->probElectron();
		pro_particle[1]  = pid->probMuon();
		pro_particle[2]  = pid->probPion();
		pro_particle[3]  = pid->probKaon();
		pro_particle[4]  = pid->probProton();
		double max=0;
		for(int j=0;j<5;j++){
			if(pro_particle[j]>max) {
				max=pro_particle[j];
				if(mdcKalTrk->charge() > 0){idparticle[i]=2*j;}
				else{idparticle[i]=2*j+1;}
			}
		}
	}
	return true;
}
// ************* Assign 4-momentum to each charged track*****************************************
void omegaRecoil::Assign_Trkmomen(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Vint iGood,Int_t nGood,int idparticle[], Vp4 *pAtrk){
	for (int i = 0; i < nGood; ++i)
	{   
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		RecMdcKalTrack::setPidType  (RecMdcKalTrack::PidType(int (idparticle[i]/2)));
		HepLorentzVector ptrk;
		ptrk.setPx(mdcKalTrk->px());
		ptrk.setPy(mdcKalTrk->py());
		ptrk.setPz(mdcKalTrk->pz());
		double p3 = ptrk.mag();
		ptrk.setE(sqrt(p3*p3+mpi*mpi));
		(*pAtrk).push_back(ptrk);
	}               
}
// ************* vertex fit  *************
void omegaRecoil::Vtxfit(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Vint iPip,Int_t npip,Vint iPim,Int_t npim,Vint *viPip,Vint *viPim){
        //Test vertex fit
    HepPoint3D vx(0., 0., 0.);
    HepSymMatrix Evx(3, 0);
    double bx = 1E+6;
    double by = 1E+6;
    double bz = 1E+6;
    Evx[0][0] = bx*bx;
    Evx[1][1] = by*by;
    Evx[2][2] = bz*bz;
    VertexParameter vxpar;
    vxpar.setVx(vx);
    vxpar.setEvx(Evx);
    VertexFit* vtxfit = VertexFit::instance();
    vtxfit->init();
    Double_t vxchisq;
    WTrackParameter wvallTrk[2];
    for (int i=0; i<npip; i++) {
        for (int j=0;j<npim;j++) {
            RecMdcKalTrack *PipTrk = (*(evtRecTrkCol->begin()+iPip[i]))->mdcKalTrack();
            RecMdcKalTrack *PimTrk = (*(evtRecTrkCol->begin()+iPim[j]))->mdcKalTrack();
            PipTrk->setPidType(RecMdcKalTrack::PidType(2));
            PimTrk->setPidType(RecMdcKalTrack::PidType(2));
            wvallTrk[0] = WTrackParameter(xmass[2], itTrk->getZHelix(), itTrk->getZError());
            wvallTrk[1] = WTrackParameter(xmass[2], itTrk->getZHelix(), itTrk->getZError());
            vtxfit->AddTrack(0,wvallTrk[0]);
            vtxfit->AddTrack(1,wvallTrk[1]);
            vtxfit->AddVertex(0,vxpar,0,1);
            if(!vtxfit->Fit(0)) continue;
            vxchisq = vtxfit->chisq(0);
            vtxfit->Swim(0);
            (*viPip).push_back(iPip[i]);
            (*viPim).push_back(iPim[j]);
        }
    }
}
//**********************************
void omegaRecoil::Reconstruct_2gamma(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,Vint iGamma,Int_t nGamma,Vp4 pGamma,double Rec_m,Vint *gamma1,Vint *gamma2){
	double mass_cut=0.02;
	vector<double> vkmchisq; //for find the least chisq
	vkmchisq.clear();
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	for(int i1=0;i1<nGamma-1;i1++){
		RecEmcShower* gammaTrk1 = (*(evtRecTrkCol->begin()+iGamma[i1]))->emcShower();
		for(int i2=i1+1;i2<nGamma;i2++){
			RecEmcShower* gammaTrk2 = (*(evtRecTrkCol->begin()+iGamma[i2]))->emcShower();
			if(fabs((pGamma[i1]+pGamma[i2]).m()-Rec_m) > mass_cut) continue;
			kmfit->init();
			kmfit->AddTrack(0, 0.0, gammaTrk1);
			kmfit->AddTrack(1, 0.0, gammaTrk2);
			kmfit->AddResonance(0, 0.1349766, 0, 1);
			if(!kmfit->Fit(0)) continue;
			(*gamma1).push_back(iGamma[i1]);
			(*gamma2).push_back(iGamma[i2]);
			vkmchisq.push_back(kmfit->chisq());
                //HepLorentzVector p_pi0 = pGamma[i1]+pGamma[i2]；
                //TLorentzVector pi0;
                //pi0.SetPxPyPzE(p_pi0.px(),p_pi0.py(),p_pi0.pz(),p_pi0.e());
                //new ((*Pi0)[n]) TLorentzVector(pi0);n++;
            
		}
	}
}




