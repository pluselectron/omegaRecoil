#ifndef Physics_Analysis_omegaRecoil_H
#define Physics_Analysis_omegaRecoil_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"    //No NTuple!
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/IJobOptionsSvc.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EmcRecEventModel/RecEmcHit.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "EmcRawEvent/EmcDigi.h"
#include "EvTimeEvent/RecEsTime.h"
#include "EventNavigator/EventNavigator.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "MdcRecEvent/RecMdcKalTrack.h" 
#include "DstEvent/TofHitStatus.h"
#include "RootCnvSvc/RootInterface.h"
#include "RootCnvSvc/RootCnvSvc.h"
#include "ParticleID/ParticleID.h"
#include "PartPropSvc/PartPropSvc.h"

//#include "VertexFit/ReadBeamParFromDb.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
 
#include "Identifier/Identifier.h"
#include "McTruth/McParticle.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/MucMcHit.h"  
#include "McTruth/McEvent.h"

#include "TMath.h"
#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include "TBenchmark.h"

using CLHEP::Hep3Vector;
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<TLorentzVector> Vt4;

class omegaRecoil : public Algorithm {

	public:
		omegaRecoil(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  

	private:
		bool Primary_Vertex(Hep3Vector*,HepSymMatrix*);
        void Photon_Select(SmartDataPtr<EvtRecEvent>,SmartDataPtr<EvtRecTrackCol> ,Hep3Vector,Vint *);
        void Assign_Phomomen(SmartDataPtr<EvtRecTrackCol>,Int_t,Vint,Hep3Vector,Vp4 *);
        void Charged_Select(SmartDataPtr<EvtRecEvent>,SmartDataPtr<EvtRecTrackCol>,Hep3Vector,Vint *,Int_t *);
        bool Particle_Id(SmartDataPtr<EvtRecTrackCol>,Vint,Int_t,int *);
        void Assign_Trkmomen(SmartDataPtr<EvtRecTrackCol>,Vint,Int_t,int *, Vp4 *);
        void Vtxfit(SmartDataPtr<EvtRecTrackCol>,Vint,Int_t,Vint,Int_t ,Vint *,Vint *)
        void Reconstruct_2gamma(SmartDataPtr<EvtRecTrackCol> ,Vint ,Int_t ,Vp4 ,double ,Vint *,Vint *);
      
        //Declare r0, z0 and cos cut for charged tracks
                Double_t m_vr0cut;
                Double_t m_vz0cut;
                Double_t m_ccoscut;
                //Declare energy, dphi, dthe cuts for fake gamma's
                Double_t m_isoAngleCut;		
		//Declare flag for analysis
                Int_t m_readSig;                
                //Declare name of output file and data structure
                TFile *saveFile;
				std::string m_OutputFileName;
                TClonesArray *Gamma;
                TClonesArray *Electronp;
                TClonesArray *Electronm;
                TClonesArray *Muonp;
                TClonesArray *Muonm;
                TClonesArray *Pip;
                TClonesArray *Pim;
                TClonesArray *Kp;
                TClonesArray *Km;
                TClonesArray *Protonp;
                TClonesArray *Protonm;

                //Define Tree structure 
				TTree *TreeAna;
                Int_t runid;
                Int_t evtid;
                //Int_t nevt;
                Int_t nGamma;
                Int_t nGood;
                Int_t nPip;
                Int_t nPim;
                Int_t nKp;
                Int_t nKm;
                Int_t nProtonp;
                Int_t nProtonm;
                Int_t nElectronp;
                Int_t nElectronm;
                Int_t nMuonp;
                Int_t nMuonm;
                Int_t igamma1_omega;
                Int_t igamma2_omega;
                Int_t ipip_omega;
                Int_t ipim_omega;
                Double_t m_omega;
                
};

#endif 


