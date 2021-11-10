/* PiPiPi0TagsAlg | 02/11/2021 | Author: Dr. Alex Gilman (Many functions written by Dr. Hajime Muramatsu)
 * This code identifies D0D0 pairs where one D0 decays to pi+pi-pi0 and the other decays to a selected number of final states.
 * A single candidate is chosen based on the the invariant mass recoiling against the D^{*+} (recoil mass or mRec).
 * PDG codes are from https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
 *
 * What's saved in the ntuple for each event: 
 * 
 */

//Dependencies

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"
#include "EvTimeEvent/RecEsTime.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecEtaToGG.h"

//-----------------------------------------
// MC stuff
#include "GaudiKernel/IPartPropSvc.h"

// EventNavigator
#include "EventNavigator/EventNavigator.h"

// MC data
#include "McTruth/McParticle.h"

// MDC reconstructed data
#include "MdcRecEvent/RecMdcKalTrack.h"

// EMC reconstructed data
#include "EmcRecEventModel/RecEmcShower.h"

// Muon
#include "MucRecEvent/RecMucTrack.h"
 //-----------------------------------------

// Math and Tupling
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "CLHEP/Geometry/Point3D.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D < double > HepPoint3D;
#endif

// Vertex (and K-FITS)
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/WTrackParameter.h"
#include "ParticleID/ParticleID.h"
#include "MucRecEvent/RecMucTrack.h"
#include "VertexFit/KalmanKinematicFit.h"

// PID used in DTag
#include "SimplePIDSvc/ISimplePIDSvc.h"

#include "MeasuredEcmsSvc/IMeasuredEcmsSvc.h" // <<==== Run-dependent Ebeam for 4180 data
 //-----------------------------------------
#include "AIDA/AIDA.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IAxis.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/IHistogramFactory.h"

using AIDA::IHistogram1D;
using AIDA::IHistogram2D;
using AIDA::IHistogram3D;
using AIDA::ICloud1D;
//std::auto_ptr<AIDA::IAnalysisFactory>  af( AIDA_createAnalysisFactory() );
//std::auto_ptr<AIDA::ITreeFactory> tf( af -> createTreeFactory() ); 
//std::auto_ptr<AIDA::ITree> tree( tf -> create() ); 
//std::auto_ptr<AIDA::IHistogramFactory> hf( af->createHistogramFactory( *tree ) ); 
//-----------------------------------------
#include <iostream>
#include <fstream>

using namespace std;
//-------------------------------------------

// needed to convert the Cell ID
#include "Identifier/Identifier.h"

#include "/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics/PiPiPi0TagsAlg/00-00-01/PiPiPi0TagsAlg/PiPiPi0Tags.h"
#include <vector>

// Define Constants
const double m_photonmass = 0.;
const double m_electmass = 0.000510998910;
const double m_muonmass = 0.105658367;
const double m_pionmass = 0.13957018;
const double m_pi0mass = 0.1349766;
const double m_kaonmass = 0.493677;
const double m_ksmass = 0.497614;
const double m_etamass = 0.547853;
const double m_omgmass = 0.78265;
const double m_protonmass = 0.938272013;
const double m_phimass = 1.019460;
const double m_D0mass = 1.86480;
const double m_dpmass = 1.86961;
const double m_dsmass = 1.96830;
const double m_dstzero = 2.00696;
const double m_dstplus = 2.01026;
const double m_dsstmass = 2.1121;
const double m_jpsimass = 3.096916;

const unsigned int maxNumPhotons=10;
const unsigned int maxNumFSR=5;
// const double velc = 299.792458;   // tof path unit in mm
typedef std::vector < int > Vint;
typedef std::vector < HepLorentzVector > Vp4;

// counter
static long m_cout_all(0);

/////////////////////////////////////////////////////////////////////////////

PiPiPi0Tags::PiPiPi0Tags(const std::string & name, ISvcLocator * pSvcLocator):
  Algorithm(name, pSvcLocator) {
    //Declare the properties  

  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PiPiPi0Tags::initialize() {
  cout << "Initialize" << endl;
  MsgStream logmess(msgSvc(), name());

  logmess << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;

  NTuplePtr nt1(ntupleSvc(), "FILE1/ups");
  if (nt1)
    m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc() -> book("FILE1/ups", CLID_ColumnWiseTuple, "My first example");

    if (m_tuple1) {
      status = m_tuple1 -> addItem("EVT", m_EVT);
      status = m_tuple1 -> addItem("run", m_run);
      status = m_tuple1 -> addItem("Ecm", m_Ecm);
      status = m_tuple1 -> addItem("m_betavec", 3,m_betavec);

      status = m_tuple1 -> addItem("DStar_gamma",    m_DStar_gamma);
      status = m_tuple1 -> addItem("DStar_pi0",      m_DStar_pi0);
      status = m_tuple1 -> addItem("DStarBar_gamma", m_DStarBar_gamma);
      status = m_tuple1 -> addItem("DStarBar_pi0",   m_DStarBar_pi0);

      status = m_tuple1 -> addItem("D0_KK_gen",       m_D0_KK_gen);
      status = m_tuple1 -> addItem("D0Bar_KK_gen",    m_D0Bar_KK_gen);
      status = m_tuple1 -> addItem("D0_KSPi0_gen",       m_D0_KSPi0_gen);
      status = m_tuple1 -> addItem("D0Bar_KSPi0_gen",    m_D0Bar_KSPi0_gen);
      status = m_tuple1 -> addItem("D0_PiPiPi0_gen",       m_D0_PiPiPi0_gen);
      status = m_tuple1 -> addItem("D0Bar_PiPiPi0_gen",    m_D0Bar_PiPiPi0_gen);

      status = m_tuple1 -> addItem("D0_KK_rec",       m_D0_KK_rec);
      status = m_tuple1 -> addItem("D0Bar_KK_rec",    m_D0Bar_KK_rec);
      status = m_tuple1 -> addItem("D0_KSPi0_rec",       m_D0_KSPi0_rec);
      status = m_tuple1 -> addItem("D0Bar_KSPi0_rec",    m_D0Bar_KSPi0_rec);
      status = m_tuple1 -> addItem("D0_PiPiPi0_rec",       m_D0_PiPiPi0_rec);
      status = m_tuple1 -> addItem("D0Bar_PiPiPi0_rec",    m_D0Bar_PiPiPi0_rec);

      status = m_tuple1 -> addItem("D2_KSPi0_mInv",    m_D2_KSPi0_mInv);
      status = m_tuple1 -> addItem("D2_KSPi0_mRec",    m_D2_KSPi0_mRec);
      status = m_tuple1 -> addItem("D1_PiPiPi0_mInv",    m_D1_PiPiPi0_mInv);
      status = m_tuple1 -> addItem("D1_PiPiPi0_mRec",    m_D1_PiPiPi0_mRec);
      status = m_tuple1 -> addItem("D2_PiPiPi0_mInv",    m_D2_PiPiPi0_mInv);
      status = m_tuple1 -> addItem("D2_PiPiPi0_mRec",    m_D2_PiPiPi0_mRec);
      status = m_tuple1 -> addItem("D2_KSPi0_mInv",    m_D2_KSPi0_mInv);

      status = m_tuple1 -> addItem("D1_SubPiPi_mInv",    m_D1_SubPiPi_mInv);// Invariant mass of the Pi+Pi- pair in PiPiPi0/KSPi0
      status = m_tuple1 -> addItem("D1_pi0_mInv",       m_D1_pi0_mInv);// Invariant mass of the pi0 pair in PiPiPi0/KSPi0
      status = m_tuple1 -> addItem("D2_SubPiPi_mInv",    m_D2_SubPiPi_mInv);//Invariant mass of the Pi+Pi- pair in PiPiPi0/KSPi0
      status = m_tuple1 -> addItem("D2_pi0_mInv",    m_D2_pi0_mInv);

      status = m_tuple1 -> addItem("Emiss",  m_Emiss);// Invariant mass of the pi0 pair in PiPiPi0/KSPi0
      status = m_tuple1 -> addItem("mmiss2", m_mmiss2);


    } else {
      logmess << MSG::ERROR << "    Cannot book N-nTuple:" << long(m_tuple1) << endmsg;

      return StatusCode::FAILURE;
    }

  }

  logmess << MSG::INFO << "successfully return from initialize()" << endmsg;

  return StatusCode::SUCCESS;
  // **********************************************************************

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PiPiPi0Tags::execute() {

  m_myevt++;

  StatusCode sc = StatusCode::SUCCESS;
  //save the events passed selection to a new file
  setFilterPassed(false);

  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < Event::EventHeader > eventHeader(eventSvc(), "/Event/EventHeader");
  SmartDataPtr<EvtRecPi0Col> recPi0Col(eventSvc(), "/Event/EvtRec/EvtRecPi0Col");
  int runNo = eventHeader -> runNumber(); //Record run number (negative for MC events)
  int eventNo = eventHeader -> eventNumber(); // Record event number

  m_run = runNo;
  m_EVT = eventNo;

  m_ebeam = 4.17861 / 2; // Default Beam Energies
  if (fabs(runNo) >= 47543 && fabs(runNo) <= 48170) {
    m_ebeam = 4.18877 / 2.0;
  } //4190
  else if (fabs(runNo) >= 48172 && fabs(runNo) <= 48713) {
    m_ebeam = 4.19890 / 2.0;
  } //4200
  else if (fabs(runNo) >= 48714 && fabs(runNo) <= 49239) {
    m_ebeam = 4.20921 / 2.0;
  } //4210
  else if (fabs(runNo) >= 49270 && fabs(runNo) <= 49787) {
    m_ebeam = 4.21874 / 2.0;
  } //4220
  else if (fabs(runNo) >= 32239 && fabs(runNo) <= 32849) {
    m_ebeam = (4320.34 - 0.00287 * fabs(runNo)) * 0.5 / 1000.;
  } //4230
  else if (fabs(runNo) >= 32850 && fabs(runNo) <= 33484) {
    m_ebeam = 4.22554 / 2.0;
  } //4230

  Hep3Vector m_beta; //Default Beta Vector
  m_beta.setX(0.011);
  m_beta.setY(0);
  m_beta.setZ(0);

  Hep3Vector xorigin(0, 0, 0); // Read in the Interaction Point vertex from the database
  IVertexDbSvc * vtxsvc;
  HepSymMatrix xorigin_err(3,0);
  Gaudi::svcLocator() -> service("VertexDbSvc", vtxsvc);
  if (vtxsvc -> isVertexValid()) {

    double * vertex = vtxsvc -> PrimaryVertex();
    double* vsigma = vtxsvc->SigmaPrimaryVertex();
    xorigin.setX(vertex[0]);
    xorigin.setY(vertex[1]);
    xorigin.setZ(vertex[2]);

    xorigin_err[0][0] = vsigma[0]*vsigma[0];
    xorigin_err[1][1] = vsigma[1]*vsigma[1];
    xorigin_err[2][2] = vsigma[2]*vsigma[2];
  }

  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < RecMucTrackCol > mucTracks(eventSvc(), "/Event/Recon/RecMucTrackCol");
  SmartDataPtr < EvtRecVeeVertexCol > evtRecVeeVertexCol(eventSvc(), "/Event/EvtRec/EvtRecVeeVertexCol");
  SmartDataPtr < EvtRecPi0Col > recPi0Col(eventSvc(), "/Event/EvtRec/EvtRecPi0Col");
  SmartDataPtr < EvtRecEtaToGGCol > recEtaCol(eventSvc(), "/Event/EvtRec/EvtRecEtaToGGCol");

  /////////////////////// Fill Defaults ////////////////////////////

  HepLorentzVector m_cm_p4 = 2 * m_ebeam * HepLorentzVector(0, 0, 0, 1); // in the cm frame

m_D0_kpi_gen=0;
m_D0Bar_kpi_gen=0;
m_D0_kspipi_gen=0;
m_D0Bar_kspipi_gen=0;
m_D0_kpi_rec=0;
m_D0Bar_kpi_rec=0;
m_D0_kspipi_rec=0;
m_D0Bar_kspipi_rec=0;

  ////////////////////////// Generator Truth Information ////////////////////////////////////
  /*
   * Loop over the generated decay tree to find the processes of interest: DST->Pi D0 and D0->KPi.
   */
  if (runNo < 0) {
 bool kaon_found=false;
 bool pion_found=false;
 int phot_index=0;
 int FSR_index=0;
    SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(), EventModel::MC::McParticleCol);
    if (mcParticleCol) {
      Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
        }
        if (( * iter_mc) -> particleProperty() == 421) //D0
        {
          HepLorentzVector D_p4=( * iter_mc)->initialFourMomentum();
      status = m_tuple1 -> addItem("D0_PiPiPi0_pip_px",       m_D0_PiPiPi0_pip_px);
      status = m_tuple1 -> addItem("D0_PiPiPi0_pip_py",       m_D0_PiPiPi0_pip_py);
      status = m_tuple1 -> addItem("D0_PiPiPi0_pip_pz",       m_D0_PiPiPi0_pip_pz);
      status = m_tuple1 -> addItem("D0_PiPiPi0_pim_px",       m_D0_PiPi0_pim_px);
      status = m_tuple1 -> addItem("D0_PiPiPi0_pim_py",       m_D0_PiPi0_pim_py);
      status = m_tuple1 -> addItem("D0_PiPiPi0_pim_pz",       m_D0_PiPi0_pim_pz);
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int km = 0; km < gc.size(); km++) {
              if (gc[km] -> particleProperty() == -321) //K-
              {
                for (unsigned int kp = 0; kp < gc.size(); kp++) {
                  if (gc[kp] -> particleProperty() == 321) //K+
                  {
                    m_D0_KK_gen = 1;
                  } //kp found
                } //kp loop
              } //km found
            } //km loop
          } // if 2 daughters
          
        } //if D0
        if (( * iter_mc) -> particleProperty() == -421) //D0Bar
        {
          HepLorentzVector D_p4=( * iter_mc)->initialFourMomentum();
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2 ) {
            for (unsigned int kp = 0; kp < gc.size(); kp++) {
              if (gc[kp] -> particleProperty() == 321) //K+
              {
                for (unsigned int kp = 0; km < gc.size(); km++) {
		  if (kp==km) {continue;}
                  if (gc[pi] -> particleProperty() == -211) //K-
                  { 
                    m_D0Bar_KK_gen = 1;
                  } //kp found
                } //kp loop
              } //km found
	      
            } //km loop
          } // if 2 daughters
        } //if D0Bar

        if (( * iter_mc) -> particleProperty() == 423) //DStar0
        {
          HepLorentzVector DStar_p4=( * iter_mc)->initialFourMomentum();
          m_DStarMass= DStar_p4.mag();
          m_DStarMom= DStar_p4.vect().mag();
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int k = 0; k < gc.size(); k++) {
              if (gc[k] -> particleProperty() == 22){m_DStar_gamma=1;} //Gamma
              if (gc[k] -> particleProperty() == 111){m_DStar_pi0=1;} //Pi0
            } //transition loop
          } // if 2 daughters
        } //if D0Star

        if (( * iter_mc) -> particleProperty() == -423) //DStar0Bar
        {
          HepLorentzVector DStar_p4=( * iter_mc)->initialFourMomentum();
          m_DStarMass= DStar_p4.mag();
          m_DStarMom= DStar_p4.vect().mag();
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int k = 0; k < gc.size(); k++) {
              if (gc[k] -> particleProperty() == 22){m_DStarBar_gamma=1;} //Gamma
              if (gc[k] -> particleProperty() == 111){m_DStarBar_pi0=1;} //Pi0
            } //transition loop
          } // if 2 daughters
        } //if D0StarBar

      } //generated particle loop
    } //mc particle col

     if(countPar( 421, 211)==1&& countPar( 421,-211)==1 &&
       countPar( 421, 321)==0 && countPar( 421,-321)==0 &&
       countPar( 421,111)==1&& countPar( 421, 221)==0&&
       countPar( 421,310)==0 && countPar( 421, 113)==0&&
       countPar( 310,22)==0 && countPar( 421,331)==0 ) {m_D0_PiPiPi0_gen=true;}

     if(countPar( -421, 211)==1&& countPar( -421,-211)==1 &&
       countPar( -421, 321)==0&& countPar( -421,-321)==0 &&
       countPar( -421,111)==1&& countPar( -421, 221)==0&&
       countPar( -421,310)==0 && countPar( -421, 113)==0&&
       countPar(310,22)==0 && countPar( -421,331)==0 ) {m_D0Bar_PiPiPi0_gen=true;}

     if(countPar( 421, 211)==1&& countPar( 421,-211)==1 &&
       countPar( 421, 321)==0 && countPar( 421,-321)==0 &&
       countPar( 421,111)==1&& countPar( 421, 221)==0&&
       countPar( 421,310)==1 && countPar( 421, 113)==0&&
       countPar( 310,22)==0 && countPar( 421,331)==0 ) {m_D0_Kspi0_gen=true;}

     if(countPar( -421, 211)==1&& countPar( -421,-211)==1 &&
       countPar( -421, 321)==0&& countPar( -421,-321)==0 &&
       countPar( -421,111)==1&& countPar( -421, 221)==0&&
       countPar( -421,310)==1 && countPar( -421, 113)==0&&
       countPar(310,22)==0 && countPar( -421,331)==0 ) {m_D0Bar_Kspi0_gen=true;}

  } //if MC

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Begin DTag

  DTagTool dtagTool;
  if (dtagTool.isDTagListEmpty()) {}
  dtagTool.setPID(true);
  DTagToolIterator iter_begin = dtagTool.modes_begin();
  DTagToolIterator iter_end = dtagTool.modes_end();
  vector < int > d0itindex = dtagTool.D0modes();

  double min_delta_DStar_Mass=1000;
  double min_delta_D_Mass=1000;

  for (int D1_index = 0; D1_index < d0itindex.size(); D1_index++) // Loop over D0 Candidates
  {
    DTagToolIterator D1_iter_dtag = dtagTool.modes_begin() + d0itindex[D1_index];
    m_ebeam= (*D1_iter_dtag)->beamE();
    m_Ecm= 2* m_ebeam;
    m_cm_p4 = 2 * m_ebeam * HepLorentzVector(0, 0, 0, 1); // in the cm frame
    int D1_mode = ( * D1_iter_dtag) -> decayMode();
    if (D1_mode != 101) {continue;} //must be PiPiPi0
    // Default for charged kaons and pions are just be "good tracks".
    //     // If PID is needed, require pdg=1 (=EvtRecDTag::Tight)
    if (( * D1_iter_dtag) -> type() != 1) {
      continue;
    } //Require PID be applied to Kaon/Pion daughters of the D candidate
    //cout << "PID" << endl;  
    if (( * D1_iter_dtag) -> type() != EvtRecDTag::Tight) {
      continue;
    }
    //        //cout << "EvtRecDTag" << endl;    

    HepLorentzVector D1_p4 = ( * D1_iter_dtag) -> p4();
    D1_p4.boost(-m_beta);

    SmartRefVector < EvtRecTrack > D1_tracks = ( * D1_iter_dtag) -> tracks(); //List of tracks used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > D1_showers = ( * D1_iter_dtag) -> showers(); //List of showers used in the D0 reconstruction
    Vint D1_pi0Index=dtagTool.pi0Id(D1_iter_dtag, 1);
    for(EvtRecPi0Col::iterator pi0Itr = recPi0Col->begin(); pi0Itr < recPi0Col->end(); pi0Itr++){
      int pi0idx=pi0Itr-recPi0Col->begin();
      if(pi0Itr - recPi0Col->begin() != D1_pi0Index[0])continue;
    }

    // D1 Candidate Particles
    RecMdcKalTrack * D1_pip_mdcTrk = (D1_tracks[0]) -> mdcKalTrack();
    D1_pip_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    RecMdcKalTrack * D1_pim_mdcTrk = (D1_tracks[1]) -> mdcKalTrack();
    D1_pim_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    EvtRecPi0Col::iterator D1_pi0= D1_pi0Index[0] + recPi0Col->begin();
    HepLorentzVector D1_pi0_p4 = D1_pi0->hiPfit() + D1_pi0->loPfit();
    RecEmcShower *  D1_hiEn_shw =  const_cast<EvtRecTrack*>((*D1_pi0)->hiEnGamma())->emcShower();
    RecEmcShower *  D1_loEn_shw =  const_cast<EvtRecTrack*>((*pi0)->loEnGamma())->emcShower();
    RecEmcShower *  D1_g1_shw = (D1_showers[0]) -> emcShower();
    RecEmcShower *  D1_g2_shw = (D1_showers[1]) -> emcShower();

    if(D1_hiEn_shw->trackID()!= D1_g1_shw->trackId() && D1_hiEn_shw->trackID()!= D1_g2_shw->trackId() )
    {
      cout<<"Uh-oh! You've got the wrong pi0 there, friend."<<endl;
    }
    // Truth Matching: Compare simulated detector hits with generated trajectories to match a track with a generated particle
    // The MCTPiPiPi0DCHG function can return the PDG code associated with a track or the code associated its parent. It can also be used to require that the PDG codes of the particle (2nd argument) and its parent (3rd argument) match certain values. See the MCTPiPiPi0DCHG function for more details.

    int this_D0_PiPiPi0_rec=0;
    if (runNo < 0) {
    }
  if (m_D0_PiPiPi0_gen
      && MCSHPIDCHG(D1_g1_shw -> trackId(), 22, 0, 421) &&  MCSHPIDCHG(D1_g1_shw -> trackId(), 22, 0, 421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D1_pip_mdcTrk -> trackId(), 211, 0, 421) && MCTKPIDCHG(D1_pim_mdcTrk -> trackId(), -211, 0, 421) // Are Pi+/- from D0?
      && !MCTKPIDCHG(D1_pip_mdcTrk -> trackId(), 211, 310, 421) && !MCTKPIDCHG(D1_pim_mdcTrk -> trackId(), -211, 310, 421)  ) // But not from K0S
     {m_D0_PiPiPi0_rec=1;
      this_D0_PiPiPi0_rec=1;
}
  if (m_D0Bar_PiPiPi0_gen
      && MCSHPIDCHG(D1_g1_shw -> trackId(), 22, 0, -421) &&  MCSHPIDCHG(D1_g1_shw -> trackId(), 22, 0, -421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D1_pip_mdcTrk -> trackId(), 211, 0, -421) && MCTKPIDCHG(D1_pim_mdcTrk -> trackId(), -211, 0, -421) // Are Pi+/- from D0Bar?
      && !MCTKPIDCHG(D1_pip_mdcTrk -> trackId(), 211, 310, -421) && !MCTKPIDCHG(D1_pim_mdcTrk -> trackId(), -211, 310, -421)  ) // But not from K0S
     {m_D0Bar_PiPiPi0_rec=1;
      this_D0Bar_PiPiPi0_rec=1;
     }

      HepLorentzVector D1_con_tag_mom = HepLorentzVector(D1_p4.vect(),pow(D1_p4.vect().dot(D1_p4.vect())+m_D0mass*m_D0mass,0.5));
      double D1_rec_mass = (m_db_lab - D1_con_tag_mom).mag();

  // Set D2 variables
  double D2_rec_mass;
  HepLorentzVector D2_p4;

      for (int D2_index = 0; D2_index < d0itindex.size(); D2_index++) // Loop over D0 Candidates
  {
    DTagToolIterator D2_iter_dtag = dtagTool.modes_begin() + d0itindex[D2_index];
    int D2_mode = ( * D2_iter_dtag) -> decayMode();
    if (D2_mode == 101) { // If D2 is PiPiPi0
    // Default for charged kaons and pions are just be "good tracks".
    //     // If PID is needed, require pdg=1 (=EvtRecDTag::Tight)
    if (( * D2_iter_dtag) -> type() != 1) {
      continue;
    } //Require PID be applied to Kaon/Pion daughters of the D candidate
    //cout << "PID" << endl;  
    if (( * D2_iter_dtag) -> type() != EvtRecDTag::Tight) {
      continue;
    }
    //        //cout << "EvtRecDTag" << endl;    

    D2_p4 = ( * D2_iter_dtag) -> p4();

    SmartRefVector < EvtRecTrack > D2_tracks = ( * D2_iter_dtag) -> tracks(); //List of tracks used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > D2_showers = ( * D2_iter_dtag) -> showers(); //List of showers used in the D0 reconstruction

    D1_p4.boost(-m_beta);


    // D2 Candidate Particles
    RecMdcKalTrack * D2_pip_mdcTrk = (D2_tracks[0]) -> mdcKalTrack();
    D2_pip_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    RecMdcKalTrack * D2_pim_mdcTrk = (D2_tracks[1]) -> mdcKalTrack();
    D2_pim_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    RecEmcShower *  D2_g1_shw = (D2_showers[0]) -> emcShower();
    RecEmcShower *  D2_g2_shw = (D2_showers[1]) -> emcShower();

    //Overlap Checking
    
    if (D2_pip_mdcTrk->trackId()==D1_pip_mdcTrk->trackId()
        || D2_pim_mdcTrk->trackId()==D1_pim_mdcTrk->trackId()
        || D2_g1_shw->trackId() == D2_g1_shw->trackId()
        || D2_g1_shw->trackId() == D2_g2_shw->trackId()
        || D2_g2_shw->trackId() == D2_g1_shw->trackId()
        || D2_g2_shw->trackId() == D2_g2_shw->trackId())
    {continue;}


    // Truth Matching: Compare simulated detector hits with generated trajectories to match a track with a generated particle
    // The MCT0DCHG function can return the PDG code associated with a track or the code associated its parent. It can also be used to require that the PDG codes of the particle (2nd argument) and its parent (3rd argument) match certain values. See the MCTPiPiPi0DCHG function for more details.

    int this_D0_PiPiPi0_rec=0;
    if (runNo < 0) {
    }
  if (m_D0_PiPiPi0_gen
      && MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, 421) &&  MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, 421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 0, 421) && MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 0, 421) // Are Pi+/- from D0?
      && !MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 310, 421) && !MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 310, 421)  ) // But not from K0S
     {m_D0_PiPiPi0_rec=1;
      this_D0_PiPiPi0_rec=1;
}
  if (m_D0Bar_PiPiPi0_gen
      && MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, -421) &&  MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, -421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 0, -421) && MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 0, -421) // Are Pi+/- from D0Bar?
      && !MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 310, -421) && !MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 310, -421)  ) // But not from K0S
     {m_D0Bar_PiPiPi0_rec=1;
      this_D0Bar_PiPiPi0_rec=1;
}

      HepLorentzVector D2_con_tag_mom = HepLorentzVector(D2_p4.vect(),pow(D2_p4.vect().dot(D2_p4.vect())+m_D0mass*m_D0mass,0.5));
      D2_rec_mass = (m_db_lab - D2_con_tag_mom).mag();
}//If D2 is PiPiPi0

    if (D2_mode == 102) { // If D2 is KSPi0
    // Default for charged kaons and pions are just be "good tracks".
    //     // If PID is needed, require pdg=1 (=EvtRecDTag::Tight)
    if (( * D2_iter_dtag) -> type() != 1) {
      continue;
    } //Require PID be applied to Kaon/Pion daughters of the D candidate
    //cout << "PID" << endl;  
    if (( * D2_iter_dtag) -> type() != EvtRecDTag::Tight) {
      continue;
    }
    //        //cout << "EvtRecDTag" << endl;    

    HepLorentzVector D2_unfitted_p4 = ( * D2_iter_dtag) -> p4();
    D2_unfitted_p4.boost(-m_beta);

    SmartRefVector < EvtRecTrack > D2_tracks = ( * D2_iter_dtag) -> tracks(); //List of tracks used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > D2_showers = ( * D2_iter_dtag) -> showers(); //List of showers used in the D0 reconstruction

    Vint D2_pi0Index=dtagTool.pi0Id(D2_iter_dtag, 1);
    for(EvtRecPi0Col::iterator pi0Itr = recPi0Col->begin(); pi0Itr < recPi0Col->end(); pi0Itr++){
      int pi0idx=pi0Itr-recPi0Col->begin();
      if(pi0Itr - recPi0Col->begin() != D2_pi0Index[0])continue;
    }

    // D2 Candidate Particles
    RecMdcKalTrack * D2_pip_mdcTrk = (D2_tracks[0]) -> mdcKalTrack();
    D2_pip_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    RecMdcKalTrack * D2_pim_mdcTrk = (D2_tracks[1]) -> mdcKalTrack();
    D2_pim_mdcTrk -> setPidType(RecMdcKalTrack::pion);
    EvtRecPi0Col::iterator D2_pi0= D2_pi0Index[0] + recPi0Col->begin();
    HepLorentzVector D2_pi0_p4 = D2_pi0->hiPfit() + D2_pi0->loPfit();
    RecEmcShower *  D2_hiEn_shw =  const_cast<EvtRecTrack*>((*D2_pi0)->hiEnGamma())->emcShower();
    RecEmcShower *  D2_loEn_shw =  const_cast<EvtRecTrack*>((*D2_pi0)->loEnGamma())->emcShower();
    RecEmcShower *  D2_g1_shw = (D2_showers[0]) -> emcShower();
    RecEmcShower *  D2_g2_shw = (D2_showers[1]) -> emcShower();

    if(D2_hiEn_shw->trackID()!= D2_g1_shw->trackId() && D2_hiEn_shw->trackID()!= D2_g2_shw->trackId() )
    {
      cout<<"Uh-oh! You've got the wrong pi0 there, friend."<<endl;
    }
    //Overlap Checking
    
    if (D2_pip_mdcTrk->trackId()==D1_pip_mdcTrk->trackId()
        || D2_pim_mdcTrk->trackId()==D1_pim_mdcTrk->trackId()
        || D2_g1_shw->trackId() == D2_g1_shw->trackId()
        || D2_g1_shw->trackId() == D2_g2_shw->trackId()
        || D2_g2_shw->trackId() == D2_g1_shw->trackId()
        || D2_g2_shw->trackId() == D2_g2_shw->trackId())
    {continue;}

      // Ks Vertex Fit
      std::vector<double> kschisq_vals = KsVertexFit(kspipi_tracks[0],kspipi_tracks[1],xorigin,xorigin_err);

      double decayLength = secVtxFit->decayLength();
      double decayLengthError = secVtxFit->decayLengthError();
      m_ks_flightSig=decayLength / decayLengthError;
      HepLorentzVector ks_vtxfit_p4 = secVtxFit->p4par();// Fitted Ks Momentum
      // Get VertexFitted D0-Momentum
      D2_p4=ks_vtxfit_p4+D2_pi0_p4;
      D2_p4.boost(-m_beta);

    // Truth Matching: Compare simulated detector hits with generated trajectories to match a track with a generated particle
    // The MCTKSPi0DCHG function can return the PDG code associated with a track or the code associated its parent. It can also be used to require that the PDG codes of the particle (2nd argument) and its parent (3rd argument) match certain values. See the MCTKSPi0DCHG function for more details.

    int this_D0_KSPi0_rec=0;
    if (runNo < 0) {
    }
  if (m_D0_KSPi0_gen
      && MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, 421) &&  MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, 421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 310, 421) && MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 310, 421)  ) // But not from K0S
     {m_D0_KSPi0_rec=1;
      this_D0_KSPi0_rec=1;
}
  if (m_D0Bar_KSPi0_gen
      && MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, -421) &&  MCSHPIDCHG(D2_g1_shw -> trackId(), 22, 0, -421) // Need to Test Pi0 Truth Eff.
      && MCTKPIDCHG(D2_pip_mdcTrk -> trackId(), 211, 310, -421) && MCTKPIDCHG(D2_pim_mdcTrk -> trackId(), -211, 310, -421)  ) 
     {m_D0Bar_KSPi0_rec=1;
      this_D0Bar_KSPi0_rec=1;
}
      HepLorentzVector D2_con_tag_mom = HepLorentzVector(D2_p4.vect(),pow(D2_p4.vect().dot(D2_p4.vect())+m_D0mass*m_D0mass,0.5));
      D2_rec_mass = (m_db_lab - D2_con_tag_mom).mag();
}//If D2 is KSPi0

    if (D2_mode == 105) { // If D2 is KK
    // Default for charged kaons and pions are just be "good tracks".
    //     // If PID is needed, require pdg=1 (=EvtRecDTag::Tight)
    if (( * D2_iter_dtag) -> type() != 1) {
      continue;
    } //Require PID be applied to Kaon/Pion daughters of the D candidate
    //cout << "PID" << endl;  
    if (( * D2_iter_dtag) -> type() != EvtRecDTag::Tight) {
      continue;
    }
    //        //cout << "EvtRecDTag" << endl;    

    D2_p4 = ( * D2_iter_dtag) -> p4();
    D2_p4.boost(-m_beta);

    SmartRefVector < EvtRecTrack > D2_tracks = ( * D2_iter_dtag) -> tracks(); //List of tracks used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > D2_showers = ( * D2_iter_dtag) -> showers(); //List of showers used in the D0 reconstruction


    // D2 Candidate Particles
    RecMdcKalTrack * D2_Kp_mdcTrk = (D2_tracks[0]) -> mdcKalTrack();
    D2_Kp_mdcTrk -> setPidType(RecMdcKalTrack::kaon);
    RecMdcKalTrack * D2_Km_mdcTrk = (D2_tracks[1]) -> mdcKalTrack();
    D2_Km_mdcTrk -> setPidType(RecMdcKalTrack::kaon);

    //Overlap Checking
    
    if (D2_Kp_mdcTrk->trackId()==D1_pip_mdcTrk->trackId()
        || D2_Km_mdcTrk->trackId()==D1_pim_mdcTrk->trackId()
    {continue;}


//Truth Matching
    int this_D0_KK_rec=0;
    if (runNo < 0) {
    }
  if (m_D0_KK_gen
      && MCTKPIDCHG(D2_Kp_mdcTrk -> trackId(), 321, 421, 0) && MCTKPIDCHG(D2_Km_mdcTrk -> trackId(), -321, 421, 0)  )
     {m_D0_KK_rec=1;
      this_D0_KK_rec=1;
}
  if (m_D0Bar_KK_gen
      && MCTKPIDCHG(D2_Kp_mdcTrk -> trackId(), 321, -421, 0) && MCTKPIDCHG(D2_Km_mdcTrk -> trackId(), -321, -421, 0)  ) 
     {m_D0Bar_KK_rec=1;
      this_D0Bar_KK_rec=1;
}
      HepLorentzVector D2_con_tag_mom = HepLorentzVector(D2_p4.vect(),pow(D2_p4.vect().dot(D2_p4.vect())+m_D0mass*m_D0mass,0.5));
      D2_rec_mass = (m_db_lab - D2_con_tag_mom).mag();
}//If D2 is KK


}// For D2 candidates



HepLorentzVector DTag_p4=D1_p4+D2_p4;
HepLorentzVector miss_p4=m_cm_p4-DTag_p4;
 m_Emiss=miss_p4.e();
 m_mmiss2=miss_p4.mag2();
m_ClosestToDStar_Recmass=std::min(D1_rec_mass-m_dstzero,D2_rec_mass-m_dstzero);

  } //D1 PiPiPi0 candidate loop
  /////////////////////////////////////////

m_SlowestPion_p=std::min(pion1,pion2...);
for (Pion in GoodPions)
{
  if (pion.p4().vect().mag()<m_SlowestPion_p){m_SlowestPion_p=pion.p4().vect().mag();}
}

  //Write to the mtuple: 1 entry per event
  m_tuple1 -> write();
  return sc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PiPiPi0Tags::finalize() {

  MsgStream logmess(msgSvc(), name());
  cout << "_/_/_/_/_/_/_/_/_/ all event- " << m_cout_all << endl;
  return StatusCode::SUCCESS;
}
// ************************** FUNCTIONS *****************************************************
HepLorentzVector PiPiPi0Tags::myboost(const HepLorentzVector & p1,
  const HepLorentzVector & p2) // boost particle p1 to the rest frame of p2 (boosted p1 returned as the function value)
{
  double am2(p2.m());
  double e((p1.t() * p2.t() - (p1.x() * p2.x() + p1.y() * p2.y() +
    p1.z() * p2.z())) / am2);
  double r((p1.t() + e) / (p2.t() + am2));
  return HepLorentzVector(p1.x() - r * p2.x(),
    p1.y() - r * p2.y(),
    p1.z() - r * p2.z(),
    e);
}

HepLorentzVector PiPiPi0Tags::getP4(RecEmcShower * gTrk) {

  double eraw = gTrk -> energy();
  double phi = gTrk -> phi();
  double the = gTrk -> theta();

  return HepLorentzVector(eraw * sin(the) * cos(phi),
    eraw * sin(the) * sin(phi),
    eraw * cos(the),
    eraw);
}
HepLorentzVector PiPiPi0Tags::gettrp4(RecMdcKalTrack * mdcKalTrack, int pid) {
  HepVector zhelix;
  double mass = 0;

  if (pid == 0) {
    zhelix = mdcKalTrack -> getZHelixE();
    mass = 0.000511;
  } else if (pid == 1) {
    zhelix = mdcKalTrack -> getZHelixMu();
    mass = 0.105658;
  } else if (pid == 2) {
    zhelix = mdcKalTrack -> getZHelix();
    mass = 0.139570;
  } else if (pid == 3) {
    zhelix = mdcKalTrack -> getZHelixK();
    mass = 0.493677;
  } else {
    zhelix = mdcKalTrack -> getZHelixP();
    mass = 0.938272;
  }

  double dr(0), phi0(0), kappa(0), dz(0), tanl(0);
  dr = zhelix[0];
  phi0 = zhelix[1];
  kappa = zhelix[2];
  dz = zhelix[3];
  tanl = zhelix[4];

  int charge = 0;
  if (kappa > 0.0000000001)
    charge = 1;
  else if (kappa < -0.0000000001)
    charge = -1;

  double pxy = 0;
  if (kappa != 0) pxy = 1.0 / fabs(kappa);

  double px = pxy * (-sin(phi0));
  double py = pxy * cos(phi0);
  double pz = pxy * tanl;

  double e = sqrt(pxy * pxy + pz * pz + mass * mass);

  return HepLorentzVector(px, py, pz, e);
}
HepLorentzVector PiPiPi0Tags::gettrp4fac(RecMdcKalTrack * mdcKalTrack, int pid, double fac) {
  HepVector zhelix;
  double mass = 0;

  if (pid == 0) {
    zhelix = mdcKalTrack -> getZHelixE();
    mass = 0.000511;
  } else if (pid == 1) {
    zhelix = mdcKalTrack -> getZHelixMu();
    mass = 0.105658;
  } else if (pid == 2) {
    zhelix = mdcKalTrack -> getZHelix();
    mass = 0.139570;
  } else if (pid == 3) {
    zhelix = mdcKalTrack -> getZHelixK();
    mass = 0.493677;
  } else {
    zhelix = mdcKalTrack -> getZHelixP();
    mass = 0.938272;
  }

  double dr(0), phi0(0), kappa(0), dz(0), tanl(0);
  dr = zhelix[0];
  phi0 = zhelix[1];
  kappa = zhelix[2];
  dz = zhelix[3];
  tanl = zhelix[4];

  int charge = 0;
  if (kappa > 0.0000000001)
    charge = 1;
  else if (kappa < -0.0000000001)
    charge = -1;

  double pxy = 0;
  if (kappa != 0) pxy = 1.0 / fabs(kappa);

  double px = pxy * (-sin(phi0));
  double py = pxy * cos(phi0);
  double pz = pxy * tanl;

  double e = sqrt(pxy * pxy + pz * pz + mass * mass);

  px = fac * px;
  py = fac * py;
  pz = fac * pz;
  e = sqrt(px * px + py * py + pz * pz + mass * mass);

  return HepLorentzVector(px, py, pz, e);
}
//int PiPiPi0Tags::MCGENQ(int parID, int children){
//  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), EventModel::MC::McParticleCol);
//  int numChildren=0;
//  if(mcParticleCol)
//  {
//	Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
//  SmartRefVector<Event::McParticle> gc;
//  for ( ; iter_mc != mcParticleCol->end(); iter_mc++) 
//        {
//          if((*iter_mc)->particleProperty() == parID)
//          {
//          gc=(*iter_mc)->daughterList();
//               {
//                  for(unsigned int child = 0; child < gc.size(); child++) 
//                  {
//	                  if(!(gc[child]->decayInFlight())){numChildren++}
//                    if(gc[trk1]->particleProperty()==daughterList[0])
//                    {
//                     for(unsigned int k1 = 0; k1 < gc.size(); k1++)
//		     {
//	               if(k1==pi){continue;} 
//                     if(gc[k1]->particleProperty()==321)
//		     {
//                     for(unsigned int k2 = 0; k2 < gc.size(); k2++)
//		     {
//	               if(k2==k1){continue;} 
//	               if(k2==pi){continue;} 
//                       if(gc[k2]->particleProperty()==-321)
//		       {
//			       m_Dsp_kkp=1;
//		       }//if second kaon found
//		     }//second kaon loop
//		     }//if first kaon found
//		     }//first kaon loop
//		  }//if pion found
//	       }//pion loop
//	       }//if 3 trax
//
//}

int PiPiPi0Tags::MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  int ismatched = 0;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == -1 && (GParPDG == 0 || GParPDG == thebestGParent)) {
        return thebestParent;
      }
      if (mcParPDG == -1 && GParPDG == 0) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }

  return ismatched;
}
bool PiPiPi0Tags::MCPARPID(int tkID, int mcPDG, int mcParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcPDG == 0) {
        if (thebestParent == mcParPDG) {
          ismatched = true;
        }
      } else if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}
bool PiPiPi0Tags::MCTKPID(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      int thebestGGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        if (!particles[i] -> decayFromGenerator()) {
          continue;
        }
        if (particles[i] -> decayInFlight()) {
          continue;
        }
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG < 0) {
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      } else if (mcParPDG == 0) {
        if (fabs(thebestmatchedPDG) == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG ||
            thebestGGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      } else {
        if (fabs(thebestmatchedPDG) == mcPDG && fabs(thebestParent) == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG ||
            thebestGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}

bool PiPiPi0Tags::MCSHPID(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == 0) {
        if (fabs(thebestmatchedPDG) == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          //if (temp_hits>20) {ismatched = true;}
          ismatched = true;
        }
      } else {
        if (fabs(thebestmatchedPDG) == mcPDG && fabs(thebestParent) == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}
int PiPiPi0Tags::MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  int ismatched = -1;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }

      if (mcParPDG == -1) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}
bool PiPiPi0Tags::Par2SH(int GParPDG, int mcPDG1, int mcPDG2, int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int matchedPDG = -2;
      int Parent = -2;
      int ParentMCID = -2;
      int MC_sh1 = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          matchedPDG = particles[i] -> particleProperty();
          Parent = particles[i] -> mother().particleProperty();
          ParentMCID = particles[i] -> mother().trackIndex();
          MC_sh1 = particles[i] -> trackIndex();
        }
      }
      if (matchedPDG == mcPDG1 && Parent == GParPDG) {
        for (iter_mc = mcParticles -> begin(); iter_mc != mcParticles -> end(); iter_mc++) {
          if (ParentMCID == ( * iter_mc) -> trackIndex()) {
            const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
            if (gc.size() > 0) {
              for (unsigned int ii = 0; ii < gc.size(); ii++) {
                if (gc[ii] -> particleProperty() == mcPDG2 &&
                  MC_sh1 != gc[ii] -> trackIndex()) {
                  int MCtkID2 = -2;
                  RecEmcShowerVector tracks2 = navigator -> getEmcRecShowers(gc[ii]);
                  int temp_hits2 = 0;
                  for (unsigned int z = 0; z < tracks2.size(); z++) {
                    int tkHITS = navigator -> getMcParticleRelevance(tracks2[z], gc[ii]);
                    if (tkHITS > temp_hits2) {
                      temp_hits2 = tkHITS;
                      for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
                        EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
                        if (!( * itTrk) -> isEmcShowerValid()) {
                          continue;
                        }
                        RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
                        if (emcTrk -> getShowerId().get_value() == tracks2[z] -> getShowerId().get_value()) {
                          MCtkID2 = ( * itTrk) -> trackId();
                        }
                      }
                    }
                  }
                  if (MCtkID2 == tkID2) {
                    ismatched = true;
                  }
                } // matched the 2nd daughter
              } // Loop over the 2nd daughter
            }
          }
        }
      }
    }
  }

  return ismatched;
}
bool PiPiPi0Tags::ETA3PI(int tkID1, int tkID2, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!Par2SH(111, 22, 22, shID1, shID2)) {
    return ismatched;
  }

  int eta_tk1 = -2;
  int eta_tk2 = -2;
  int eta_sh1 = -2;
  int eta_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          eta_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          eta_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          eta_sh1 = particles[i] -> mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          eta_sh2 = particles[i] -> mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (eta_tk1 != -2 &&
    eta_tk1 == eta_tk2 && eta_tk1 == eta_sh1 && eta_tk1 == eta_sh2 &&
    eta_tk2 == eta_sh1 && eta_tk2 == eta_sh2 &&
    eta_sh1 == eta_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool PiPiPi0Tags::ETAPGG(int tkID1, int tkID2, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!Par2SH(221, 22, 22, shID1, shID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_sh1 = -2;
  int etap_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh2 = particles[i] -> mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_sh1 && etap_tk1 == etap_sh2 &&
    etap_tk2 == etap_sh1 && etap_tk2 == etap_sh2 &&
    etap_sh1 == etap_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool PiPiPi0Tags::ETAP3P(int tkID1, int tkID2, int tkID3, int tkID4, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!ETA3PI(tkID3, tkID4, shID1, shID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_tk3 = -2;
  int etap_tk4 = -2;
  int etap_sh1 = -2;
  int etap_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_tk3 = -2;
  int MCID_tk4 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID3) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> trackIndex() != MCID_tk2 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk3 = particles[i] -> mother().mother().trackIndex();
          MCID_tk3 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID4) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> trackIndex() != MCID_tk2 &&
          particles[i] -> trackIndex() != MCID_tk3 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk4 = particles[i] -> mother().mother().trackIndex();
          MCID_tk4 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh2 = particles[i] -> mother().mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_tk3 && etap_tk1 == etap_tk4 && etap_tk1 == etap_sh1 && etap_tk1 == etap_sh2 &&
    etap_tk2 == etap_tk3 && etap_tk2 == etap_tk4 && etap_tk2 == etap_sh1 && etap_tk2 == etap_sh2 &&
    etap_tk3 == etap_tk4 && etap_tk3 == etap_sh1 && etap_tk3 == etap_sh2 &&
    etap_tk4 == etap_sh1 && etap_tk4 == etap_sh2 &&
    etap_sh1 == etap_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool PiPiPi0Tags::RHOPP(int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;

  int rho_tk1 = -2;
  int rho_tk2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          rho_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          rho_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (rho_tk1 != -2 && rho_tk1 == rho_tk2) {
    ismatched = true;
  }

  return ismatched;
}
bool PiPiPi0Tags::ETAPRG(int tkID1, int tkID2, int shID1) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!RHOPP(tkID1, tkID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_sh1 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_sh1 &&
    etap_tk2 == etap_sh1) {
    ismatched = true;
  }

  return ismatched;
}
bool PiPiPi0Tags::KSPP(int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;

  int ks_tk1 = -2;
  int ks_tk2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 310 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          ks_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 310 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          ks_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  if (ks_tk1 != -2 && ks_tk1 == ks_tk2) {
    ismatched = true;
  }

  return ismatched;
}

bool PiPiPi0Tags::MatchedTK(int tkID, int mcPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        if (!particles[i] -> decayFromGenerator()) {
          continue;
        }
        if (particles[i] -> decayInFlight()) {
          continue;
        }
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
        }
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}

int PiPiPi0Tags::MatchedTKID(int mcID, RecMdcKalTrackVector tracks) {
  int tkID = 0;
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  if (!navigator) {
    //cout << "EventNavigator not found" << endreq;
    tkID = -1;
    return int(tkID);
  }

  for (; iter_mc != mcParticles -> end(); iter_mc++) {
    if (mcID == ( * iter_mc) -> trackIndex()) {
      int temp_tkHITS = 0;
      for (unsigned int z = 0; z < tracks.size(); z++) {
        int tkHITS = navigator -> getMcParticleRelevance(tracks[z], * iter_mc);
        if (tkHITS > temp_tkHITS) {
          temp_tkHITS = tkHITS;
          tkID = tracks[z] -> trackId();
        }
      } // matched tracks
    }
  } // iter_mc

  return int(tkID);
}
int PiPiPi0Tags::MatchedSHID(int mcID, RecEmcShowerVector showers) {
  int shID = 0;
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  if (!navigator) {
    //cout << "EventNavigator not found" << endreq;
    shID = -1;
    return int(shID);
  }
  int shIDtemp = -1;
  for (; iter_mc != mcParticles -> end(); iter_mc++) {
    if (mcID == ( * iter_mc) -> trackIndex()) {
      int temp_shHITS = 0;
      for (unsigned int z = 0; z < showers.size(); z++) {
        int shHITS = navigator -> getMcParticleRelevance(showers[z], * iter_mc);
        if (shHITS > temp_shHITS) {
          temp_shHITS = shHITS;
          shIDtemp = showers[z] -> getShowerId().get_value();
        }
      } // matched showers
    }
  } // iter_mc

  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (emcTrk -> getShowerId().get_value() == shIDtemp) {
      shID = ( * itTrk) -> trackId();
    }
  }

  return int(shID);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
std::string PiPiPi0Tags::ParticleName(int pdg) {
  IPartPropSvc * p_PartPropSvc;
  static
  const bool CREATEIFNOTTHERE(true);
  StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
  if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
    std::cout << " Could not initialize Particle Properties Service" << std::endl;
    return "0";
  }
  m_particleTable = p_PartPropSvc -> PDT();

  std::string name;
  if (m_particleTable -> particle(pdg))
    name = m_particleTable -> particle(pdg) -> name();
  else if (m_particleTable -> particle(-pdg))
    name = "anti " + m_particleTable -> particle(-pdg) -> name();

  return name;
}
void PiPiPi0Tags::DumpTree() {

  SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(),
    EventModel::MC::McParticleCol);
  if (mcParticleCol) {
    cout << "-------------------------" << endl;
    //cout << "Run- " << m_run << ", Event- " << m_event << endl;

    //////////////////////////////
    /// Dump vertices
    cout << "Vertices- " << endl;

    // "Cluster" is first particle before BesEvtGen particle
    bool foundClusterAsMother = false;
    //if(!m_BesEvtGenOnly) foundClusterAsMother = true;
    Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
    for (; iter_mc != mcParticleCol -> end(); iter_mc++) {

      if (!( * iter_mc) -> primaryParticle()) {
        if (( * iter_mc) -> mother().particleProperty() == 91) foundClusterAsMother = true;
      }
      // took this out for qqbar MC
      //if(!foundClusterAsMother) continue;

      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();

      if (gc.size() > 0) {
        cout << " " << ParticleName(( * iter_mc) -> particleProperty()) << " [" <<
          ( * iter_mc) -> trackIndex() << "] -> ";

        for (unsigned int ii = 0; ii < gc.size(); ii++) {
          if (ii != (gc.size() - 1))
            cout << ParticleName(gc[ii] -> particleProperty()) << " [" <<
            gc[ii] -> trackIndex() << "], ";
          else
            cout << ParticleName(gc[ii] -> particleProperty()) <<
            " [" << gc[ii] -> trackIndex() << "]" << endl;
        }

      } // End of "gc.size() > 0" IF

    } // End of "McParticleCol" FOR LOOP

    //////////////////////////////////////
    /// Dump particles

    cout << endl << "Particles-  [#Children, primParticle, leafParticle, decayFromGen, decayInFlight] " <<
      endl;

    foundClusterAsMother = false;
    //if(!m_BesEvtGenOnly) foundClusterAsMother = true;

    bool firstDecayInFlight = true;
    iter_mc = mcParticleCol -> begin();
    for (; iter_mc != mcParticleCol -> end(); iter_mc++) {

      if (!( * iter_mc) -> primaryParticle()) {
        if (( * iter_mc) -> mother().particleProperty() == 91) foundClusterAsMother = true;
      }
      // took this out for qqbar MC
      //if(!foundClusterAsMother) continue;

      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
      int numChildren = gc.size();

      string primaryParticle = "F";
      if (( * iter_mc) -> primaryParticle()) primaryParticle = "T";

      string leafParticle = "F";
      if (( * iter_mc) -> leafParticle()) leafParticle = "T";

      string decayFromGen = "F";
      if (( * iter_mc) -> decayFromGenerator()) decayFromGen = "T";

      string decayInFlight = "F";
      if (( * iter_mc) -> decayInFlight()) {
        decayInFlight = "T";

        if (firstDecayInFlight) {
          cout << endl;
          firstDecayInFlight = false;
        }
      }

      //if(!(*iter_mc)->decayFromGenerator()) 
      //{
      cout << " " << ( * iter_mc) -> trackIndex() << "- " <<
        ParticleName(( * iter_mc) -> particleProperty()) <<
        "  ID = " << ( * iter_mc) -> particleProperty() <<
        " p4 = " << ( * iter_mc) -> initialFourMomentum() << " [" <<
        numChildren << ", " <<
        primaryParticle << ", " << leafParticle << ", " <<
        decayFromGen << ", " << decayInFlight << "]" <<
        endl;
      //}
    } // End of "McParticleCol" FOR LOOP

  }
  cout << endl << endl;
}
int PiPiPi0Tags::countPar(int parentID, int targetID) {
  SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(),
    EventModel::MC::McParticleCol);
  int numcount = 0;
  Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
  for (; iter_mc != mcParticleCol -> end(); iter_mc++) {
    if (( * iter_mc) -> particleProperty() == parentID &&
      ( * iter_mc) -> decayFromGenerator() && !( * iter_mc) -> decayInFlight()) {
      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
      if (gc.size() > 0) {
        for (unsigned int ii = 0; ii < gc.size(); ii++) {
          if (gc[ii] -> particleProperty() == targetID &&
            gc[ii] -> decayFromGenerator() && !gc[ii] -> decayInFlight()) {
            numcount++;
          }
          if (gc[ii] -> daughterList().size() > 0) {
            const SmartRefVector < Event::McParticle > & ggc = gc[ii] -> daughterList();
            for (unsigned int iii = 0; iii < ggc.size(); iii++) {
              if (ggc[iii] -> particleProperty() == targetID &&
                ggc[iii] -> decayFromGenerator() && !ggc[iii] -> decayInFlight()) {
                numcount++;
              }
              if (ggc[iii] -> daughterList().size() > 0) {
                const SmartRefVector < Event::McParticle > & gggc = ggc[iii] -> daughterList();
                for (unsigned int iiii = 0; iiii < gggc.size(); iiii++) {
                  if (gggc[iiii] -> particleProperty() == targetID &&
                    gggc[iiii] -> decayFromGenerator() && !gggc[iiii] -> decayInFlight()) {
                    numcount++;
                  }
                  if (gggc[iiii] -> daughterList().size() > 0) {
                    const SmartRefVector < Event::McParticle > & ggggc = gggc[iiii] -> daughterList();
                    for (unsigned int iiiii = 0; iiiii < ggggc.size(); iiiii++) {
                      if (ggggc[iiiii] -> particleProperty() == targetID &&
                        ggggc[iiiii] -> decayFromGenerator() && !ggggc[iiiii] -> decayInFlight()) {
                        numcount++;
                      }
                      if (ggggc[iiiii] -> daughterList().size() > 0) {
                        const SmartRefVector < Event::McParticle > & gggggc = ggggc[iiiii] -> daughterList();
                        for (unsigned int iiiiii = 0; iiiiii < gggggc.size(); iiiiii++) {
                          if (gggggc[iiiiii] -> particleProperty() == targetID &&
                            gggggc[iiiiii] -> decayFromGenerator() && !gggggc[iiiiii] -> decayInFlight()) {
                            numcount++;
                          }
                          if (gggggc[iiiiii] -> daughterList().size() > 0) {
                            const SmartRefVector < Event::McParticle > & ggggggc = gggggc[iiiiii] -> daughterList();
                            for (unsigned int iiiiiii = 0; iiiiiii < ggggggc.size(); iiiiiii++) {
                              if (ggggggc[iiiiiii] -> particleProperty() == targetID &&
                                ggggggc[iiiiiii] -> decayFromGenerator() && !ggggggc[iiiiiii] -> decayInFlight()) {
                                numcount++;
                              }
                            } // iiiiiii
                          }
                        } // iiiiii
                      }
                    } // iiiii
                  }
                } // iiii
              }
            } // iii
          }
        } // ii
      }
    } // parent found
  }
  return numcount;
}
bool PiPiPi0Tags::isTagTrue(DTagToolIterator iter_dtag) {
    bool isfound = false;
    SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
    if (!navigator) {
      return isfound;
    }
    int m_mode = ( * iter_dtag) -> decayMode();
    int m_charge = ( * iter_dtag) -> charge();
    SmartRefVector < EvtRecTrack > tracks = ( * iter_dtag) -> tracks();
    SmartRefVector < EvtRecTrack > showers = ( * iter_dtag) -> showers();
    vector < int > numchan;
    numchan.clear();
    //K+/- :1 
    //Pi+/-:2
    //Pi0  :3
    //Eta  :4
    //Ks   :5
    //eta'(pipieta)         :6
    //eta'(rhogamma)        :7
    //eta'(pipieta(PiPiPi0)):8
    //eta(PiPiPi0)          :9

    //cout << " *********** the mode = " << m_mode << endl;

    //if (m_mode != 404) {return isfound;}

    if (m_mode == 400) {
      numchan.push_back(5);
      numchan.push_back(1);
    } // DstoKsK
    else if (m_mode == 401) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
    } // DstoKDStar
    else if (m_mode == 402) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(3);
    } // DstoKsDStar0
    else if (m_mode == 403) {
      numchan.push_back(5);
      numchan.push_back(5);
      numchan.push_back(2);
    } // DstoKsKsPi
    else if (m_mode == 404) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(3);;
    } // DstoKDStarPi0
    else if (m_mode == 405) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKsKplusPiPi
    else if (m_mode == 406) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKsKminusPiPi
    else if (m_mode == 407) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKDStarPiPi
    else if (m_mode == 420) {
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPi0
    else if (m_mode == 421) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoPiPiPi
    else if (m_mode == 422) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPiPiPi0
    else if (m_mode == 423) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoPiPiPiPiPi
    else if (m_mode == 424) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPiPiPiPiPi0
    else if (m_mode == 425) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(3);
    } // DstoPiPiPiPi0Pi0 <-- this is wrong
    else if (m_mode == 440) {
      numchan.push_back(2);
      numchan.push_back(4);
    } // DstoPiEta
    else if (m_mode == 441) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(4);
    } // DstoPiPi0Eta
    else if (m_mode == 442) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(4);
    } // DstoPiPiPiEta
    else if (m_mode == 450) {
      numchan.push_back(2);
      numchan.push_back(9);
    } // DstoPiEtaPiPiPi0
    else if (m_mode == 451) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(9);
    } // DstoPiPi0EtaPiPiPi0
    else if (m_mode == 452) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(9);
    } // DstoPiPiPiEtaPiPiPi0
    else if (m_mode == 460) {
      numchan.push_back(2);
      numchan.push_back(6);
    } // DstoPiEPPiPiEta
    else if (m_mode == 461) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(6);
    } // DstoPiPi0EPPiPiEta
    else if (m_mode == 470) {
      numchan.push_back(2);
      numchan.push_back(8);
    } // DstoPiEPPiPiEtaPiPiPi0
    else if (m_mode == 471) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(8);
    } // DstoPiPi0EPPiPiEtaPiPiPi0
    else if (m_mode == 480) {
      numchan.push_back(2);
      numchan.push_back(7);
    } // DstoPiEPRhoGam
    else if (m_mode == 481) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(7);
    } // DstoPiPi0EPRhoGam
    else if (m_mode == 500) {
      numchan.push_back(5);
      numchan.push_back(2);
    } // DstoKsPi
    else if (m_mode == 501) {
      numchan.push_back(5);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoKsPiPi0
    else if (m_mode == 502) {
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoDStarPi
    else if (m_mode == 503) {
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoDStarPiPi0
    else if (m_mode == 504) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(1);
    } // DstoKKK
    bool isparmatched = true;
    int npi = 0;
    int nka = 0;
    int np0 = 0;
    int net = 0;
    int nks = 0;
    int nep6 = 0;
    int nep7 = 0;
    int nep8 = 0;
    int net9 = 0;
    for (int i = 0; i < ( * iter_dtag) -> numOfChildren(); i++) {
      //cout << " EACH PARTICLES = " << numchan[i] << endl;
      if (numchan[i] == 1) { // Kaon
        nka++;
        if (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 321, 0, 431)) &&
          !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 321, 0, -431))) {
          isparmatched = false;
        }
        //cout << " K1 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),321,0, 431) << endl
        //   << " K2 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),321,0,-431) << endl
        //   << " K matched = " << isparmatched << endl;
      } else if (numchan[i] == 2) { // pion
        npi++;
        if (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 0, 431)) &&
          !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 0, -431))) {
          isparmatched = false;
        }
        //cout << " P1 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,0, 431) << endl
        //   << " P2 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,0,-431) << endl
        //   << " P matched = " << isparmatched << endl;
      } else if (numchan[i] == 3) { // pi0->gg
        np0++;
        //if ( (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ||
        //   (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) ) {isparmatched=false;}

        //cout << " P01 = " << (m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) << endl
        //   << " P02 = " << (m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431)) << endl
        //   << " P03 = " << (m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) << endl
        //   << " P04 = " << (m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431)) << endl
        //   << " P0 matched = " << isparmatched << endl;
        //cout << " Try1 = " << (m_charge>0 &&
        //		       (MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431) ||
        //			MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431))) << endl
        //   << " Try2 = " << (m_charge<0 &&
        //		       (MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431) ||
        //			MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) << endl;

        if (!(m_charge > 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 111, 431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 111, 431))) &&
          !(m_charge < 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 111, -431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 111, -431)))) {
          isparmatched = false;
        }

      } else if (numchan[i] == 4) { // eta->gg
        net++;
        //if ( (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,221, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,221,-431))) ||
        //   (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,221, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,221,-431))) ) {isparmatched=false;} 

        if (!(m_charge > 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, 431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, 431))) &&
          !(m_charge < 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, -431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 5) { // Ks
        nks++;
        if ((!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 310, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 310, -431))) ||
          (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 310, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 310, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 6) { // eta' -> pipi eta(->gg)
        nep6++;
        if ((!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 331, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 331, -431))) ||
          (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 331, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 331, -431))) ||
          (!(m_charge > 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, 431)) &&
            !(m_charge < 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, -431))) ||
          (!(m_charge > 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, 431)) &&
            !(m_charge < 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 7) { // eta' -> rho gamma
                nep7++;
        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,113, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,113,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,113, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,113,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,331, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,331,-431))) ) {isparmatched=false;} }
      else if(numchan[i]==8){ // eta' -> pipi eta(->PiPiPi0)
        nep8++;

        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-4]->trackId(),211,331, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-4]->trackId(),211,331,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-3]->trackId(),211,331, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-3]->trackId(),211,331,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ) {isparmatched=false;} }
      else if(numchan[i]==9){ // eta -> PiPiPi0
        net9++;
        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) ) {isparmatched=false;} }
      //if (m_mode == 404)
      // {
      //   cout << " EACH PARTICLES = " << numchan[i] << endl;
      // }
    }
      if (!isparmatched) {return isfound;}

  //cout << "============ B "<< isparmatched << endl;

  int totalPi = npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9;
  int totalKa = nka;
  int totalP0 = np0+nep8+net9;
  int totalET = net+nep6+nep8+net9;
  int totalKs = nks;
  int totalRo = nep7;
  int totalEP = nep6+nep7+nep8;
  if ( m_charge>0 &&
       (countPar( 431, 211)+countPar( 431,-211))==totalPi &&
       (countPar( 431, 321)+countPar( 431,-321))==totalKa &&
       countPar( 431,111)==totalP0&&countPar( 431, 221)==totalET&&
       countPar( 431,310)==totalKs&&countPar( 431, 113)==totalRo&&countPar( 431,331)==totalEP ) {isfound=true;}
  if ( m_charge<0 &&
       (countPar(-431, 211)+countPar(-431,-211))==totalPi &&
       (countPar(-431, 321)+countPar(-431,-321))==totalKa &&
       countPar(-431,111)==totalP0&&countPar(-431, 221)==totalET&&
       countPar(-431,310)==totalKs&&countPar(-431, 113)==totalRo&&countPar(-431,331)==totalEP ) {isfound=true;}

  //cout << "np0 = " << np0 << " nep8 = " << nep8 << " net0 = " << net9 << endl;
  //cout << "============ C "<< isfound << endl;

  //if(m_mode==450 && m_charge>0)
  //if(m_mode==401)
  //if (m_charge>0)
  //{
  //  cout << "------------------------" << endl;
  //  cout << " mode = " << m_mode << endl;
  //  cout << " pi+-= " << countPar( 431, 211)+countPar( 431,-211) << " expected = " << totalPi << endl
  //       << " K+- = " << countPar( 431, 321)+countPar( 431,-321) << " expected = " << totalKa << endl
  //       << " pi0 = " << countPar( 431, 111) << " expected = " << totalP0 << endl
  //       << " eta = " << countPar( 431, 221) << " expected = " << totalET << endl
  //       << " Ks  = " << countPar( 431, 310) << " expected = " << totalKs << endl
  //       << " Ro  = " << countPar( 431, 113) << " expected = " << totalRo << endl
  //       << " EP  = " << countPar( 431, 331) << " expected = " << totalEP << endl;
  //}

  return isfound;
}

std::vector<double> PiPiPi0Tags::KsVertexFit(EvtRecTrack* trk1,EvtRecTrack* trk2,Hep3Vector beamIP,HepSymMatrix beamIPerror)
{
  std::vector<double> chisq_vals;
  double kschisq = -9.0;
  VertexParameter IPVtxPar;

  IPVtxPar.setVx(beamIP);
  IPVtxPar.setEvx(beamIPerror);

  vtxFit = VertexFit::instance();
  secVtxFit = SecondVertexFit::instance();

  HepPoint3D PointOrigin(0., 0., 0.);
  HepSymMatrix PointOrigin_err(3, 0);
  PointOrigin_err[0][0] = 1E+6;
  PointOrigin_err[1][1] = 1E+6;
  PointOrigin_err[2][2] = 1E+6;

  VertexParameter vtxPar;
  vtxPar.setVx(PointOrigin);
  vtxPar.setEvx(PointOrigin_err);

  RecMdcKalTrack *trkPi1 = trk1->mdcKalTrack();
  trkPi1->setPidType(RecMdcKalTrack::pion);
  WTrackParameter wtpPi1 = WTrackParameter(m_pionmass, trkPi1->getZHelix(), trkPi1->getZError());

  RecMdcKalTrack *trkPi2 = trk2->mdcKalTrack();
  trkPi2->setPidType(RecMdcKalTrack::pion);
  WTrackParameter wtpPi2 = WTrackParameter(m_pionmass, trkPi2->getZHelix(), trkPi2->getZError());

  //First fit to obtain decay vertex
    vtxFit->init();
      vtxFit->AddTrack(0, wtpPi1);
        vtxFit->AddTrack(1, wtpPi2);
          vtxFit->AddVertex(0, vtxPar, 0, 1);
  
            vtxFit->Fit(0);
              if (vtxFit->Fit(0))
                {
                    vtxFit->Swim(0);
       kschisq = vtxFit->chisq(0);
    chisq_vals.push_back(kschisq);
  }
  // build the virtual track
     vtxFit->BuildVirtualParticle(0);
  
       //Obtain origin vertex
         secVtxFit->init();
           secVtxFit->AddTrack(0, vtxFit->wVirtualTrack(0));
             secVtxFit->setVpar(vtxFit->vpar(0));
               secVtxFit->setPrimaryVertex(IPVtxPar);
                 secVtxFit->Fit();
                   wKstrk = secVtxFit->wpar();
  
                     chisq_vals.push_back(secVtxFit->chisq());
  
                       wTrk_Ks.clear();
                         for (int i = 0; i < 2; i++)
                           {
                               WTrackParameter itrk = vtxFit->wtrk(i);
                                   wTrk_Ks.push_back(itrk);
                                     }
  
                                       return chisq_vals;
  
                                       }                          
