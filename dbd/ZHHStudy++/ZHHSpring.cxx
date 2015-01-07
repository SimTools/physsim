//*****************************************************************************
//* =====================
//*  ZHHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZHH generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "ZHHSpring.h"

#include "TSystem.h"
#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZHHSpring)
ClassImp(ZHHSpringBuf)
ClassImp(ZHHBases)

#ifdef WITH_DBD_STANDARD
extern "C" {
  //  extern void userin_();
  extern void whizard_spectrum_(int& isrbm, float& roots, double& frbm1, double& frbm2, 
				float& ebm, double& dpdebm, float& embm, float& epbm);
  extern void isr_whizard_spectrum_(double& frisrm, double& frisrp, double& frisr1,
				    double& frisr2, double& frisr3, double& frist4,
				    float& alpha0, int& isr_llr_order, float& roots,
				    float& embm, float& epbm, float qv[4],
				    double pisrdbl[2][4], double& wat, float& rs0);
}
#endif

//-----------------------------------------------------------------------------
// ==============================
//  class ZHHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHHSpring::ZHHSpring(const char      *name,
                     const char      *title,
                            ZHHBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new ZHHSpringBuf("ZHHSpringBuf",
                               "ZHHSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new ZHHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZHHSpring::~ZHHSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZHHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZHHBases *bs = static_cast<ZHHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZHHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZHHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZHHBases     *bases   = (ZHHBases*)((ZHHSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
#ifdef WITH_DBD_STANDARD
  fNparton = 7;
#else
  fNparton = 5;
#endif

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // h1
  pv[1] = bases->fP[1]; // h2
  pv[3] = bases->fP[2]; // f1 from Z
  pv[4] = bases->fP[3]; // f2 from Z
  pv[2] = pv[3] + pv[4];
#ifdef WITH_DBD_STANDARD
  pv[5] = bases->fP[4];
  pv[6] = bases->fP[5];
#endif

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2ZHH        = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2Z          = bases->GetQ2Z();
  fCosThetaF    = bases->GetCosThetaF();
  fPhiF         = bases->GetPhiF();
#ifdef WITH_DBD_STANDARD
  ANL4DVector *fQZHH = &(bases->fQV);
  TVector3 boostv = fQZHH->BoostVector();
#else
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);
#endif

#ifdef WITH_DBD_STANDARD
  for (Int_t i=0; i<fNparton-2; i++) pv[i].Boost(boostv);
#else
  for (Int_t i=0; i<fNparton; i++) pv[i].Boost(boostv);
#endif

  // ----------------------------------------------
  //  Set final state parton infomation
  // ----------------------------------------------

  static TVector **qp = 0;
  if (!qp) {
    qp = new TVector* [fNparton];
    for (Int_t i=0; i<fNparton; i++) qp[i] = new TVector(4);
  }
  for (Int_t i=0; i<fNparton; i++) {
    TVector &q = *qp[i];
    q(0) = pv[i].E ();
    q(1) = pv[i].Px();
    q(2) = pv[i].Py();
    q(3) = pv[i].Pz();
  }
  Int_t    idh     = 25;                          // PDG code for H
  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    idf     = bases->f3Ptr->GetPID   ();   // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge();   // F charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // F mass
  Double_t m4      = m3;                          // F mass
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	  // shower level
  Int_t    icf     = 2;                           // color flux id
  Double_t rq2z    = pv[2].Mag();
  Int_t    hel1    = bases->fHelFinal[2];
  Int_t    hel2    = bases->fHelFinal[3];

  Double_t mass    = bases->GetMass();

#ifdef WITH_DBD_STANDARD
  // ISR
  Int_t   idgam    = 22;
#endif

  //                                No. PID  Mass  Charge   pv  Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, idh, mass,    0., *qp[0], 0, 0,  0,    0,   0,     0, 0, 204);
  new (partons[1]) JSFSpringParton(2, idh, mass,    0., *qp[1], 0, 0,  0,    0,   0,     0, 0, 205);
  new (partons[2]) JSFSpringParton(3, idz, rq2z,    0., *qp[2], 2, 4,  0,    0,   0,     0, 0, -20);
  new (partons[3]) JSFSpringParton(4, idf,   m3,  chrg, *qp[3], 0, 0,  3, hel1, icf, islev, 0, 202);
  new (partons[4]) JSFSpringParton(5,-idf,   m4, -chrg, *qp[4], 0, 0,  3, hel2, icf, islev, 0, 203);
#ifdef WITH_DBD_STANDARD								       
  new (partons[5]) JSFSpringParton(6, idgam, 0.,    0., *qp[5], 0, 0,  0,    1,    0,    0, 0, 201);
  new (partons[6]) JSFSpringParton(7, idgam, 0.,    0., *qp[6], 0, 0,  0,    1,    0,    0, 0, 200);
#endif


  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHHBases::ZHHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 120.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fISRBM      ( 1),
	   fISR_LLA_Order (0),
	   fStoreRemnants (1),
	   fStoreBeams (0),
           fPole      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2ZHH     (0.),
           fQ2Z       (0.),
           fQ2HH      (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fXQ2Z      (0.),
           fCosThetaF (0.),
           fPhiF      (0.),
           fXQ2HH     (0.),
           fCosThetaH (0.),
           fPhiH      (0.),
           fR_BW_m    (0),
           fR_BW_p    (0),
           fR_BS_m    (0),
           fR_BS_p    (0),
           fR_ISR_var (0),
           fR_ISR_side(0),
	   fR_BM_m    (0),
	   fR_BM_p    (0),
	   fR_ISR_m   (0),
	   fR_ISR_p   (0),
	   fR_ISR_1   (0),
	   fR_ISR_2   (0),
	   fR_ISR_3   (0),
	   fR_ISR_4   (0)
{
  //  Constructor of bases.  Default parameter should be initialized here
  //
  // --------------------------------------------
  //  Get parameters from jsf.conf, if specified
  // --------------------------------------------

  cout << "Init zhbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("ZHHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.CosthHRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.PhiHOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.MassH","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Ecm","500."));        // E_cm (1TeV)
  ins >> fEcmInit;

#ifdef WITH_DBD_STANDARD
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.ISRBM","3"));
  ins >> fISRBM;
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.ISR_LLA_Order","0"));
  ins >> fISR_LLA_Order;
// if fISR_LLA_Order is gt 0, LLA form with order, LLA_Order, is used 
// to generate ISR spectrum
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Store_Remnants","1"));
  ins >> fStoreRemnants;
// (0, 1)=(not save, save) ISR photon info.
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Store_Beams","0"));
  ins >> fStoreBeams;
// (0,1)=(not save, save) initial e+/e- beam after beamstrahlung

  if( fISRBM < 100 &&  ( fStoreRemnants != 0 || fISR_LLA_Order != 0 )) {
    std::cout << "Fatal error in ZHHBases .. " << std::endl;
    std::cout << "  Input parameter .. fISRBM=" << fISRBM << " fStoreRemnants=" << fStoreRemnants 
              << " fISR_LLA_Order=" << fISR_LLA_Order << std::endl;
    std::cout << "  But fISRBM should be > 100 if fStoreRemnants or fISR_LLA_Order .ne. 0 " << std::endl;
    exit(-1);
  }      
  fLumiFileDirectory=gJSF->Env()->GetValue("ZHHBases.LumiFileDirectory","");
#else
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.BeamWidth","0.002")); // Beam spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;
#endif

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHHBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

#ifdef WITH_DBD_STANDARD
  if( fISRBM > 100 ) {
    if( fLumiFileDirectory.size() < 1 ) {
      std::cout << "Error!! ZHHBases.LumiFile_Directory is not set, though ISRBM> 100" << std::endl;
      std::cout << "  ISRBM=" << fISRBM << std::endl;
      exit(-1);
    }
    gSystem->Setenv("LUMI_LINKER",(fLumiFileDirectory+std::string("/lumi_linker_000")).c_str());
    gSystem->Setenv("PHOTONS_B1",(fLumiFileDirectory+std::string("/photons_beam1_linker_000")).c_str());
    gSystem->Setenv("PHOTONS_B2",(fLumiFileDirectory+std::string("/photons_beam2_linker_000")).c_str());
    gSystem->Setenv("EBEAM",(fLumiFileDirectory+std::string("/ebeam_in_linker_000")).c_str());
    gSystem->Setenv("PBEAM",(fLumiFileDirectory+std::string("/pbeam_in_linker_000")).c_str());
  }
#endif

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fXQ2HH         , 0., 1., 1, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2Z          , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaF, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiF     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaH, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiH     , fXL[5], fXU[5], 1, 1);

#ifndef WITH_DBD_STANDARD
  //--
  //  Beamstrahlung and natural Ebeam spread
  //--
  if (fBeamStr == 1) {
    DefineVariable(fR_BW_m, 0., 1., 0, 1);
    DefineVariable(fR_BW_p, 0., 1., 0, 1);
    DefineVariable(fR_BS_m, 0., 1., 0, 1);
    DefineVariable(fR_BS_p, 0., 1., 0, 1);
  }
 
  //--
  //  ISR
  //--
  if (fISR==1) {
    DefineVariable(fR_ISR_var,  0., 1., 0, 1);
    DefineVariable(fR_ISR_side, 0., 1., 0, 1);
  }
#else
    DefineVariable(fR_BM_m, 0., 1., 0, 1);
    DefineVariable(fR_BM_p, 0., 1., 0, 1);
    DefineVariable(fR_ISR_m, 0., 1., 0, 1);
    DefineVariable(fR_ISR_p, 0., 1., 0, 1);
    DefineVariable(fR_ISR_1, 0., 1., 0, 1);
    DefineVariable(fR_ISR_2, 0., 1., 0, 1);
    DefineVariable(fR_ISR_3, 0., 1., 0, 1);
    DefineVariable(fR_ISR_4, 0., 1., 0, 1);
#endif

  // --------------------------------------------
  //  Set Bases integration parameters
  // --------------------------------------------
  SetNoOfSample(80000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 10);
  SetIteration2(0.05, 20);

}
// --------------------------
//  D-tor
// --------------------------
ZHHBases::~ZHHBases()
{
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZHHBases::Func()
{
  //  Bases Integrand
  //

  Double_t bsWeight = 1.; // Jacobian factor
  Double_t eminus;        // E_e- after ISR and beamstrahlung
  Double_t eplus;         // E_e+ after ISR and beamstrahlung

  // --------------------------------------------
  //  Beamstrahlung
  // --------------------------------------------
#ifndef WITH_DBD_STANDARD
  if (fBeamStr == 1) {
    bsWeight = fBM->GetWeight(fR_BW_m, fR_BW_p,
                              fR_BS_m, fR_BS_p,
                              eminus, eplus);
#if 0
    eplus   *= fBM->GetNominalEnergy();
    eminus  *= fBM->GetNominalEnergy();
#else
    eplus   *= fEcmInit/2.;
    eminus  *= fEcmInit/2.;
#endif
    fEcmIP   = 2.*TMath::Sqrt(eplus*eminus);     // ignore electron mass 
  } else {
    fEcmIP   = fEcmInit;      
    eplus    = fEcmInit/2.;
    eminus   = eplus;
  }
#endif
   
  // --------------------------------------------
  //  Initial State Radiation
  // --------------------------------------------
#ifndef WITH_DBD_STANDARD
  static const Double_t kFisr = 2. * kAlpha0 / TMath::Pi();

  Double_t beta_isr = kFisr * (2.*(TMath::Log(fEcmIP/kM_e) - 1.));

  if (fISR == 1) {
    Double_t xisr    = TMath::Power(fR_ISR_var,1./beta_isr); // Ephoton / Ebeam 
             fEcmIP *= TMath::Sqrt(1. - xisr);                // reduced Ecm
    
    if (fR_ISR_side > 0.5) {
      eplus  *= 1. - xisr;
    } else {
      eminus *= 1. - xisr;
    }
  } // if Bremsstrahlung


  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung
#endif

#ifdef WITH_DBD_STANDARD
  fEcmIP = fEcmInit;
  // beamstrahlung
  //  bool iBSDebug = true;
  bool iBSDebug = false;
  if (iBSDebug) cerr << "------Print the information for BS and ISR------" << endl;
  double wat = 1.;
  int isrbm = fISRBM;
  float roots = fEcmIP;
  double frbm1 = fR_BM_m;
  double frbm2 = fR_BM_p;
  float ebm = roots/2;
  double dpdebm = 1.;
  float embm = roots/2;
  float epbm = roots/2;
  float rs0 = roots;
  if (iBSDebug)  cerr << "Ebeam0: " << embm << " " << epbm << endl;
  whizard_spectrum_(isrbm,roots,frbm1,frbm2,ebm,dpdebm,embm,epbm);
  wat *= dpdebm;
  if (iBSDebug)  cerr << "Ebeam1: " << embm << " " << epbm << " DPDE: " << dpdebm << " RS: " << roots << " Weight: " << wat << endl;
  rs0 = 2*TMath::Sqrt(embm*epbm);
  if (rs0 < 5.) return 0.;

  // ISR
  double frisrm = fR_ISR_m;
  double frisrp = fR_ISR_p;
  double frisr1 = fR_ISR_1;
  double frisr2 = fR_ISR_2;
  double frisr3 = fR_ISR_3;
  double frisr4 = fR_ISR_4;
  float alpha0 = kAlpha0;
  int isr_lla_order = fISR_LLA_Order;
  if (iBSDebug) cerr << "Alpha: " << alpha0 << " LLA: " << isr_lla_order << endl;
  float qv[4];
  double pisrdbl[2][4];
  isr_whizard_spectrum_(frisrm,frisrp,frisr1,frisr2,frisr3,frisr4,
			alpha0,isr_lla_order,roots,embm,epbm,
			qv,pisrdbl,wat,rs0);
  if (iBSDebug) {
    cerr << "E_ISR: " << rs0 << " Weight: " << wat << endl;
    cerr << "R_ISR: " << qv[0] << " " << qv[3] << endl;
    cerr << "P_ISR: " << pisrdbl[0][0] << " " << pisrdbl[1][0] << endl;
  }
  bsWeight *= wat;
  fEcmIP = rs0;
  fP[4] = ANL4DVector(pisrdbl[0][0],pisrdbl[0][1],pisrdbl[0][2],pisrdbl[0][3]);
  fP[5] = ANL4DVector(pisrdbl[1][0],pisrdbl[1][1],pisrdbl[1][2],pisrdbl[1][3]);
  fQV   = ANL4DVector(qv[0],qv[1],qv[2],qv[3]);
  if (iBSDebug) {
    cerr << "QZHH: " << fQV.E() << " " << fQV.Px() << " " << fQV.Py() << " " << fQV.Pz() << endl;
    cerr << "PISR1: " << fP[4].E() << " " << fP[4].Px() << " " << fP[4].Py() << " " << fP[4].Pz() << endl;
    cerr << "PISR2: " << fP[5].E() << " " << fP[5].Px() << " " << fP[5].Py() << " " << fP[5].Pz() << endl;
  }
  if (rs0 < 0.) return 0.;
#endif

  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  Double_t weight = 1;

  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < 2*fMass + m3 + m4) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average for e- is taken
  //  care of here
  SelectHelicities(weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  fQ2ZHH = fEcmIP*fEcmIP;

  // Z
  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - 2*fMass;
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // HH system
  Double_t qhh2min = TMath::Power(2*fMass,2);
  Double_t qhh2max = TMath::Power(rs - TMath::Sqrt(fQ2Z),2);
           fQ2HH   = qhh2min + (qhh2max - qhh2min)*fXQ2HH;
  bsWeight *= qhh2max - qhh2min;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch zbranch (fQ2Z  , fCosThetaF, fPhiF, m3*m3, m4*m4);
  GENBranch hhbranch(fQ2HH , fCosThetaH, fPhiH, fMass*fMass, fMass*fMass);
  GENBranch cmbranch(fQ2ZHH, fCosTheta , fPhi , &hhbranch, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill( 5, fCosThetaF        , (bsWeight*sigma));
  Xh_fill( 6, fPhiF             , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2HH), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaH        , (bsWeight*sigma));
  Xh_fill( 9, fPhiH             , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(12, (Double_t)fZMode  , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZHHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2hh  = cmbranch.GetM12();
  Double_t q2z   = cmbranch.GetM22();
  Double_t coshh = cmbranch.GetCosTheta();
  Double_t phihh = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2hh, q2z, cmframe, coshh, phihh, 0);
  ANL4DVector phh = phaseCM.GetFourMomentum(0);
  ANL4DVector pz  = phaseCM.GetFourMomentum(1);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  GENBranch &hhbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosh = hhbranch.GetCosTheta();
  Double_t phih = hhbranch.GetPhi     ();
  Double_t mh12 = hhbranch.GetM12();
  Double_t mh22 = hhbranch.GetM22();
  GENFrame hhframe = phaseCM.GetFrame();
  GENPhase2 phaseHH(phh, mh12, mh22, hhframe, cosh, phih, 1);
  ANL4DVector ph1 = phaseHH.GetFourMomentum(0);
  ANL4DVector ph2 = phaseHH.GetFourMomentum(1);
  fP[0] = ph1;
  fM[0] = TMath::Sqrt(mh12);
  fP[1] = ph2;
  fM[1] = TMath::Sqrt(mh22);
  Double_t betah = phaseHH.GetBetaBar();
  if (betah <= 0.) return 0.;

  GENBranch &zbranch = * cmbranch.GetBranchPtr(1);
  Double_t cosf = zbranch.GetCosTheta();
  Double_t phif = zbranch.GetPhi     ();
  Double_t m32  = zbranch.GetM12();
  Double_t m42  = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m32, m42, cmframe, cosf, phif, 1);
  fP[2] = phaseZ.GetFourMomentum(0);
  fP[3] = phaseZ.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betaf = phaseZ.GetBetaBar();
  if (betaf <= 0.) return 0.;

  Double_t eb     = rs/2.;
  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
  Double_t s      = rs * rs;
#ifdef __DEBUG__
  cerr << " fP[0] = (" 
       << fP[0].E () << ","
       << fP[0].Px() << ","
       << fP[0].Py() << ","
       << fP[0].Pz() << ")" << endl;
  cerr << " fP[1] = (" 
       << fP[1].E () << ","
       << fP[1].Px() << ","
       << fP[1].Py() << ","
       << fP[1].Pz() << ")" << endl;
  cerr << " fP[2] = (" 
       << fP[2].E () << ","
       << fP[2].Px() << ","
       << fP[2].Py() << ","
       << fP[2].Pz() << ")" << endl;
  cerr << " fP[3] = (" 
       << fP[3].E () << ","
       << fP[3].Px() << ","
       << fP[3].Py() << ","
       << fP[3].Pz() << ")" << endl;
  ANL4DVector qz = fP[2] + fP[3];
  cerr << " qz = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector pcm = fP[0] + fP[1] + qz;
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
#endif

  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared(cmbranch);

  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1./2.;                         // identical particle factor
  Double_t dPhase = kFact * betax * betaf * betah; // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor
  Double_t spin   = 1./2.;                         // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ZHHBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

#ifdef __NODECAY__
  Double_t q2z    = cmbranch.GetM22();
  Double_t rq2z   = TMath::Sqrt(q2z);
  Double_t gamz   = fZBosonPtr->GetWidth();
           amp2  *= (2*kM_z*gamz)/(TMath::Power((rq2z - kM_z)*(rq2z + kM_z),2)
			                          + TMath::Power(kM_z*gamz,2));
  amp2 *= k8Pi; // branch factor correction
#endif
  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ZHHBases::FullAmplitude()
{
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf     = f3Ptr->GetCharge();
   Double_t t3f    = f3Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  h1(fP[0]);
   HELScalar  h2(fP[1]);

   HELFermion f (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
   HELFermion fb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoZHH(em, ep, h1, h2, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZHH()
// --------------------------
Complex_t ZHHBases::AmpEEtoZHH(const HELFermion &em,
                               const HELFermion &ep,
                               const HELScalar  &h1,
                               const HELScalar  &h2,
                               const HELVector  &zf)
{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamz  = fZBosonPtr->GetWidth();

   //--------------------------------------------------
   // Calculate Amplitudes
   //--------------------------------------------------
#ifdef __PHASESPACE__
   //---------------------------
   // Just Phase Space
   //---------------------------
   Complex_t amp = 1;
#else
   //---------------------------
   // Higgs Production Amplitude
   //---------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);

   Double_t v     = 2.*kM_w/kGw;
   Double_t ghhh  = -TMath::Power(fMass,2)/v*3.; 
   Double_t gzzh  = kGz*kM_z;
   Double_t gzzhh = kGz*kGz/2.; 

   HELScalar hh(h1, h2, ghhh, fMass, 0.);
   HELVertex amp1(zs, zf, hh, gzzh);         // HHH self-coupling

   HELVertex amp2(zs, zf, h1, h2, gzzhh);    // ZZHH 4-point

   HELVector vz1(zf, h1, gzzh, kM_z, gamz);
   HELVertex amp3(zs, vz1, h2, gzzh);        // double H-strahlung

   HELVector vz2(zf, h2, gzzh, kM_z, gamz);
   HELVertex amp4(zs, vz2, h1, gzzh);        // double H-strahlung

   Complex_t amp  = amp1 + amp2 + amp3 + amp4;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZHHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZHHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZHHBases.BeamstrahlungFilename","trc500"));
    ins >> bsfilename;
    
    TFile *fBeamFile = new TFile((bsfiledir+bsfilename+".root").data());
    if (!fBeamFile) {
      cerr << " Unable to open a file for beamstrahlung" << endl;
      abort();
    }
    
  // --------------------------------------------
  //  Initialize beam generator
  // --------------------------------------------

    fBM = (JSFBeamGenerationCain*)fBeamFile->Get(bsfilename.data());
    fBM->SetIBParameters(fBeamWidth);

    fBM->MakeBSMap();
    fBM->Print(); 

    cout << " Nominal Energy (beamstrahlung) = " << fBM->GetNominalEnergy() 
         << " beam width     (beamstrahlung) ="  << fBM->GetIBWidth() 
         << endl;

    if (!(fBM->GetNominalEnergy() == (fEcmInit/2))) {
      cout << "Nominal energy from beamstrahung is " << fBM->GetNominalEnergy() 
           << " which is different from fEcm/2: "    << (fEcmInit/2)
           << endl;        
    } // check the energy homogeneity
  } // if beamstrahlung is on

  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZBosonPtr->GetMode(m); 
     if (mp && (m<fZModesLo || m>fZModesHi)) {
        mp->Lock();
     }
  }
  fZBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh" );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init( 4,    70.,   110.,       50, "Mz"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthF");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF"  );
  Xh_init( 7,     0., fEcmInit,     50, "Mhh"   );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthH");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiH"  );
  Xh_init(10,     0.,     2.,        2, "Helin ");
  Xh_init(11,     0.,     2.,        2, "Helot ");
  Xh_init(12,     0.,    12.,       12, "Z mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZHHBases::Userout()
{
  cout << "End of ZHHBases----------------------------------- "   << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremsstrahlung       = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void ZHHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][4] = {{0, 0, -1, +1},
                                           {0, 0, +1, -1}};
   Double_t helm = (1. - fPole)/2.;
   if (fHelCombInitial < helm) {
      fJCombI = 0;
   } else {
      fJCombI = 1;
   }
   fHelInitial[0] = kIHelComb[fJCombI][0];
   fHelInitial[1] = kIHelComb[fJCombI][1];
   fJCombF = (Int_t)(fHelCombFinal*kNf);
   fJCombF = TMath::Min(fJCombF, kNf-1);
   fHelFinal  [0] = kFHelComb[fJCombF][0];
   fHelFinal  [1] = kFHelComb[fJCombF][1];
   fHelFinal  [2] = kFHelComb[fJCombF][2];
   fHelFinal  [3] = kFHelComb[fJCombF][3];
   weight = kNf;
}
