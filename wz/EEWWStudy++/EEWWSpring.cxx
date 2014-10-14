//*****************************************************************************
//* =====================
//*  EEWWSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> EEWW generator
//*
//* (Update Record)
//*    2014/10/13  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "EEWWSpring.h"

#include <sstream>
#include <iomanip>
//#define __ZEROWIDTH__
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(EEWWSpring)
ClassImp(EEWWSpringBuf)
ClassImp(EEWWBases)

//-----------------------------------------------------------------------------
// ==============================
//  class EEWWSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
EEWWSpring::EEWWSpring(const char      *name,
                       const char      *title,
                             EEWWBases *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new EEWWSpringBuf("EEWWSpringBuf",
                                "EEWWSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new EEWWBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
EEWWSpring::~EEWWSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t EEWWSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    EEWWBases *bs = static_cast<EEWWBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> EEWWBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class EEWWSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t EEWWSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  EEWWBases    *bases   = static_cast<EEWWBases *>(
                          static_cast<EEWWSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 8;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // e-
  pv[1] = bases->fP[1]; // e+
  pv[2] = bases->fP[2]; // fub
  pv[3] = bases->fP[3]; // fd
  pv[4] = bases->fP[4]; // fu
  pv[5] = bases->fP[5]; // fdb
  pv[6] = pv[2] + pv[3];// W-
  pv[7] = pv[4] + pv[5];// W+

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------
  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fXi           = bases->GetXi();
  fEta1         = bases->GetEta1();
  fEta2         = bases->GetEta2();
  fPhi1         = bases->GetPhi1();
  fPhi2         = bases->GetPhi2();
  fQ2WW         = bases->GetQ2WW();
  fCosWm        = bases->GetCosWm ();
  fPhiWm        = bases->GetPhiWm ();
  fQ2Wm         = bases->GetQ2Wm  ();
  fCosWmF       = bases->GetCosWmF();
  fPhiWmF       = bases->GetPhiWmF();
  fQ2Wp         = bases->GetQ2Wp  ();
  fCosWpF       = bases->GetCosWpF();
  fPhiWpF       = bases->GetPhiWpF();
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);

  for (Int_t i=0; i<fNparton; i++) pv[i].Boost(boostv);

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
  Int_t    idw     = 24;                          // PDG code for W
  Int_t    ide     = 11;                          // PDG code for e
  Double_t me      = kM_e;                        // electron mass
  Double_t qe      = -1;                          // electron charge

  // W-
  Int_t    idf3    = -bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = -bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      =  bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    =  bases->fHelFinal[2];         // f3 helicity
  Double_t color3  =  bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = -bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = -bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      =  bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    =  bases->fHelFinal[3];         // f4 helicity

  Int_t    islevwm = color3 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfwm   = 2;                           // color flux id
  Double_t rq2wm   = pv[6].Mag();

  // W+
  Int_t    idf5    = bases->f5Ptr->GetPID   ();   // PDG code for f5
  Double_t chrg5   = bases->f5Ptr->GetCharge();   // f5 charge
  Double_t m5      = bases->f5Ptr->GetMass  ();   // f5 mass
  Int_t    hel5    = bases->fHelFinal[4];         // f5 helicity
  Double_t color5  = bases->f5Ptr->GetColor();    // color factor for f5

  Int_t    idf6    = bases->f6Ptr->GetPID   ();   // PDG code for f6
  Double_t chrg6   = bases->f6Ptr->GetCharge();   // f6 charge
  Double_t m6      = bases->f6Ptr->GetMass  ();   // f6 mass
  Int_t    hel6    = bases->fHelFinal[5];         // f6 helicity

  Int_t    islevwp = color5 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfwp   = 3;                           // color flux id
  Double_t rq2wp   = pv[7].Mag();
#ifdef __DEBUG__
  cerr << endl;
  ANL4DVector qcm;
  for (Int_t ip=0; ip<fNparton; ip++) {
    cerr << "pv[" << ip << "] = (" 
         <<  pv[ip].E() << ", "
         <<  pv[ip].Px() << ", "
         <<  pv[ip].Py() << ", "
         <<  pv[ip].Pz() << ") "
	 << "m = " << pv[ip].GetMass() << endl;
    if (ip<fNparton-2) qcm += pv[ip];
  }
  cerr << "qcm = ("
       <<  qcm.E() << ", "
       <<  qcm.Px() << ", "
       <<  qcm.Py() << ", "
       <<  qcm.Pz() << ") " << endl;
  cerr << "----" << endl;
#endif

  //                              No. PID   Mass  Charge   pv   Nd 1st Mom  hel   col  shower
  new (partons[0]) JSFSpringParton(1,  ide,   me,    qe, *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton(2, -ide,   me,   -qe, *qp[1], 0, 0,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton(3, -idw,rq2wm,   -1., *qp[6], 2, 5,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton(4,  idw,rq2wp,   +1., *qp[7], 2, 7,  0,    0,    0,      0);
  new (partons[4]) JSFSpringParton(5, idf3,   m3, chrg3, *qp[2], 0, 0,  3, hel3,icfwm,islevwm);
  new (partons[5]) JSFSpringParton(6, idf4,   m4, chrg4, *qp[3], 0, 0,  3, hel4,icfwm,islevwm);
  new (partons[6]) JSFSpringParton(7, idf5,   m5, chrg5, *qp[4], 0, 0,  4, hel5,icfwp,islevwp);
  new (partons[7]) JSFSpringParton(8, idf6,   m6, chrg6, *qp[5], 0, 0,  4, hel6,icfwp,islevwp);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class EEWWBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
EEWWBases::EEWWBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 125.),
           fWidth     ( 0.004),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPolem     (0.),
           fPolep     (0.),
           fWmModesLo ( 1),
           fWmModesHi (12),
           fWpModesLo ( 1),
           fWpModesHi (12),
	   fACC1      (0.05),
	   fACC2      (0.05),
	   fITMX1     (20),
	   fITMX2     (40),
           fWmBosonPtr( 0),
           fWpBosonPtr( 0),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2WW      (0.),
	   fCosWm     (0.),
	   fPhiWm     (0.),
	   fQ2Wm      (0.),
	   fCosWmF    (0.),
	   fPhiWmF    (0.),
	   fQ2Wp      (0.),
	   fCosWpF    (0.),
	   fPhiWpF    (0.),
           fWmModePtr (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fWpModePtr (0),
           f5Ptr      (0),
           f6Ptr      (0),
	   fSh1       (0.),
	   fCh1       (0.),
	   fSh2       (0.),
	   fCh2       (0.),
           fR_BW_m    (0),
           fR_BW_p    (0),
           fR_BS_m    (0),
           fR_BS_p    (0),
           fR_ISR_var (0),
           fR_ISR_side(0)
{
  //  Constructor of bases.  Default parameter should be initialized here
  //
  // --------------------------------------------
  //  Get parameters from jsf.conf, if specified
  // --------------------------------------------

  cout << "Init eewwbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("EEWWBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.WidthH","0.006")); // M_x [GeV]
  ins >> fWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.Polem","0."));       // electron polarization
  ins >> fPolem;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.Polep","0."));       // positron polarization
  ins >> fPolep;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.CosthWmRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.PhiWmOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.CosthWmFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.PhiWmFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.CosthWpFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.PhiWpFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.WmModesLo","1"));      // Wm decay mode lo
  ins >> fWmModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.WmModesHi","12"));     // Wm decay mode hi
  ins >> fWmModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.WpModesLo","1"));      // Wp decay mode lo
  ins >> fWpModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.WpModesHi","12"));     // Wp decay mode hi
  ins >> fWpModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEWWBases.NCALL","80000"));
  ins >> fNCALL;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Beamstrahlung and natural Ebeam spread
  //--
  if (fBeamStr == 1) {
    DefineVariable(fR_BW_m, 0., 1., 1, 1);
    DefineVariable(fR_BW_p, 0., 1., 1, 1);
    DefineVariable(fR_BS_m, 0., 1., 1, 1);
    DefineVariable(fR_BS_p, 0., 1., 1, 1);
  }

  //--
  //  ISR
  //--
  if (fISR==1) {
    DefineVariable(fR_ISR_var,  0., 1., 1, 1);
    DefineVariable(fR_ISR_side, 0., 1., 1, 1);
  }

  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);

  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1, q2ww
  //--
  DefineVariable(fXXi    , 0., 1., 1, 1);
  DefineVariable(fXEta1  , 0., 1., 1, 1);
  DefineVariable(fXEta2  , 0., 1., 1, 1);
  DefineVariable(fXQ2WW  , 0., 1., 1, 1);
  DefineVariable(fXPhi21 , 0., 1., 0, 1);
  DefineVariable(fXPhi1  , 0., 1., 0, 1);

  //--
  //  Q2 of Wm and Wp
  //--
  DefineVariable(fXQ2Wm  , 0., 1., 0, 1);
  DefineVariable(fXQ2Wp  , 0., 1., 0, 1);

  //--
  //  cos(theta) and phi
  //--
  fXL[1] *= TMath::Pi();
  fXU[1] *= TMath::Pi();
  fXL[3] *= TMath::Pi();
  fXU[3] *= TMath::Pi();
  fXL[5] *= TMath::Pi();
  fXU[5] *= TMath::Pi();

  DefineVariable(fCosWm , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhiWm , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosWmF, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiWmF, fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosWpF, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiWpF, fXL[5], fXU[5], 0, 1);

  //--
  //  Final states
  //--
  DefineVariable(fWmDecayMode   , 0., 1., 0, 1);
  DefineVariable(fWpDecayMode   , 0., 1., 0, 1);

  // --------------------------------------------
  //  Set Bases integration parameters
  // --------------------------------------------
  SetNoOfSample(fNCALL);

  SetTuneValue (1.5);
  SetIteration1(fACC1, fITMX1);
  SetIteration2(fACC2, fITMX2);
}
// --------------------------
//  D-tor
// --------------------------
EEWWBases::~EEWWBases()
{
  delete fWmBosonPtr;
  delete fWpBosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t EEWWBases::Func()
{
  //  Bases Integrand
  //

  Double_t bsWeight = 1.; // Jacobian factor
  Double_t eminus;        // E_e- after ISR and beamstrahlung
  Double_t eplus;         // E_e+ after ISR and beamstrahlung

  // --------------------------------------------
  //  Beamstrahlung
  // --------------------------------------------
  if (fBeamStr == 1) {
    bsWeight = fBM->GetWeight(fR_BW_m, fR_BW_p,
                              fR_BS_m, fR_BS_p,
                              eminus, eplus);
    eplus   *= fEcmInit/2.;
    eminus  *= fEcmInit/2.;
    fEcmIP   = 2.*TMath::Sqrt(eplus*eminus);     // ignore electron mass 
  } else {
    fEcmIP   = fEcmInit;      
    eplus    = fEcmInit/2.;
    eminus   = eplus;
  }
   
  // --------------------------------------------
  //  Initial State Radiation
  // --------------------------------------------

  static const Double_t kFisr = 2. * kAlpha0 / TMath::Pi();

  Double_t beta_isr = kFisr * (2.*(TMath::Log(fEcmIP/kM_e) - 1.));

  Double_t qed = 1.;
  if (fISR == 1) {
    qed     = (1. + 3.*beta_isr/4.)
              *(1. + kFisr * ( TMath::Pi()*TMath::Pi()/6. - 1./4. ));
    Double_t xisr    = TMath::Power(fR_ISR_var,1./beta_isr); // Ephoton / Ebeam 
             fEcmIP *= TMath::Sqrt(1. - xisr);               // reduced Ecm
    if (fR_ISR_side > 0.5) {
      eplus  *= 1. - xisr;
    } else {
      eminus *= 1. - xisr;
    }
  } // if Bremsstrahlung


  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung

  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  Double_t weight = 1;

  GENDecayMode *fWmModePtr = fWmBosonPtr->PickMode(fWmDecayMode, weight, fWmMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  GENDecayMode *fWpModePtr = fWpBosonPtr->PickMode(fWpDecayMode, weight, fWpMode);
  bsWeight *= weight;
  f5Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(0));
  f6Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(1));
  Double_t m5   = f5Ptr->GetMass();
  Double_t m6   = f6Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = kM_e; // Me
  Double_t m2   = kM_e; // Meb
  if (fEcmIP < m1 + m2 + m3 + m4 + m5 + m6) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average is taken care of here
  SelectHelicities(weight);
  if (weight == 0.) return 0.;
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  Double_t s      = fEcmIP*fEcmIP;
  Double_t rs     = fEcmIP;

  // Wm
  Double_t qwmmin = m3 + m4;
  Double_t qwmmax = rs - (m1 + m2 + m5 + m6);
#ifndef __ZEROWIDTH__
  fQ2Wm = fWmBosonPtr->GetQ2BW(qwmmin, qwmmax, fXQ2Wm, weight);
#else
  fQ2Wm = TMath::Power(fWmBosonPtr->GetMass(),2);
  weight = kPi*fWmBosonPtr->GetMass()*fWmBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2wm  = TMath::Sqrt(fQ2Wm);

  // Wp
  Double_t qwpmin = m5 + m6;
  Double_t qwpmax = rs - (m1 + m2 + rq2wm);
#ifndef __ZEROWIDTH__
  fQ2Wp = fWpBosonPtr->GetQ2BW(qwpmin, qwpmax, fXQ2Wp, weight);
#else
  fQ2Wp = TMath::Power(fWpBosonPtr->GetMass(),2);
  weight = kPi*fWpBosonPtr->GetMass()*fWpBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2wp  = TMath::Sqrt(fQ2Wp);

  // Q2 of WW
  Double_t qww2mn = TMath::Power(rq2wm + rq2wp, 2);
  Double_t qww2mx = TMath::Power(rs - (m1+m2),2);
  fQ2WW     = qww2mn + (qww2mx-qww2mn)*fXQ2WW;
  bsWeight *= qww2mx - qww2mn;
  Double_t rq2ww  = TMath::Sqrt(fQ2WW);


  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t xilo  = TMath::Log(rq2ww*(rq2ww+m1+m2)/s);
  Double_t xihi  = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);
  fXi       = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t rxi  = TMath::Exp(fXi);

  Double_t dm1  = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1  = (m1*m1/s);
  Double_t dm2  = dp1;
  Double_t dp2  = dm1;

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;
  fEta1     = etlo + (ethi-etlo)*fXEta1;
  bsWeight *= ethi - etlo;

           etlo = -TMath::Log((1.+dm2)/dp2)/2.;
           ethi =  TMath::Log((1.+dp2)/dm2)/2.;
  fEta2     = etlo + (ethi-etlo)*fXEta2;
  bsWeight *= ethi - etlo;

  fPhi1  = fXPhi1*k2Pi;
  fPhi2  = fXPhi21*k2Pi + fPhi1;
  fPhi2  = fPhi2 > k2Pi ? fPhi2 - k2Pi : fPhi2;
  bsWeight *= k2Pi*k2Pi;

  fM[0] = m1;
  fM[1] = m2;

  GENBranch wmbranch(fQ2Wm, fCosWmF, fPhiWmF, m3*m3    , m4*m4);
  GENBranch wpbranch(fQ2Wp, fCosWpF, fPhiWpF, m5*m5    , m6*m6);
  GENBranch wwbranch(fQ2WW, fCosWm , fPhiWm , &wmbranch, &wpbranch);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX(wwbranch) * qed;

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fXi               , (bsWeight*sigma));
  Xh_fill( 3, fEta1             , (bsWeight*sigma));
  Xh_fill( 4, fEta2             , (bsWeight*sigma));
  Xh_fill( 5, fPhi1             , (bsWeight*sigma));
  Xh_fill( 6, fPhi2             , (bsWeight*sigma));
  Xh_fill( 7, k2Pi*fXPhi21      , (bsWeight*sigma));
  Xh_fill( 8, TMath::Sqrt(fQ2WW), (bsWeight*sigma));
  Xh_fill( 9, fCosWm            , (bsWeight*sigma));
  Xh_fill(10, fPhiWm            , (bsWeight*sigma));
  Xh_fill(11, TMath::Sqrt(fQ2Wm), (bsWeight*sigma));
  Xh_fill(12, fCosWmF           , (bsWeight*sigma));
  Xh_fill(13, fPhiWmF           , (bsWeight*sigma));
  Xh_fill(14, TMath::Sqrt(fQ2Wp), (bsWeight*sigma));
  Xh_fill(15, fCosWpF           , (bsWeight*sigma));
  Xh_fill(16, fPhiWpF           , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fWmMode , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fWpMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t EEWWBases::DSigmaDX(GENBranch &wwbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs  = fEcmIP;
  Double_t s   = rs*rs; 
  Double_t eb  = rs/2.;

  Double_t qww2 = wwbranch.GetQ2();  // q^2_{WW}
  Double_t qww  = TMath::Sqrt(qww2);
  Double_t q2wm = wwbranch.GetM12(); // Wm
  Double_t q2wp = wwbranch.GetM22(); // Wp

  Double_t rxi = TMath::Exp(fXi);
  weight *= rxi;

  Double_t m1  = fM[0]; // Me
  Double_t m2  = fM[1]; // Meb
  Double_t dm1 = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1 = (m1*m1/s);
  Double_t dm2 = dp1;
  Double_t dp2 = dm1;

  Double_t exp2et1 = TMath::Exp(2.*fEta1);
  Double_t sh1     = TMath::Min(TMath::Sqrt((1.+dp1+dm1)/(1.+exp2et1)   - dm1), 1.);
  Double_t ch1     = TMath::Min(TMath::Sqrt((1.+dp1+dm1)/(1.+1./exp2et1)- dp1), 1.);
  Double_t sn1     = 2.*sh1*ch1;
  Double_t cs1     = (1.+dp1+dm1)*TMath::TanH(fEta1) - dp1 + dm1;
  Double_t fi1     = fPhi1;
  weight *= 4.* (1.+dp1+dm1)/((1.+exp2et1)*(1+1/exp2et1));

  Double_t exp2et2 = TMath::Exp(2.*fEta2);
  Double_t sh2     = TMath::Min(TMath::Sqrt((1.+dp2+dm2)/(1.+exp2et2)   - dm2), 1.);
  Double_t ch2     = TMath::Min(TMath::Sqrt((1.+dp2+dm2)/(1.+1./exp2et2)- dp2), 1.);
  Double_t sn2     = 2.*sh2*ch2;
  Double_t cs2     = (1.+dp2+dm2)*TMath::TanH(fEta2) - dp2 + dm2;
  Double_t fi2     = fPhi2;
  weight *= 4.* (1.+dp2+dm2)/((1.+exp2et2)*(1+1/exp2et2));

  Double_t cs12    = cs1*cs2 + sn1*sn2*TMath::Cos(fi2-fi1);

  Double_t x1, x2;
  if (fEta1 > -fEta2) {
    Double_t rximn = (qww-m1+m2)*(qww+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - qww2/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - qww2/s);
  } else {
    Double_t rximn = (qww+m1-m2)*(qww+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - qww2/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - qww2/s);
  }

  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

  Double_t e1 = eb*x1;
  if (e1 < m1) return 0.;
  Double_t ap1 = TMath::Sqrt((e1-m1)*(e1+m1));
  fP[0].SetXYZT(ap1*sn1*TMath::Cos(fi1), ap1*sn1*TMath::Sin(fi1), ap1*cs1, e1);

  Double_t e2 = eb*x2;
  if (e2 < m2) return 0.;
  Double_t ap2 = TMath::Sqrt((e2-m2)*(e2+m2));
  fP[1].SetXYZT(ap2*sn2*TMath::Cos(fi2), ap2*sn2*TMath::Sin(fi2), ap2*cs2, e2);

  ANL4DVector qcm(rs, 0.,0.,0.);
  ANL4DVector qvww = qcm - fP[0] - fP[1];
  GENFrame    cmframe;
  GENPhase2   phaseWW(qvww, q2wm, q2wp, cmframe, fCosWm, fPhiWm, 1);
  ANL4DVector pwm   = phaseWW.GetFourMomentum(0);
  ANL4DVector pwp   = phaseWW.GetFourMomentum(1);
  Double_t    betaz = phaseWW.GetBetaBar();
  if (betaz <= 0.) return 0.;
  weight *= betaz;

  GENBranch &wmbranch = *wwbranch.GetBranchPtr(0);
  Double_t coswmf  = wmbranch.GetCosTheta();
  Double_t phiwmf  = wmbranch.GetPhi     ();
  Double_t m32     = wmbranch.GetM12();
  Double_t m42     = wmbranch.GetM22();
  GENFrame wwframe = phaseWW.GetFrame(1);
  GENPhase2 phaseWm(pwm, m32, m42, wwframe, coswmf, phiwmf, 1);
  fP[2] = phaseWm.GetFourMomentum(0);
  fP[3] = phaseWm.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betawmf = phaseWm.GetBetaBar();
  if (betawmf <= 0.) return 0.;
  weight *= betawmf;

  GENBranch &wpbranch = *wwbranch.GetBranchPtr(1);
  Double_t coswpf = wpbranch.GetCosTheta();
  Double_t phiwpf = wpbranch.GetPhi     ();
  Double_t m52    = wpbranch.GetM12();
  Double_t m62    = wpbranch.GetM22();
  GENPhase2 phaseWp(pwp, m52, m62, wwframe, coswpf, phiwpf, 1);
  fP[4] = phaseWp.GetFourMomentum(0);
  fP[5] = phaseWp.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(m52);
  fM[5] = TMath::Sqrt(m62);
  Double_t betawpf = phaseWp.GetBetaBar();
  if (betawpf <= 0.) return 0.;
  weight *= betawpf;

  fSh1 = sh1;
  fCh1 = ch1;
  fSh2 = sh2;
  fCh2 = ch2;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
#ifdef __DEBUG__
  for (Int_t i=0; i<6; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qvwm = fP[2] + fP[3];
  cerr << " pwm = (" 
       << qvwm.E () << ","
       << qvwm.Px() << ","
       << qvwm.Py() << ","
       << qvwm.Pz() << ")" << endl;
  ANL4DVector qvwp = fP[4] + fP[5];
  cerr << " pwp = (" 
       << qvwp.E () << ","
       << qvwp.Px() << ","
       << qvwp.Py() << ","
       << qvwp.Pz() << ")" << endl;
  ANL4DVector pcm = qvwm + qvwp + fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
  cerr << " wmmass = " << qvwm.GetMass() << endl;
  cerr << " wpmass = " << qvwp.GetMass() << endl;
#endif
  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();
  
  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 5;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * weight;                // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor

  Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                       // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t EEWWBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor() * f5Ptr->GetColor();
  Int_t     ig3   = f3Ptr->GetGenNo() - 1;
  Int_t     ig4   = f4Ptr->GetGenNo() - 1;
  Int_t     ig5   = f5Ptr->GetGenNo() - 1;
  Int_t     ig6   = f6Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig3][ig4]*kVkm[ig5][ig6],2);

  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t EEWWBases::FullAmplitude()
{
   Double_t mw     = fWmBosonPtr->GetMass();
   Double_t gamw   = fWmBosonPtr->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em (fK[0], kM_e,  fHelInitial[0], +1, kIsIncoming);
   HELFermion ep (fK[1], kM_e,  fHelInitial[1], -1, kIsOutgoing);

   HELFermion emf(fP[0], fM[0], fHelFinal  [0], +1, kIsOutgoing); // final e-
   HELFermion epf(fP[1], fM[1], fHelFinal  [1], -1, kIsIncoming); // final e+

   HELFermion f3b(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming); // fubar
   HELFermion f4 (fP[3], fM[3], fHelFinal [3], +1, kIsOutgoing); // fd
   HELVector  wm(f3b, f4, glw, grw, mw, gamw);                 // W-

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // fu
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // fdbar
   HELVector  wp(f6b, f5, glw, grw, mw, gamw);                 // W+

   Complex_t amp = AmpEEtoEEWW(em, ep, emf, epf, wm, wp);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoEEWW()
// --------------------------
Complex_t EEWWBases::AmpEEtoEEWW(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &emf,
                                 const HELFermion &epf,
                                 const HELVector  &wm,
                                 const HELVector  &wp)
{
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
   // WW Production Amplitudes
   //---------------------------
   Double_t  mn    = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamn  = 0.;
   Double_t  game  = 0.;

   Double_t  mw    = fWmBosonPtr->GetMass();
   Double_t  gamw  = fWmBosonPtr->GetWidth();
   Double_t  mz    = fZBosonPtr->GetMass();
   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  mh    = fMass;
   Double_t  gamh  = fWidth;

   Double_t  glw   = -kGw*kSqh;
   Double_t  grw   = 0.;

   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);
   Double_t  glae  = -kGe*qe;
   Double_t  grae  = -kGe*qe;

   Double_t  t3n   = +1./2.;
   Double_t  glzn  = -kGz*t3n;
   Double_t  grzn  = 0.;

   Double_t  gwwa  = kGe;     
   Double_t  gwwz  = kGw*kCosW;
   Double_t  gw    = kGw;

   Double_t  gwwh  = kGw*mw;
   Double_t  gzzh  = kGz*mz;

   Double_t  ebm   = TMath::Abs(em. GetFourMomentum()(0));
   Double_t  e1    = TMath::Abs(emf.GetFourMomentum()(0));
   Double_t  e2    = TMath::Abs(epf.GetFourMomentum()(0));

   //-------------------------
   // Calculate internal lines
   //-------------------------
   // virtual almost real A/Z
   // (e- leg)
   HELVector wrk01(ebm, e1, fSh1, fCh1, fPhi1, 
                   fHelInitial[0], fHelFinal[0],+1, kGe, me);
   HELVector wrk02(em, emf, glze, grze, mz, gamz);

   // virtual almost real A/Z
   // (e+ leg)
   HELVector wrk03(ebm, e2, fSh2, fCh2, fPhi2, 
                   fHelInitial[1], fHelFinal[1],-1, kGe, me);
   HELVector wrk04(epf, ep, glze, grze, mz, gamz);

   HELFermion wrk05(emf  , wp   , glw  , grw  , mn, gamn);
   HELVector  wrk06(em   , wrk05, glw  , grw  , mw, gamw);

   HELFermion wrk07(epf  , wm   , glw  , grw  , mn, gamn);
   HELVector  wrk08(wrk07, ep   , glw  , grw  , mw, gamw);

   HELFermion wrk09(em   , wm   , glw  , grw  , mn, gamn);
   HELVector  wrk10(wrk09, emf  , glw  , grw  , mw, gamw);

   HELFermion wrk11(ep   , wp   , glw  , grw  , mn, gamn);
   HELVector  wrk12(epf  , wrk11, glw  , grw  , mw, gamw);

   HELFermion wrk13(em   , wrk03, glae , grae , me, game);
   HELFermion tmp01(em   , wrk04, glze , grze , me, game);
   wrk13[0] = wrk13[0] + tmp01[0];
   wrk13[1] = wrk13[1] + tmp01[1];
   wrk13[2] = wrk13[2] + tmp01[2];
   wrk13[3] = wrk13[3] + tmp01[3];

   HELFermion wrk14(ep   , wrk01, glae , grae , me, game);
   HELFermion tmp02(ep   , wrk02, glze , grze , me, game);
   wrk14[0] = wrk14[0] + tmp02[0];
   wrk14[1] = wrk14[1] + tmp02[1];
   wrk14[2] = wrk14[2] + tmp02[2];
   wrk14[3] = wrk14[3] + tmp02[3];

   HELFermion wrk15(emf  , wrk03, glae , grae , me, game);
   HELFermion tmp03(emf  , wrk04, glze , grze , me, game);
   wrk15[0] = wrk15[0] + tmp03[0];
   wrk15[1] = wrk15[1] + tmp03[1];
   wrk15[2] = wrk15[2] + tmp03[2];
   wrk15[3] = wrk15[3] + tmp03[3];

   HELFermion wrk16(epf  , wrk01, glae , grae , me, game);
   HELFermion tmp04(epf  , wrk02, glze , grze , me, game);
   wrk16[0] = wrk16[0] + tmp04[0];
   wrk16[1] = wrk16[1] + tmp04[1];
   wrk16[2] = wrk16[2] + tmp04[2];
   wrk16[3] = wrk16[3] + tmp04[3];

   HELVector wrk17 (wrk01[0]*(gwwa/gw) + wrk02[0]*(gwwz/gw),
                    wrk01[1]*(gwwa/gw) + wrk02[1]*(gwwz/gw),
                    wrk01[2]*(gwwa/gw) + wrk02[2]*(gwwz/gw),
                    wrk01[3]*(gwwa/gw) + wrk02[3]*(gwwz/gw),
                    wrk01.GetFourMomentum());

   HELVector wrk18 (wrk03[0]*(gwwa/gw) + wrk04[0]*(gwwz/gw),
                    wrk03[1]*(gwwa/gw) + wrk04[1]*(gwwz/gw),
                    wrk03[2]*(gwwa/gw) + wrk04[2]*(gwwz/gw),
                    wrk03[3]*(gwwa/gw) + wrk04[3]*(gwwz/gw),
                    wrk03.GetFourMomentum());

   //---------------------------
   // Now calculate amplitudes
   //---------------------------
   // Non-fusion diagrams
   // ( 1)
   HELVertex amp01 (wrk07, ep   , wrk06, glw , grw );
   // ( 2)
   HELVertex amp02 (wrk07, wrk11, wrk02, glzn, grzn);
   // ( 3)
   HELVertex amp03 (wrk09, wrk05, wrk04, glzn, grzn);
   // ( 4)
   HELVertex amp04 (epf  , wrk11, wrk10, glw , grw );
   // ( 5)
   HELVector tmp05 (wrk13, emf  , glae, grae, glze, grze, mz, gamz);
   HELVertex amp05 (wm, wp, tmp05, gw);
   // ( 6)
   HELVector tmp06 (em   , wrk15, glae, grae, glze, grze, mz, gamz);
   HELVertex amp06 (wm, wp, tmp06, gw);
   // ( 7)
   HELVector tmp07 (wrk16, ep   , glae, grae, glze, grze, mz, gamz);
   HELVertex amp07 (wm, wp, tmp07, gw);
   // ( 8)
   HELVector tmp08 (epf  , wrk14, glae, grae, glze, grze, mz, gamz);
   HELVertex amp08 (wm, wp, tmp08, gw);
   // ( 9)
   HELVertex amp09 (wrk13, wrk05, wm, glw, grw);
   // (10)
   HELVertex amp10 (wrk09, wrk15, wp, glw, grw);
   // (11)
   HELVertex amp11 (wrk16, wrk11, wm, glw, grw);
   // (12)
   HELVertex amp12 (wrk07, wrk14, wp, glw, grw);
   // (13)
   HELVertex amp13 (wm   , wrk06, wrk18, gw);
   // (14)
   HELVertex amp14 (wm   , wrk12, wrk17, gw);
   // (15)
   HELVertex amp15 (wrk10, wp   , wrk18, gw);
   // (16)
   HELVertex amp16 (wrk08, wp   , wrk17, gw);

   //---------------------------
   // Fusion diagrams
   // (17-19)
   HELVertex amp17 (wm, wrk17, wp, wrk18, gw, gw, mw, gamw,kTRUE);

   // Higgs
   HELScalar tmp20 (wm, wp, gwwh, mh, gamh);
   HELVertex amp18 (wrk02, wrk04, tmp20, gzzh);

   //---------------------------
   // Sum up amplitudes
   //---------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 
	         + amp05 + amp06 + amp07 + amp08 
		 + amp09 + amp10 + amp11 + amp12
		 + amp13 + amp14 + amp15 + amp16 
		 + amp17 
		 + amp18
		 ;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void EEWWBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("EEWWBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("EEWWBases.BeamstrahlungFilename","trc500"));
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
  if (!fWmBosonPtr) fWmBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWmBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWmBosonPtr->GetMode(m); 
     if (mp && (m<fWmModesLo || m>fWmModesHi)) {
        mp->Lock();
     }
  }
  fWmBosonPtr->DebugPrint();

  if (!fWpBosonPtr) fWpBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWpBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWpBosonPtr->GetMode(m); 
     if (mp && (m<fWpModesLo || m>fWpModesHi)) {
        mp->Lock();
     }
  }
  fWpBosonPtr->DebugPrint();

  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  fZBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs   = fEcmInit;
  Double_t s    = rs*rs;

  Double_t qww  = 2*kM_z;
  Double_t m1   = kM_e;
  Double_t m2   = kM_e;

  Double_t xilo = TMath::Log(qww*(qww+m1+m2)/s);
  Double_t xihi = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);

  Double_t rxi  = TMath::Exp(xilo);
  Double_t dm1  = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1  = (m1*m1/s);

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;
  Double_t qwwlo = 0.; //2*mw - 40.;
  Double_t qwwhi = rs;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,  qwwlo,  qwwhi,       50, "qww"   );
  Xh_init( 9,    -1.,    +1.,       50, "coswm" );
  Xh_init(10,     0.,   k2Pi,       50, "phiwm" );
  Xh_init(11,    60.,   100.,       50, "Mwm"   );
  Xh_init(12,    -1.,    +1.,       50, "CosF3" );
  Xh_init(13,     0.,   k2Pi,       50, "PhiF3" );
  Xh_init(14,    60.,   100.,       50, "Mwp"   );
  Xh_init(15,    -1.,    +1.,       50, "CosF5" );
  Xh_init(16,     0.,   k2Pi,       50, "PhiF5" );
  Xh_init(17,     0.,     4.,        4, "Helin ");
  Xh_init(18,     0.,     4.,        4, "Helot ");
  Xh_init(19,     0.,    12.,       12, "Wm mode");
  Xh_init(20,     0.,    12.,       12, "Wp mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void EEWWBases::Userout()
{
  cout << "End of EEWWBases----------------------------------- "  << endl
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
void EEWWBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 4;
   static const Int_t kIHelComb[kNi][2] = {{-1, -1},
                                           {-1, +1},
                                           {+1, -1},
                                           {+1, +1}};
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][6] = {{-1, -1, +1, -1, -1, +1},
                                           {-1, +1, +1, -1, -1, +1},
                                           {+1, -1, +1, -1, -1, +1},
                                           {+1, +1, +1, -1, -1, +1}};

   Double_t helm = (1. - fPolem)/2.;
   Double_t help = (1. - fPolep)/2.;
   if (fHelCombInitial < helm || helm == 1.) {
      Double_t helcombi = fHelCombInitial/helm;
      if (helcombi < help) fJCombI = 0;
      else                 fJCombI = 1;
   } else {
      Double_t helcombi = (fHelCombInitial-helm)/(1.-helm);
      if (helcombi < help) fJCombI = 2;
      else                 fJCombI = 3;
   }
   fHelInitial[0] = kIHelComb[fJCombI][0];
   fHelInitial[1] = kIHelComb[fJCombI][1];
   fJCombF = (Int_t)(fHelCombFinal*kNf);
   fJCombF = TMath::Min(fJCombF, kNf-1);

   fHelFinal  [0] = kFHelComb[fJCombF][0];
   fHelFinal  [1] = kFHelComb[fJCombF][1];
   fHelFinal  [2] = kFHelComb[fJCombF][2];
   fHelFinal  [3] = kFHelComb[fJCombF][3];
   fHelFinal  [4] = kFHelComb[fJCombF][4];
   fHelFinal  [5] = kFHelComb[fJCombF][5];
   weight = kNf;
}
