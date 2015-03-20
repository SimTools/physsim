//*****************************************************************************
//* =====================
//*  NNWWSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> NNWW generator
//*
//* (Update Record)
//*    2014/09/19  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "NNWWSpring.h"
#include "HBoson.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__
#define TEMP_H

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(NNWWSpring)
ClassImp(NNWWSpringBuf)
ClassImp(NNWWBases)

//-----------------------------------------------------------------------------
// ==============================
//  class NNWWSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNWWSpring::NNWWSpring(const char      *name,
                       const char      *title,
                             NNWWBases *bases)
           : JSFSpring(name, title, bases)
{
  fEventBuf = new NNWWSpringBuf("NNWWSpringBuf",
                                "NNWWSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new NNWWBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
NNWWSpring::~NNWWSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t NNWWSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    NNWWBases *bs = static_cast<NNWWBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> NNWWBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class NNWWSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t NNWWSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  NNWWBases    *bases   = static_cast<NNWWBases *>(
                          static_cast<NNWWSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 8;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // ne
  pv[1] = bases->fP[1]; // neb
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
  Int_t    idw     = 24;                          // PDG code for W+
  Int_t    idne    = 12;                          // PDG code for nu_e

  // W-
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f3 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f4 helicity

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
    qcm += pv[ip];
  }
  cerr << "qcm = ("
       <<  qcm.E() << ", "
       <<  qcm.Px() << ", "
       <<  qcm.Py() << ", "
       <<  qcm.Pz() << ") " << endl;
  cerr << "----" << endl;
#endif

  //                              No. PID   Mass  Charge   pv    Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, idne,   0.,    0., *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton(2,-idne,   0.,    0., *qp[1], 0, 0,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton(3, -idw,rq2wm,   -1., *qp[6], 2, 5,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton(4,  idw,rq2wp,   +1., *qp[7], 2, 7,  0,    0,    0,      0);
  new (partons[4]) JSFSpringParton(5,-idf3,   m3,-chrg3, *qp[2], 0, 0,  3, hel3,icfwm,islevwm);
  new (partons[5]) JSFSpringParton(6,-idf4,   m4,-chrg4, *qp[3], 0, 0,  3, hel4,icfwm,islevwm);
  new (partons[6]) JSFSpringParton(7, idf5,   m5, chrg5, *qp[4], 0, 0,  4, hel5,icfwp,islevwp);
  new (partons[7]) JSFSpringParton(8, idf6,   m6, chrg6, *qp[5], 0, 0,  4, hel6,icfwp,islevwp);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class NNWWBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNWWBases::NNWWBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 125.),
           fLambda     (1000.),
           fA          (   0.),
           fB          (   0.),
           fBtilde     (   0.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fWmModesLo ( 1),
           fWmModesHi (12),
           fWpModesLo ( 1),
           fWpModesHi (12),
	   fNCALL     (80000),
	   fACC1      (0.05),
	   fACC2      (0.05),
	   fITMX1     (20),
	   fITMX2     (40),
           fHBosonPtr  ( 0),
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

  cout << "Init nnhbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("NNWWBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Lambda","1000.")); 	 // Lambda [GeV]
  ins >> fLambda;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.A","0.")); 	 // a
  ins >> fA;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.B","0.")); 	 // b
  ins >> fB;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Btilde","0.")); 	 // btilde
  ins >> fBtilde;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.WmModesLo","1"));      // W- decay mode lo
  ins >> fWmModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.WmModesHi","12"));     // W- decay mode hi
  ins >> fWmModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.WpModesLo","1"));      // W+ decay mode lo
  ins >> fWpModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.WpModesHi","12"));     // W+ decay mode hi
  ins >> fWpModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNWWBases.NCALL","80000"));
  ins >> fNCALL;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fWmDecayMode   , 0., 1., 0, 1);
  DefineVariable(fWpDecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2WW         , 0., 1., 1, 1);
  DefineVariable(fXQ2Wm         , 0., 1., 0, 1);
  DefineVariable(fXQ2Wp         , 0., 1., 0, 1);
  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1
  //--
  DefineVariable(fXXi    , 0., 1., 1, 1);
  DefineVariable(fXEta1  , 0., 1., 1, 1);
  DefineVariable(fXEta2  , 0., 1., 1, 1);
  DefineVariable(fXPhi1  , 0., 1., 1, 1);
  DefineVariable(fXPhi21 , 0., 1., 1, 1);
  DefineVariable(fXQ2WW  , 0., 1., 1, 1);
  DefineVariable(fXCosWm , 0., 1., 0, 1);
  DefineVariable(fXPhiWm , 0., 1., 0, 1);
  DefineVariable(fXQ2Wm  , 0., 1., 1, 1);
  DefineVariable(fXCosWmF, 0., 1., 0, 1);
  DefineVariable(fXPhiWmF, 0., 1., 0, 1);
  DefineVariable(fXQ2Wp  , 0., 1., 1, 1);
  DefineVariable(fXCosWpF, 0., 1., 0, 1);
  DefineVariable(fXPhiWpF, 0., 1., 0, 1);

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
NNWWBases::~NNWWBases()
{
  delete fHBosonPtr;
  delete fWmBosonPtr;
  delete fWpBosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t NNWWBases::Func()
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
  Double_t m1   = 0.; // Mne
  Double_t m2   = 0.; // Mneb
  if (fEcmIP < m1 + m2 + m3 + m4 + m5 + m6) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average for e- is taken
  //  care of here
  SelectHelicities(weight);
  if (weight == 0.) return 0.;
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  Double_t s      = fEcmIP*fEcmIP;
  Double_t rs     = fEcmIP;

#ifdef HIGGS_ONLY
  // H
#if 0
  // BW line shape
  Double_t qwwmn = m3 + m4 + m5 + m6;
  Double_t qwwmx = rs - (m1+m2);
  fQ2WW  = fHBosonPtr->GetQ2BW(qwwmn, qwwmx, fXQ2WW, weight);
  Double_t qww   = TMath::Sqrt(fQ2WW);
#else
  // Zero width approx.
  fQ2WW  = TMath::Power(fHBosonPtr->GetMass(),2);
  weight = kPi*fHBosonPtr->GetMass()*fHBosonPtr->GetWidth();
  Double_t qww   = fMass;
#endif
  bsWeight *= weight;
#endif
  // W-
  Double_t qwmmin = m3 + m4;
#ifndef HIGGS_ONLY
  Double_t qwmmax = rs - (m1 + m2 + m5 + m6);
#else
  Double_t qwmmax = qww - (m5 + m6);
#endif
#ifndef __ZEROWIDTH__
  fQ2Wm = fWmBosonPtr->GetQ2BW(qwmmin, qwmmax, fXQ2Wm, weight);
#else
  fQ2Wm = TMath::Power(fWmBosonPtr->GetMass(),2);
  weight = kPi*fWmBosonPtr->GetMass()*fWmBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // W+
  Double_t rq2wm  = TMath::Sqrt(fQ2Wm);
  Double_t qwpmin = m5 + m6;
#ifndef HIGGS_ONLY
  Double_t qwpmax = rs - (rq2wm + m1 + m2);
#else
  Double_t qwpmax = qww - rq2wm;
#endif
#ifndef __ZEROWIDTH__
  fQ2Wp = fWpBosonPtr->GetQ2BW(qwpmin, qwpmax, fXQ2Wp, weight);
#else
  fQ2Wp = TMath::Power(fWpBosonPtr->GetMass(),2);
  weight = kPi*fWpBosonPtr->GetMass()*fWpBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2wp  = TMath::Sqrt(fQ2Wp);

#ifndef HIGGS_ONLY
  // H
#if 0
  Double_t qww2mn = TMath::Power(rq2wm + rq2wp,2);
  Double_t qww2mx = TMath::Power(rs - (m1+m2),2);
  fQ2WW     = qww2mn + (qww2mx-qww2mn)*fXQ2WW;
  bsWeight *= qww2mx - qww2mn;
#else
  Double_t qwwmn = rq2wm + rq2wp;
  Double_t qwwmx = rs - (m1+m2);
  fQ2WW     = fHBosonPtr->GetQ2BW(qwwmn, qwwmx, fXQ2WW, weight);
  bsWeight *= weight;
#endif
  Double_t qww   = TMath::Sqrt(fQ2WW);
#endif

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t xilo  = TMath::Log(qww*(qww+m1+m2)/s);
  Double_t xihi  = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);
  fXi       = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t mw   = fWmBosonPtr->GetMass();
  Double_t dm1  = mw*mw/s;
  Double_t dp1  = 1.;
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

  fPhi1  = k2Pi * fXPhi1;
  fPhi2  = k2Pi * fXPhi21 + fPhi1;
  fPhi2  = fPhi2 > k2Pi ? fPhi2 - k2Pi : fPhi2;
  bsWeight *= k2Pi*k2Pi;

  fCosWm = -1. + 2.*fXCosWm;
  fPhiWm = k2Pi * fXPhiWm;
  bsWeight *= 2.*k2Pi;

  fCosWmF = -1. + 2.*fXCosWmF;
  fPhiWmF = k2Pi * fXPhiWmF;
  bsWeight *= 2.*k2Pi;

  fCosWpF = -1. + 2.*fXCosWpF;
  fPhiWpF = k2Pi * fXPhiWpF;
  bsWeight *= 2.*k2Pi;

  fM[0] = m1;
  fM[1] = m2;

  GENBranch wmbranch(fQ2Wm, fCosWmF, fPhiWmF, m3*m3    , m4*m4);
  GENBranch wpbranch(fQ2Wp, fCosWpF, fPhiWpF, m5*m5    , m6*m6);
  GENBranch wwbranch(fQ2WW, fCosWm , fPhiWm , &wmbranch, &wpbranch);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX(wwbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------
#ifdef TEMP_H
  H1Fill("hMww"     , TMath::Sqrt(fQ2WW), (bsWeight*sigma));
  H1Fill("hCosWm"   , fCosWm            , (bsWeight*sigma));
  H1Fill("hRSH"     , fEcmIP            , (bsWeight*sigma));
#endif

#if 1
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
#endif

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t NNWWBases::DSigmaDX(GENBranch &wwbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs   = fEcmIP;
  Double_t s    = rs*rs; 
  Double_t eb   = rs/2.;

  Double_t qww2 = wwbranch.GetQ2();  // q^2_{WW}
  Double_t qww  = TMath::Sqrt(qww2);
  Double_t q2wm = wwbranch.GetM12(); // W-
  Double_t q2wp = wwbranch.GetM22(); // W+

  Double_t rxi  = TMath::Exp(fXi);
  weight *= rxi;

  Double_t m1  = fM[0]; // Mne
  Double_t m2  = fM[1]; // Mneb
  Double_t mw  = fWmBosonPtr->GetMass();
  Double_t dm1 = mw*mw/s;
  Double_t dp1 = 1.;
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
  Double_t    betaw = phaseWW.GetBetaBar();
  if (betaw <= 0.) return 0.;
  weight *= betaw;

  GENBranch &wmbranch = *wwbranch.GetBranchPtr(0);
  Double_t coswmf  = wmbranch.GetCosTheta();
  Double_t phiwmf  = wmbranch.GetPhi     ();
  Double_t m32     = wmbranch.GetM12();
  Double_t m42     = wmbranch.GetM22();
  GENFrame wwframe = phaseWW.GetFrame();
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
  ANL4DVector qwm = fP[2] + fP[3];
  ANL4DVector qwp = fP[4] + fP[5];
  cerr << " qwm = (" 
       << qwm.E () << ","
       << qwm.Px() << ","
       << qwm.Py() << ","
       << qwm.Pz() << ")" << endl;
  cerr << " qwp = (" 
       << qwp.E () << ","
       << qwp.Px() << ","
       << qwp.Py() << ","
       << qwp.Pz() << ")" << endl;
  ANL4DVector pcm = qwm + qwp + fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
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
  Double_t spin   = 1./2.;                         // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t NNWWBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor() * f5Ptr->GetColor();
  Int_t     ig3   = f3Ptr->GetGenNo() - 1;
  Int_t     ig4   = f4Ptr->GetGenNo() - 1;
  Int_t     ig5   = f5Ptr->GetGenNo() - 1;
  Int_t     ig6   = f6Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig3][ig4]*kVkm[ig5][ig6],2);

  Complex_t amp   = FullAmplitude();
  Double_t  amp2   = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t NNWWBases::FullAmplitude()
{
   Double_t gamw   = fWmBosonPtr->GetWidth();
   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELFermion ne (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing);
   HELFermion neb(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);

   HELFermion f3b(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming); // fubar
   HELFermion f4 (fP[3], fM[3], fHelFinal [3], +1, kIsOutgoing); // fd
   HELVector  wm(f3b, f4, glw, grw, kM_w, gamw);                 // W-

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // fu
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // fdbar
   HELVector  wp(f6b, f5, glw, grw, kM_w, gamw);                 // W+

   Complex_t amp = AmpEEtoNNWW(em, ep, ne, neb, wm, wp);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNWW()
// --------------------------
Complex_t NNWWBases::AmpEEtoNNWW(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &ne,
                                 const HELFermion &neb,
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
   // Higgs Production Amplitude
   //---------------------------
   Double_t  mw     = fWmBosonPtr->GetMass();
   Double_t  gamw   = fWmBosonPtr->GetWidth();
   Double_t  mz     = fZBosonPtr->GetMass();
   Double_t  gamz   = fZBosonPtr->GetWidth();
   Double_t  ma     = 0.;
   Double_t  gama   = 0.;

   Double_t  glwf  = -kGw*kSqh;
   Double_t  grwf  = 0.;

   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  t3n   = +1./2.;
   Double_t  glzn  = -kGz*t3n;
   Double_t  grzn  = 0.;

   Double_t  mne   = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamne = 0.;
   Double_t  game  = 0.;

   //---------------
   // Internal lines
   //---------------
   const HELVector  w1   (em , ne, glwf, grwf, mw , gamw );
   const HELFermion neb07(ep , w1, glwf, grwf, mne, gamne);
   const HELFermion ep09 (neb, w1, glwf, grwf, me , game );

   const HELVector  w2   (neb, ep, glwf, grwf, mw , gamw );
   const HELFermion ne08 (em , w2, glwf, grwf, mne, gamne);
   const HELFermion em10 (ne , w2, glwf, grwf, me , game );

   const HELFermion ne03 (em  , wm, glwf, grwf, mne, gamne);
   const HELVector  z11  (ne03, ne, glzn, grzn, mz , gamz );

   const HELFermion neb04(ep  , wp   , glwf, grwf, mne, gamne);
   const HELVector  z12  (neb , neb04, glzn, grzn, mz , gamz );

   const HELFermion em05 (ne  , wm   , glwf, grwf, me , game );

   const HELFermion ep06 (neb , wp   , glwf, grwf, me , game );

   Double_t  gwwa  = kGe;     
   Double_t  gwwz  = kGw*kCosW;
   Double_t  gw    = kGw;

   Double_t mh     = fMass;
   Double_t gamh   = fHBosonPtr->GetWidth();
   Double_t gwwh   = kGw*mw;

   //-----------
   // Non-fusion
   //-----------
   // ( 1)
   HELVector atee(em, em05, glae, grae, ma, gama);
   Complex_t amp01a = HELVertex(ep06, ep, atee, glae, grae);
   HELVector ztee(em, em05, glze, grze, mz, gamz);
   Complex_t amp01z = HELVertex(ep06, ep, ztee, glze, grze);
   Complex_t amp01 = amp01a + amp01z;

   // ( 2)
   Complex_t amp02 = HELVertex(ep06, ep, z11, glze, grze);

   // ( 3)
   Complex_t amp03 = HELVertex(em, em05, z12, glze, grze);

   // ( 4)
   Complex_t amp04 = HELVertex(neb, neb04, z11, glzn, grzn);

   // ( 5)
   HELVector zvne(ne08, ne, glzn, grzn, mz, gamz);
   Complex_t amp05 = HELVertex(wm, wp, zvne, gwwz);

   // ( 6)
   HELVector w3em(em, em10, glae, grae, glze, grze, mz, gamz);
   Complex_t amp06 = HELVertex(wm, wp, w3em, gw);

   // ( 7)
   HELVector w3ep(ep09, ep, glae, grae, glze, grze, mz, gamz);
   Complex_t amp07 = HELVertex(wm, wp, w3ep, gw);

   // ( 8)
   HELVector zvneb(neb, neb07, glzn, grzn, mz, gamz);
   Complex_t amp08 = HELVertex(wm, wp, zvneb, gwwz);

   // ( 9)
   Complex_t amp09 = HELVertex(ne08, em05, wp, glwf, grwf);

   // (10)
   Complex_t amp10 = HELVertex(ne03, em10, wp, glwf, grwf);

   // (11)
   Complex_t amp11 = HELVertex(ep09, neb04, wm, glwf, grwf);

   // (12)
   Complex_t amp12 = HELVertex(ep06, neb07, wm, glwf, grwf);

   // (13)
   HELVector w3emt(em, em05, glae, grae, glze, grze, mz, gamz);
   Complex_t amp13 = HELVertex(w2, wp, w3emt, gw);

   // (14)
   Complex_t amp14 = HELVertex(w2, wp, z11, gwwz);

   // (15)
   Complex_t amp15 = HELVertex(wm, w1, z12, gwwz);

   // (16)
   HELVector w3ept(ep06, ep, glae, grae, glze, grze, mz, gamz);
   Complex_t amp16 = HELVertex(wm, w1, w3ept, gw);

   //-----------
   // Fusion
   //-----------
   // (17-19)
   Complex_t amp1719 = HELVertex(wm,wp,w2,w1,gwwa,gwwz,mz,gamz);

   //-----------
   // Higgs
   //-----------
#ifndef ANOM_WWH
   // (20)
   HELScalar   ht(wm, w1, gwwh, mh, gamh);
   Complex_t amp20 = HELVertex(w2, wp, ht, gwwh);

   // (21)
   HELScalar   hs(wm, wp, gwwh, mh, gamh);
   Complex_t amp21 = HELVertex(w1, w2, hs, gwwh);
#else
   // (20)
   Double_t g1     = gwwh + 2 * kM_w * kM_w * (fA/fLambda);
   Double_t g2     = -2 * (fB/fLambda);
   Double_t g3     = -4 * (fBtilde/fLambda);
   HELScalar   ht(wm, w1, g1, g2, g3, mh, gamh);
   Complex_t amp20 = HELVertex(w2, wp, ht, g1, g2, g3);

   // (21)
   HELScalar   hs(wm, wp, g1, g2, g3, mh, gamh);
   Complex_t amp21 = HELVertex(w1, w2, hs, g1, g2, g3);
#endif

   //--------------------------
   // Sum up all the amplitudes
   //--------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 + amp05
                 + amp06 + amp07 + amp08 + amp09 + amp10
                 + amp11 + amp12 + amp13 + amp14 + amp15 + amp16
		 + amp1719
		 + amp20 + amp21;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void NNWWBases::Userin()
{
  TDirectory *last = gDirectory;
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("NNWWBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("NNWWBases.BeamstrahlungFilename","trc500"));
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
  //  Initialize W decay table
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

  if (!fHBosonPtr) fHBosonPtr = new HBoson(fMass,fLambda,fA,fB,fBtilde);
  fHBosonPtr->DebugPrint();

  last->cd();
  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs    = fEcmInit;
  Double_t s     = rs*rs;
  Double_t mw    = fWmBosonPtr->GetMass();
  Double_t xilo  = TMath::Log(4.*mw*mw/s);
  Double_t xihi  = 0.;
  Double_t dm1   = mw*mw/s;
  Double_t dp1   = 1.;
  Double_t etlo  = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi  =  TMath::Log((1.+dp1)/dm1)/2.;
  Double_t qwwlo = 0.; //2*mw - 40.;
  Double_t qwwhi = rs;

#if 1
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,  qwwlo,  qwwhi,       50, "qww"   );
  Xh_init( 9,    -1.,    +1.,       50, "cosw-" );
  Xh_init(10,     0.,   k2Pi,       50, "phiw-" );
  Xh_init(11,    60.,   100.,       50, "Mw-"   );
  Xh_init(12,    -1.,    +1.,       50, "CosFu" );
  Xh_init(13,     0.,   k2Pi,       50, "PhiFu" );
  Xh_init(14,    60.,   100.,       50, "Mw+"   );
  Xh_init(15,    -1.,    +1.,       50, "CosFub");
  Xh_init(16,     0.,   k2Pi,       50, "PhiFub");
  Xh_init(17,     0.,     2.,        2, "Helin ");
  Xh_init(18,     0.,     1.,        1, "Helot ");
  Xh_init(19,     0.,    12.,       12, "W- mode");
  Xh_init(20,     0.,    12.,       12, "W+ mode");
#endif
#ifdef TEMP_H
  H1Init("hMww"     ,"", 200,       0.,     fEcmInit);
  H1Init("hCosWm"   ,"", 100,      -1.,          +1.);
  H1Init("hRSH"     ,"",1100,       0., fEcmInit*1.1);
#endif
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void NNWWBases::Userout()
{
  cout << "End of NNWWBases----------------------------------- "  << endl
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
void NNWWBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][6] = {{-1, +1, +1, -1, -1, +1}};
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
   fHelFinal  [4] = kFHelComb[fJCombF][4];
   fHelFinal  [5] = kFHelComb[fJCombF][5];
   weight = kNf;
}
