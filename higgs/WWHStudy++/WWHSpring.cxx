//*****************************************************************************
//* =====================
//*  WWHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> WWH generator
//*
//* (Update Record)
//*    2010/04/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "WWHSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(WWHSpring)
ClassImp(WWHSpringBuf)
ClassImp(WWHBases)

//-----------------------------------------------------------------------------
// ==============================
//  class WWHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WWHSpring::WWHSpring(const char *name,
                     const char *title,
                     WWHBases   *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new WWHSpringBuf("WWHSpringBuf",
                               "WWHSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new WWHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
WWHSpring::~WWHSpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t WWHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    WWHBases *bs = static_cast<WWHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> WWHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class WWHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t WWHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  WWHBases     *bases   = static_cast<WWHBases *>(
		          static_cast<WWHSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 7;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[1] = bases->fP[1];
  pv[2] = bases->fP[2];
  pv[3] = bases->fP[3];
  pv[4] = bases->fP[4];
  pv[5] = pv[1] + pv[2];
  pv[6] = pv[3] + pv[4];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2WW         = bases->GetQ2WW();
  fCosThetaW    = bases->GetCosThetaW();
  fPhiW         = bases->GetPhiW();
  fQ2Wm         = bases->GetQ2Wm();
  fCosThetaWmF  = bases->GetCosThetaWmF();
  fPhiWmF       = bases->GetPhiWmF();
  fQ2Wp         = bases->GetQ2Wp();
  fCosThetaWpF  = bases->GetCosThetaWpF();
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
  Int_t    idh     = 25;                          // PDG code for H
  Double_t mass    = bases->GetMass();

  Int_t    idw     = 24;                          // PDG code for W
  // W-
  Int_t    idf1    = bases->f1Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg1   = bases->f1Ptr->GetCharge();   // f1 charge
  Double_t m1      = bases->f1Ptr->GetMass  ();   // f1 mass
  Int_t    hel1    = bases->fHelFinal[1];         // f1 helicity
  Double_t color1  = bases->f1Ptr->GetColor();    // color factor for f1

  Int_t    idf2    = bases->f2Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg2   = bases->f2Ptr->GetCharge();   // f2 charge
  Double_t m2      = bases->f2Ptr->GetMass  ();   // f2 mass
  Int_t    hel2    = bases->fHelFinal[2];         // f2 helicity

  Int_t    islevwm = color1 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfwm   = 2;                           // color flux id
  Double_t rq2wm   = pv[5].Mag();

  // W+
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f1 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f1 mass
  Int_t    hel3    = bases->fHelFinal[3];         // f1 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f1

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f2 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f2 mass
  Int_t    hel4    = bases->fHelFinal[4];         // f2 helicity

  Int_t    islevwp = color3 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfwp   = 3;                           // color flux id
  Double_t rq2wp   = pv[6].Mag();

#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idh << " m=" << mass  << " Q=" << 0.
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 4 pid=" << idf1 << " m=" << m1 << " Q=" << chrg1 
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 5 pid=" << idf2 << " m=" << m2 << " Q=" << chrg2 
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 6 pid=" << idf3 << " m=" << m3 << " Q=" << chrg3 
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  cerr << " 7 pid=" << idf4 << " m=" << m4 << " Q=" << chrg4
       << " pv=(" << (*qp[4])(0) << ", "
                  << (*qp[4])(1) << ", "
                  << (*qp[4])(2) << ", "
                  << (*qp[4])(3) << ") " << endl;
  TVector qcm(4), qh(4), qwm(4), qwp(4);
  qh  = *qp[0];
  qwm = *qp[1] + *qp[2];
  qwp = *qp[3] + *qp[4];
  qcm = qh + qwm + qwp;

  cerr << " ph=(" << qh[0] << ", "
                  << qh[1] << ", "
                  << qh[2] << ", "
                  << qh[3] << ") " << endl;
  cerr << " pwm=(" << qwm[0] << ", "
                   << qwm[1] << ", "
                   << qwm[2] << ", "
                   << qwm[3] << ") " << endl;
  cerr << " pwp=(" << qwp[0] << ", "
                   << qwp[1] << ", "
                   << qwp[2] << ", "
                   << qwp[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                                No.  PID  Mass  Charge   pv   Nd 1st Mom  hel  col  shower
  new (partons[0]) JSFSpringParton( 1,  idh, mass,    0., *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2, -idw,rq2wm,   -1., *qp[5], 2, 4,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3,  idw,rq2wp,   +1., *qp[6], 2, 6,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton( 4,-idf1,   m1,-chrg1, *qp[1], 0, 0,  2, hel1,icfwm,islevwm);
  new (partons[4]) JSFSpringParton( 5,-idf2,   m2,-chrg2, *qp[2], 0, 0,  2, hel2,icfwm,islevwm);
  new (partons[5]) JSFSpringParton( 6, idf3,   m3, chrg3, *qp[3], 0, 0,  3, hel3,icfwp,islevwp);
  new (partons[6]) JSFSpringParton( 7, idf4,   m4, chrg4, *qp[4], 0, 0,  3, hel4,icfwp,islevwp);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class WWHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WWHBases::WWHBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fMass       ( 120.),
           fEcmInit    (1000.),
           fISR        ( 1),
           fBeamStr    ( 1),
           fBeamWidth  (0.002),
           fPole       (0.),
           fWmModesLo  ( 1),
           fWmModesHi  (12),
           fWpModesLo  ( 1),
           fWpModesHi  (12),
	   fNCALL      (80000),
	   fACC1       (0.05),
	   fACC2       (0.05),
	   fITMX1      (20),
	   fITMX2      (40),
           fZBosonPtr  ( 0),
           fWmBosonPtr ( 0),
           fWpBosonPtr ( 0),
           fZBoost     (0.),
           fEcmIP      (fEcmInit),
           fQ2WW       (0.),
           fQ2Wm       (0.),
           fQ2Wp       (0.),
           fWmModePtr  (0),
           f1Ptr       (0),
           f2Ptr       (0),
           fWpModePtr  (0),
           f3Ptr       (0),
           f4Ptr       (0),
           fCosTheta   (0.),
           fPhi        (0.),
           fXQ2WW      (0.),
           fCosThetaW  (0.),
           fPhiW       (0.),
           fXQ2Wm      (0.),
           fCosThetaWmF(0.),
           fPhiWmF     (0.),
           fXQ2Wp      (0.),
           fCosThetaWpF(0.),
           fPhiWpF     (0.),
           fR_BW_m     (0),
           fR_BW_p     (0),
           fR_BS_m     (0),
           fR_BS_p     (0),
           fR_ISR_var  (0),
           fR_ISR_side (0)
{
  //  Constructor of bases.  Default parameter should be initialized here
  //
  // --------------------------------------------
  //  Get parameters from jsf.conf, if specified
  // --------------------------------------------

  cout << "Init zhbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("WWHBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.CosthWRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.PhiWOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.CosthWmFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.PhiWmFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.CosthWpFRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.PhiWpFOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.ACC1","0.05"));
  ins >> fACC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.ACC2","0.05"));
  ins >> fACC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.ITMX1","20"));
  ins >> fITMX1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.ITMX2","40"));
  ins >> fITMX2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWHBases.NCALL","80000"));
  ins >> fNCALL;

  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fWmDecayMode   , 0., 1., 0, 1);
  DefineVariable(fWpDecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2WW         , 0., 1., 1, 1);
  DefineVariable(fXQ2Wm         , 0., 1., 0, 1);
  DefineVariable(fXQ2Wp         , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] *= TMath::Pi();
  fXU[1] *= TMath::Pi();
  fXL[3] *= TMath::Pi();
  fXU[3] *= TMath::Pi();
  fXL[5] *= TMath::Pi();
  fXU[5] *= TMath::Pi();
  fXL[7] *= TMath::Pi();
  fXU[7] *= TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaW  , fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiW       , fXL[3], fXU[3], 1, 1);
  DefineVariable(fCosThetaWmF, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiWmF     , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaWpF, fXL[6], fXU[6], 1, 1);
  DefineVariable(fPhiWpF     , fXL[7], fXU[7], 0, 1);

  //--
  //  Beamstrahlung and natural Ebeam spread
  //--
  if (fBeamStr == 1) {
    DefineVariable(fR_BW_m, 0., 1., 0, 1);
    DefineVariable(fR_BW_p, 0., 1., 0, 1);
    DefineVariable(fR_BS_m, 0., 1., 1, 1);
    DefineVariable(fR_BS_p, 0., 1., 1, 1);
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
WWHBases::~WWHBases()
{
  delete fZBosonPtr;
  delete fWmBosonPtr;
  delete fWpBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t WWHBases::Func()
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
  f1Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fWpModePtr = fWpBosonPtr->PickMode(fWpDecayMode, weight, fWpMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t mh   = fMass;
  if (fEcmIP < mh + m1 + m2 + m3 + m4) {
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
  Double_t s      = fEcmIP*fEcmIP;
  Double_t rs     = fEcmIP;
  Double_t qwmmin = m1 + m2;
  Double_t qwmmax = rs - (m3 + m4 + mh);
#ifndef __ZEROWIDTH__
  fQ2Wm = fWmBosonPtr->GetQ2BW(qwmmin, qwmmax, fXQ2Wm, weight);
#else
  fQ2Wm = TMath::Power(fWmBosonPtr->GetMass(),2);
  weight = kPi*fWmBosonPtr->GetMass()*fWmBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2wm  = TMath::Sqrt(fQ2Wm);
  Double_t qwpmin = m3 + m4;
  Double_t qwpmax = rs - (rq2wm + mh);
#ifndef __ZEROWIDTH__
  fQ2Wp = fWpBosonPtr->GetQ2BW(qwpmin, qwpmax, fXQ2Wp, weight);
#else
  fQ2Wp = TMath::Power(fWpBosonPtr->GetMass(),2);
  weight = kPi*fWpBosonPtr->GetMass()*fWpBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // WW system
  Double_t rq2wp   = TMath::Sqrt(fQ2Wp);
  Double_t qww2min = TMath::Power(rq2wm + rq2wp,2);
  Double_t qww2max = TMath::Power(rs - mh,2);
           fQ2WW   = qww2min + (qww2max - qww2min)*fXQ2WW;
  bsWeight *= qww2max - qww2min;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch wmbranch(fQ2Wm, fCosThetaWmF, fPhiWmF, m1*m1    , m2*m2);
  GENBranch wpbranch(fQ2Wp, fCosThetaWpF, fPhiWpF, m3*m3    , m4*m4);
  GENBranch wwbranch(fQ2WW, fCosThetaW  , fPhiW  , &wmbranch, &wpbranch);
  GENBranch cmbranch(s    , fCosTheta   , fPhi   , mh*mh    , &wwbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2WW), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaW        , (bsWeight*sigma));
  Xh_fill( 6, fPhiW             , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Wm), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaWmF      , (bsWeight*sigma));
  Xh_fill( 9, fPhiWmF           , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2Wp), (bsWeight*sigma));
  Xh_fill(11, fCosThetaWpF      , (bsWeight*sigma));
  Xh_fill(12, fPhiWpF           , (bsWeight*sigma));
  Xh_fill(13, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(14, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(15, (Double_t)fWmMode , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fWpMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t WWHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2h  = cmbranch.GetM12();
  Double_t q2ww = cmbranch.GetM22();
  Double_t cosh = cmbranch.GetCosTheta();
  Double_t phih = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2h, q2ww, cmframe, cosh, phih, 0);
  ANL4DVector ph  = phaseCM.GetFourMomentum(0);
  ANL4DVector pww = phaseCM.GetFourMomentum(1);
  Double_t betah  = phaseCM.GetBetaBar();
  if (betah <= 0.) return 0.;
  fP[0] = ph;
  fM[0] = fMass;

  GENBranch &wwbranch = * cmbranch.GetBranchPtr(1);
  Double_t cosw = wwbranch.GetCosTheta();
  Double_t phiw = wwbranch.GetPhi     ();
  Double_t qwm2 = wwbranch.GetM12();
  Double_t qwp2 = wwbranch.GetM22();
  GENPhase2 phaseWW(pww, qwm2, qwp2, cmframe, cosw, phiw, 1);
  ANL4DVector pwm = phaseWW.GetFourMomentum(0);
  ANL4DVector pwp = phaseWW.GetFourMomentum(1);
  Double_t betaw  = phaseWW.GetBetaBar();
  if (betaw <= 0.) return 0.;

  GENBranch &wmbranch = * wwbranch.GetBranchPtr(0);
  Double_t coswmf = wmbranch.GetCosTheta();
  Double_t phiwmf = wmbranch.GetPhi     ();
  Double_t m12    = wmbranch.GetM12();
  Double_t m22    = wmbranch.GetM22();
  GENPhase2 phaseWm(pwm, m12, m22, cmframe, coswmf, phiwmf, 1);
  fP[1] = phaseWm.GetFourMomentum(0);
  fP[2] = phaseWm.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m12);
  fM[2] = TMath::Sqrt(m22);
  Double_t betawmf = phaseWm.GetBetaBar();
  if (betawmf <= 0.) return 0.;

  GENBranch &wpbranch = * wwbranch.GetBranchPtr(1);
  Double_t coswpf = wpbranch.GetCosTheta();
  Double_t phiwpf = wpbranch.GetPhi     ();
  Double_t m32    = wpbranch.GetM12();
  Double_t m42    = wpbranch.GetM22();
  GENPhase2 phaseWp(pwp, m32, m42, cmframe, coswpf, phiwpf, 1);
  fP[3] = phaseWp.GetFourMomentum(0);
  fP[4] = phaseWp.GetFourMomentum(1);
  fM[3] = TMath::Sqrt(m32);
  fM[4] = TMath::Sqrt(m42);
  Double_t betawpf = phaseWp.GetBetaBar();
  if (betawpf <= 0.) return 0.;

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
  for (Int_t i=0; i<5; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qwm = fP[1] + fP[2];
  ANL4DVector qwp = fP[3] + fP[4];
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
  ANL4DVector pcm = fP[0] + qwm + qwp;
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
  static const Int_t    kNbr  = 4;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                                        // identical particle factor
  Double_t dPhase = kFact * betah * betaw * betawmf * betawpf; // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);                       // beam flux factor
  Double_t spin   = 1./2.;                                     // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase;      // in [1/GeV^2]
           sigma *= kGeV2fb;                                   // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t WWHBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor();
  Int_t     ig1   = f1Ptr->GetGenNo() - 1;
  Int_t     ig2   = f2Ptr->GetGenNo() - 1;
  Int_t     ig3   = f3Ptr->GetGenNo() - 1;
  Int_t     ig4   = f4Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig1][ig2]*kVkm[ig3][ig4],2);

  Complex_t amp    = FullAmplitude();
  Double_t  amp2   = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t WWHBases::FullAmplitude()
{
   Double_t gamw   = fWmBosonPtr->GetWidth();
   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

   HELScalar  hs(fP[0]); // higgs

   HELFermion f1b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // fubar
   HELFermion f2 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fd
   HELVector  wm(f1b, f2, glw, grw, kM_w, gamw);                 // W-

   HELFermion f3 (fP[3], fM[3], fHelFinal [3], +1, kIsOutgoing); // fu
   HELFermion f4b(fP[4], fM[4], fHelFinal [4], -1, kIsIncoming); // fdbar
   HELVector  wp(f4b, f3, glw, grw, kM_w, gamw);                 // W+

   Complex_t amp = AmpEEtoWWH(em, ep, wm, wp, hs);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoWWH()
// --------------------------
Complex_t WWHBases::AmpEEtoWWH(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &wm,
                               const HELVector  &wp,
                               const HELScalar  &hs)
{
   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  gamw  = fWmBosonPtr->GetWidth();

   Double_t  gwwz  = kGw*kCosW;
   Double_t  gww3  = kGw;
   Double_t  glw   = -kGw*kSqh;
   Double_t  grw   = 0.;

   Double_t  gwwh  = kGw*kM_w;
   Double_t  gzzh  = kGz*kM_z;

   Double_t  mne   = kMass[0][0][0];
   Double_t  gamne = 0.;

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
   // WWH Production Amplitude
   //---------------------------
   //--
   // S-channel
   //--
   HELVector w3(em, ep, glae, grae, glze, grze, kM_z, gamz); // s-channel W3
   HELVector wmstar(wm, hs, gwwh, kM_w, gamw);
   HELVertex ampwsw3(wmstar, wp, w3, gww3);
   HELVector wpstar(wp, hs, gwwh, kM_w, gamw);
   HELVertex ampwws3(wm, wpstar, w3, gww3);

   //--
   // T-channel ne
   //--
   HELFermion nes(em, wmstar, glw, grw, mne, gamne);
   HELVertex  amptwsw(nes,ep,wp,glw,grw);
   HELFermion ne(em, wm, glw, grw, mne, gamne);
   HELVertex  amptwws(ne,ep,wpstar,glw,grw);

   //--
   // Higgs from Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel Z
   HELVector zf(wm, wp, gwwz, kM_z, gamz);
   Complex_t ampzh = HELVertex(zs, zf, hs, gzzh);

   Complex_t amp = ampwsw3 + ampwws3 
		 + amptwsw + amptwws
		 + ampzh;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void WWHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("WWHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("WWHBases.BeamstrahlungFilename","trc500"));
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
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4,     0., fEcmInit    , 50, "Qww"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosW"   );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiW"   );
  Xh_init( 7,    60.,   100.,       50, "Mw-"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosFd"  );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiFd"  );
  Xh_init(10,    60.,   100.,       50, "Mw+"    );
  Xh_init(11, fXL[6], fXU[6],       50, "CosFdb" );
  Xh_init(12, fXL[7], fXU[7],       50, "PhiFdb" );
  Xh_init(13,     0.,     2.,        2, "Helin " );
  Xh_init(14,     0.,     1.,        2, "Helot " );
  Xh_init(15,     0.,    12.,       12, "W- mode");
  Xh_init(16,     0.,    12.,       12, "W+ mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void WWHBases::Userout()
{
  cout << "End of WWHBases----------------------------------- "  << endl
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
void WWHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][5] = {{0, +1, -1, -1, +1}};
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
   weight = kNf;
}
