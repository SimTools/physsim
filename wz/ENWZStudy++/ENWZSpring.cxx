//*****************************************************************************
//* =====================
//*  ENWZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ENWZ generator
//*
//* (Update Record)
//*    2014/10/13  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ENWZSpring.h"

#include <sstream>
#include <iomanip>
//#define __ZEROWIDTH__
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ENWZSpring)
ClassImp(ENWZSpringBuf)
ClassImp(ENWZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ENWZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ENWZSpring::ENWZSpring(const char      *name,
                       const char      *title,
                             ENWZBases *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new ENWZSpringBuf("ENWZSpringBuf",
                                "ENWZSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new ENWZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ENWZSpring::~ENWZSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ENWZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ENWZBases *bs = static_cast<ENWZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ENWZBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ENWZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ENWZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ENWZBases    *bases   = static_cast<ENWZBases *>(
                          static_cast<ENWZSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 8;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // e-
  pv[1] = bases->fP[1]; // neb
  pv[2] = bases->fP[2]; // fu
  pv[3] = bases->fP[3]; // fdb
  pv[4] = bases->fP[4]; // f
  pv[5] = bases->fP[5]; // fb
  pv[6] = pv[2] + pv[3];// W+
  pv[7] = pv[4] + pv[5];// Z

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------
  fCP           = bases->GetCP();
  for (Int_t i=0; i<fNparton; i++) pv[i] *= ((Double_t)fCP);
  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fXi           = bases->GetXi();
  fEta1         = bases->GetEta1();
  fEta2         = bases->GetEta2();
  fPhi1         = bases->GetPhi1();
  fPhi2         = bases->GetPhi2();
  fQ2WZ         = bases->GetQ2WZ();
  fCosW         = bases->GetCosW ();
  fPhiW         = bases->GetPhiW ();
  fQ2W          = bases->GetQ2W  ();
  fCosWF        = bases->GetCosWF();
  fPhiWF        = bases->GetPhiWF();
  fQ2Z          = bases->GetQ2Z  ();
  fCosZF        = bases->GetCosZF();
  fPhiZF        = bases->GetPhiZF();
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
  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    idw     = fCP*24;                      // PDG code for W
  Int_t    ide     = fCP*11;                      // PDG code for e
  Int_t    idn     = fCP*(-12);                   // PDG code for e-
  Double_t mn      = 0.;                          // electron mass
  Double_t me      = kM_e;                        // electron mass
  Double_t chrgw   = fCP;                         // W charge
  Double_t chrge   = -fCP;                        // electron charge

  // W
  Int_t    idf3    = fCP*bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = fCP*bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      =     bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = fCP*bases->fHelFinal[2];         // f3 helicity
  Int_t    idf4    = fCP*bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = fCP*bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      =     bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = fCP*bases->fHelFinal[3];         // f4 helicity
  Double_t color3  =  bases->f3Ptr->GetColor();    // color factor for f3
  Int_t    ilsevw = color3 > 1. ? 101 : 0; 	   // shower level
  Int_t    icfw   = 2;                             // color flux id
  Double_t rq2w   = pv[6].Mag();

  // Z
  Int_t    idf5    =     bases->f5Ptr->GetPID   ();   // PDG code for f5
  Double_t chrg5   =     bases->f5Ptr->GetCharge();   // f5 charge
  Double_t m5      =     bases->f5Ptr->GetMass  ();   // f5 mass
  Int_t    hel5    = fCP*bases->fHelFinal[4];         // f5 helicity
  Int_t    idf6    =     bases->f6Ptr->GetPID   ();   // PDG code for f6
  Double_t chrg6   =     bases->f6Ptr->GetCharge();   // f6 charge
  Double_t m6      =     bases->f6Ptr->GetMass  ();   // f6 mass
  Int_t    hel6    = fCP*bases->fHelFinal[5];         // f6 helicity
  Double_t color5  =     bases->f5Ptr->GetColor(); // color factor for f5
  Int_t    islevz = color5 > 1. ? 201 : 0;  	   // shower level
  Int_t    icfz   = 3;                             // color flux id
  Double_t rq2z   = pv[7].Mag();
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
  new (partons[0]) JSFSpringParton(1,  ide,   me, chrge, *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton(2,  idn,   mn,    0., *qp[1], 0, 0,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton(3,  idw, rq2w, chrgw, *qp[6], 2, 5,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton(4,  idz, rq2z,    0., *qp[7], 2, 7,  0,    0,    0,      0);
  new (partons[4]) JSFSpringParton(5, idf3,   m3, chrg3, *qp[2], 0, 0,  3, hel3, icfw, ilsevw);
  new (partons[5]) JSFSpringParton(6, idf4,   m4, chrg4, *qp[3], 0, 0,  3, hel4, icfw, ilsevw);
  new (partons[6]) JSFSpringParton(7, idf5,   m5, chrg5, *qp[4], 0, 0,  4, hel5, icfz, islevz);
  new (partons[7]) JSFSpringParton(8, idf6,   m6, chrg6, *qp[5], 0, 0,  4, hel6, icfz, islevz);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ENWZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ENWZBases::ENWZBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 125.),
           fWidth     ( 0.004),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPolem     (0.),
           fPolep     (0.),
           fWModesLo  ( 1),
           fWModesHi  (12),
           fZModesLo  ( 1),
           fZModesHi  (12),
	   fACC1      (0.05),
	   fACC2      (0.05),
	   fITMX1     (20),
	   fITMX2     (40),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2WZ      (0.),
	   fCosW      (0.),
	   fPhiW      (0.),
	   fQ2W       (0.),
	   fCosWF     (0.),
	   fPhiWF     (0.),
	   fQ2Z       (0.),
	   fCosZF     (0.),
	   fPhiZF     (0.),
           fWModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fZModePtr  (0),
           f5Ptr      (0),
           f6Ptr      (0),
           fCP        ( 1),
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

  cout << "Init enwzbases " << endl;
  
  using namespace std;

  stringstream ins(gJSF->Env()->GetValue("ENWZBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.WidthH","0.006")); // M_x [GeV]
  ins >> fWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.Polem","0."));       // electron polarization
  ins >> fPolem;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.Polep","0."));       // positron polarization
  ins >> fPolep;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.CosthWRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.PhiWOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.CosthWFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.PhiWFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.CosthZFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.PhiZFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.WModesLo","1"));      // W decay mode lo
  ins >> fWModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.WModesHi","12"));     // W decay mode hi
  ins >> fWModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWZBases.NCALL","80000"));
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
  if (fISR == 1) {
    DefineVariable(fR_ISR_var , 0., 1., 1, 1);
    DefineVariable(fR_ISR_side, 0., 1., 0, 0);
  }

  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fXComb       , 0., 1., 1, 1);

  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1, q2wz
  //--
  DefineVariable(fXXi    , 0., 1., 1, 1);
  DefineVariable(fXEta1  , 0., 1., 1, 1);
  DefineVariable(fXEta2  , 0., 1., 1, 1);
  DefineVariable(fXQ2WZ  , 0., 1., 1, 1);
  //--
  //  Q2 of W and Z
  //--
  DefineVariable(fXQ2W   , 0., 1., 0, 1);
  DefineVariable(fXQ2Z   , 0., 1., 0, 1);

  DefineVariable(fXPhi21 , 0., 1., 0, 1);
  DefineVariable(fXPhi1  , 0., 1., 0, 1);

  //--
  //  cos(theta) and phi
  //--
  fXL[1] *= TMath::Pi();
  fXU[1] *= TMath::Pi();
  fXL[3] *= TMath::Pi();
  fXU[3] *= TMath::Pi();
  fXL[5] *= TMath::Pi();
  fXU[5] *= TMath::Pi();

  DefineVariable(fCosW , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhiW , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosWF, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiWF, fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosZF, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiZF, fXL[5], fXU[5], 0, 1);

  //--
  //  Final states
  //--
  DefineVariable(fWDecayMode   , 0., 1., 0, 1);
  DefineVariable(fZDecayMode   , 0., 1., 0, 1);

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
ENWZBases::~ENWZBases()
{
  delete fWBosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ENWZBases::Func()
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

  GENDecayMode *fWModePtr = fWBosonPtr->PickMode(fWDecayMode, weight, fWMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fWModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fWModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f5Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f6Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m5   = f5Ptr->GetMass();
  Double_t m6   = f6Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = kM_e; // Me
  Double_t m2   = 0.;   // Mneb
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

  // W
#ifndef __ZEROWIDTH__
  Double_t qwmin = m3 + m4;
  Double_t qwmax = rs - (m1 + m2 + m5 + m6);
  fQ2W = fWBosonPtr->GetQ2BW(qwmin, qwmax, fXQ2W, weight);
#else
  fQ2W = TMath::Power(fWBosonPtr->GetMass(),2);
  weight = kPi*fWBosonPtr->GetMass()*fWBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2w  = TMath::Sqrt(fQ2W);

  // Z
#ifndef __ZEROWIDTH__
  Double_t qzmin = m5 + m6;
  Double_t qzmax = rs - (m1 + m2 + rq2w);
  fQ2Z = fZBosonPtr->GetQ2BW(qzmin, qzmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2z  = TMath::Sqrt(fQ2Z);

  // Q2 of WZ
  Double_t qwz2mn = TMath::Power(rq2w + rq2z, 2);
  Double_t qwz2mx = TMath::Power(rs - (m1+m2),2);
  fQ2WZ     = qwz2mn + (qwz2mx-qwz2mn)*fXQ2WZ;
  bsWeight *= qwz2mx - qwz2mn;
  Double_t rq2wz  = TMath::Sqrt(fQ2WZ);


  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t xilo  = TMath::Log(rq2wz*(rq2wz+m1+m2)/s);
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

  GENBranch wbranch(fQ2W, fCosWF, fPhiWF, m3*m3    , m4*m4);
  GENBranch zbranch(fQ2Z, fCosZF, fPhiZF, m5*m5    , m6*m6);
  GENBranch wzbranch(fQ2WZ, fCosW , fPhiW , &wbranch, &zbranch);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX(wzbranch) * qed;

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
  Xh_fill( 8, TMath::Sqrt(fQ2WZ), (bsWeight*sigma));
  Xh_fill( 9, fCosW             , (bsWeight*sigma));
  Xh_fill(10, fPhiW             , (bsWeight*sigma));
  Xh_fill(11, TMath::Sqrt(fQ2W) , (bsWeight*sigma));
  Xh_fill(12, fCosWF            , (bsWeight*sigma));
  Xh_fill(13, fPhiWF            , (bsWeight*sigma));
  Xh_fill(14, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill(15, fCosZF            , (bsWeight*sigma));
  Xh_fill(16, fPhiZF            , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fWMode  , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fZMode  , (bsWeight*sigma));
  Xh_fill(21, (Double_t)fCP     , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ENWZBases::DSigmaDX(GENBranch &wzbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs  = fEcmIP;
  Double_t s   = rs*rs; 
  Double_t eb  = rs/2.;

  Double_t qwz2 = wzbranch.GetQ2();  // q^2_{WZ}
  Double_t qwz  = TMath::Sqrt(qwz2);
  Double_t q2w  = wzbranch.GetM12(); // W
  Double_t q2z  = wzbranch.GetM22(); // Z

  Double_t rxi = TMath::Exp(fXi);
  weight *= rxi;

  Double_t m1  = fM[0]; // Me
  Double_t m2  = fM[1]; // Mneb
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
    Double_t rximn = (qwz-m1+m2)*(qwz+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - qwz2/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - qwz2/s);
  } else {
    Double_t rximn = (qwz+m1-m2)*(qwz+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - qwz2/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - qwz2/s);
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
  ANL4DVector qvwz  = qcm - fP[0] - fP[1];
  GENFrame    cmframe;
  GENPhase2   phaseWZ(qvwz, q2w, q2z, cmframe, fCosW, fPhiW, 1);
  ANL4DVector pw    = phaseWZ.GetFourMomentum(0);
  ANL4DVector pz    = phaseWZ.GetFourMomentum(1);
  Double_t    betaz = phaseWZ.GetBetaBar();
  if (betaz <= 0.) return 0.;
  weight *= betaz;

  GENBranch &wbranch = *wzbranch.GetBranchPtr(0);
  Double_t coswf   = wbranch.GetCosTheta();
  Double_t phiwf   = wbranch.GetPhi     ();
  Double_t m32     = wbranch.GetM12();
  Double_t m42     = wbranch.GetM22();
  GENFrame wzframe = phaseWZ.GetFrame(1);
  GENPhase2 phaseW(pw, m32, m42, wzframe, coswf, phiwf, 1);
  fP[2] = phaseW.GetFourMomentum(0);
  fP[3] = phaseW.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betawf = phaseW.GetBetaBar();
  if (betawf <= 0.) return 0.;
  weight *= betawf;

  GENBranch &zbranch = *wzbranch.GetBranchPtr(1);
  Double_t coszf  = zbranch.GetCosTheta();
  Double_t phizf  = zbranch.GetPhi     ();
  Double_t m52    = zbranch.GetM12();
  Double_t m62    = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m52, m62, wzframe, coszf, phizf, 1);
  fP[4] = phaseZ.GetFourMomentum(0);
  fP[5] = phaseZ.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(m52);
  fM[5] = TMath::Sqrt(m62);
  Double_t betazf = phaseZ.GetBetaBar();
  if (betazf <= 0.) return 0.;
  weight *= betazf;

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
  ANL4DVector qvw = fP[2] + fP[3];
  cerr << " pw = (" 
       << qvw.E () << ","
       << qvw.Px() << ","
       << qvw.Py() << ","
       << qvw.Pz() << ")" << endl;
  ANL4DVector qvz = fP[4] + fP[5];
  cerr << " pz = (" 
       << qvz.E () << ","
       << qvz.Px() << ","
       << qvz.Py() << ","
       << qvz.Pz() << ")" << endl;
  ANL4DVector pcm = qvw + qvz + fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
  cerr << " wmass = " << qvw.GetMass() << endl;
  cerr << " zmass = " << qvz.GetMass() << endl;
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
Double_t ENWZBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor() * f5Ptr->GetColor();
  Int_t     ig3   = f3Ptr->GetGenNo() - 1;
  Int_t     ig4   = f4Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig3][ig4],2);

  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ENWZBases::FullAmplitude()
{
   Double_t mw     = fWBosonPtr->GetMass();
   Double_t gamw   = fWBosonPtr->GetWidth();
   Double_t mz     = fZBosonPtr->GetMass();
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   Double_t qf5    = f5Ptr->GetCharge();
   Double_t t3f5   = f5Ptr->GetISpin();
   Double_t glzf5  = -kGz*(t3f5 - qf5*kSin2W);
   Double_t grzf5  = -kGz*(     - qf5*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em (fK[0], kM_e,  fHelInitial[0], +1, kIsIncoming);
   HELFermion ep (fK[1], kM_e,  fHelInitial[1], -1, kIsOutgoing);

   HELFermion emf(fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // final e-
   HELFermion nbf(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // final neb

   HELFermion fu (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fu
   HELFermion fdb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // fdb
   HELVector  wf(fdb, fu, glw, grw, mw, gamw);                   // W+

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // f
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // fb
   HELVector  zf(f6b, f5, glzf5, grzf5, mz, gamz);               // Z

   Complex_t amp = AmpEEtoENWZ(em, ep, emf, nbf, wf, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoENWZ()
// --------------------------
Complex_t ENWZBases::AmpEEtoENWZ(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &emf,
                                 const HELFermion &nbf,
                                 const HELVector  &wf,
                                 const HELVector  &zf)
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
   // WZ Production Amplitudes
   //---------------------------
   Double_t  mn    = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamn  = 0.;
   Double_t  game  = 0.;

   Double_t  mw    = fWBosonPtr->GetMass();
   Double_t  gamw  = fWBosonPtr->GetWidth();
   Double_t  mz    = fZBosonPtr->GetMass();
   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  ma    = 0.;
   Double_t  gama  = 0.;
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

   //-------------------------
   // Calculate internal lines
   //-------------------------
   // virtual almost real A/Z
   // (e- leg)
   HELVector  wrk01(ebm  , e1, fSh1, fCh1, fPhi1, 
                    fHelInitial[0], fHelFinal[0],+1, kGe, me);
   HELVector  wrk02(em   , emf  , glze, grze, mz, gamz);

   HELVector  wrk03(nbf  , ep   , glw , grw , mw, gamw);

   HELFermion wrk04(emf  , zf   , glze, grze, me, game);
   HELVector  wrk05(em   , wrk04, glae, grae, ma, gama);
   HELVector  wrk06(em   , wrk04, glze, grze, mz, gamz);

   HELFermion wrk07(nbf  , zf   , glzn, grzn, mn, gamn);
   HELFermion wrk08(nbf  , wf   , glw , grw , me, game);

   HELFermion wrk09(emf  , wf   , glw , grw , mn, gamn);
   HELVector  wrk10(em   , wrk09, glw , grw , mw, gamw);

   HELFermion wrk11(ep   , zf   , glze, grze, me, game);
   HELVector  wrk12(nbf  , wrk11, glw , grw , mw, gamw);

   HELFermion wrk13(wrk07, wf   , glw , grw , me, game);
   HELFermion wrk14(wrk08, zf   , glze, grze, me, game);

   HELVector  wrk15(wf   , zf   , gwwz,       mw, gamw);
   HELFermion wrk16(nbf  , wrk15, glw , grw , me, game);

   HELVector  wrk17(wrk01[0]*(gwwa/gw)+wrk02[0]*(gwwz/gw),
                    wrk01[1]*(gwwa/gw)+wrk02[1]*(gwwz/gw),
                    wrk01[2]*(gwwa/gw)+wrk02[2]*(gwwz/gw),
                    wrk01[3]*(gwwa/gw)+wrk02[3]*(gwwz/gw),
                    wrk01.GetFourMomentum());

   HELFermion wrk18(ep   , wf   , glw , grw , mn, gamn);

   HELFermion wrk19(em   , zf   , glze, grze, me, game);

   HELFermion wrk20(em   , wrk03, glw , grw , mn, gamn);

   HELFermion wrk21(nbf  , wrk02, glzn, grzn, mn, gamn);

   HELVector  wrk22(wrk19, emf  , glae, grae, ma, gama);
   HELVector  wrk23(wrk19, emf  , glze, grze, mz, gamz);

   //---------------------------
   // Now calculate amplitudes
   //---------------------------
   // Non-fusion diagrams
   // ( 1)
   HELVertex amp01a(wrk08, wrk11, wrk01, glae, grae);
   HELVertex amp01b(wrk08, wrk11, wrk02, glze, grze);
   Complex_t amp01 = amp01a + amp01b;
   // ( 2)
   HELVertex amp02a(wrk13, ep   , wrk01, glae, grae);
   HELVertex amp02b(wrk13, ep   , wrk02, glze, grze);
   Complex_t amp02 = amp02a + amp02b;
   // ( 3)
   HELVertex amp03a(wrk14, ep   , wrk01, glae, grae);
   HELVertex amp03b(wrk14, ep   , wrk02, glze, grze);
   Complex_t amp03 = amp03a + amp03b;
   // ( 4)
   HELVertex amp04a(wrk16, ep   , wrk01, glae, grae);
   HELVertex amp04b(wrk16, ep   , wrk02, glze, grze);
   Complex_t amp04 = amp04a + amp04b;
   // ( 5)
   HELVector tmp05 (wrk07, ep   , glw, grw, mw, gamw);
   HELVertex amp05 (tmp05, wf, wrk17, gw);
   
   HELVertex amp06 (wrk12, wf, wrk17, gw);

   HELVertex amp07 (wrk07, wrk18, wrk02, glzn, grzn);

   HELVertex amp08 (wrk21, wrk18, zf   , glzn, grzn);

   HELVertex amp09 (wrk21, wrk11, wf   , glw , grw );

   HELVertex amp10 (wrk21, ep   , wrk15, glw , grw );

   HELVertex amp11 (nbf  , wrk18, wrk06, glzn, grzn);

   HELVertex amp12 (nbf  , wrk18, wrk23, glzn, grzn);

   HELVertex amp13 (wrk19, wrk09, wrk03, glw , grw );

   HELVertex amp14 (wrk20, wrk04, wf   , glw , grw );
   
   HELVertex amp15 (wrk20, wrk09, zf   , glzn, grzn);
   
   HELVertex amp16 (wrk20, emf  , wrk15, glw , grw );

   HELVector tmp17 (em   , wrk04, glae, grae, glze, grze, mz, gamz);
   HELVertex amp17 (wrk03, wf   , tmp17, gw);

   HELVector tmp18 (wrk19, emf  , glae, grae, glze, grze, mz, gamz);
   HELVertex amp18 (wrk03, wf   , tmp18, gw);

   HELVertex amp19a(wrk08, ep   , wrk05, glae, grae);
   HELVertex amp19b(wrk08, ep   , wrk06, glze, grze);
   Complex_t amp19 = amp19a + amp19b;

   HELVertex amp20a(wrk08, ep   , wrk22, glae, grae);
   HELVertex amp20b(wrk08, ep   , wrk23, glze, grze);
   Complex_t amp20 = amp20a + amp20b;

   HELVertex amp21 (wrk07, ep   , wrk10, glw, grw);

   HELVertex amp22 (nbf  , wrk11, wrk10, glw, grw);

   HELVertex amp23 (wrk03, wrk10, zf   , gwwz);

   //---------------------------
   // Fusion diagrams
   // (24-26)
   HELVertex amp24 (wrk03, wrk17, wf, zf, gw, gwwz, mw, gamw,kTRUE);

   // Higgs
   HELScalar tmp25 (wrk02, zf, gzzh, mh, gamh);
   HELVertex amp25 (wrk03, wf, tmp25, gwwh);

   //---------------------------
   // Sum up amplitudes
   //---------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 
	         + amp05 + amp06 + amp07 + amp08 
		 + amp09 + amp10 + amp11 + amp12
		 + amp13 + amp14 + amp15 + amp16 
		 + amp17 + amp18 + amp19 + amp20 
		 + amp21 + amp22 + amp23
		 + amp24
		 + amp25
		 ;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ENWZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ENWZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ENWZBases.BeamstrahlungFilename","trc500"));
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
  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWBosonPtr->GetMode(m); 
     if (mp && (m<fWModesLo || m>fWModesHi)) {
        mp->Lock();
     }
  }
  fWBosonPtr->DebugPrint();

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
  Double_t rs   = fEcmInit;
  Double_t s    = rs*rs;

  Double_t qwz  = kM_w + kM_z;
  Double_t m1   = kM_e;
  Double_t m2   = kM_e;

  Double_t xilo = TMath::Log(qwz*(qwz+m1+m2)/s);
  Double_t xihi = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);

  Double_t rxi  = TMath::Exp(xilo);
  Double_t dm1  = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1  = (m1*m1/s);

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;
  Double_t qwzlo = 0.; //2*mw - 40.;
  Double_t qwzhi = rs;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,  qwzlo,  qwzhi,       50, "qwz"   );
  Xh_init( 9,    -1.,    +1.,       50, "cosw"  );
  Xh_init(10,     0.,   k2Pi,       50, "phiw"  );
  Xh_init(11,    60.,   100.,       50, "Mw"    );
  Xh_init(12,    -1.,    +1.,       50, "CosF3" );
  Xh_init(13,     0.,   k2Pi,       50, "PhiF3" );
  Xh_init(14,    70.,   110.,       50, "Mz"    );
  Xh_init(15,    -1.,    +1.,       50, "CosF5" );
  Xh_init(16,     0.,   k2Pi,       50, "PhiF5" );
  Xh_init(17,     0.,     4.,        4, "Helin ");
  Xh_init(18,     0.,     4.,        4, "Helot ");
  Xh_init(19,     0.,    12.,       12, "W mode");
  Xh_init(20,     0.,    12.,       12, "Z mode");
  Xh_init(21,    -1.,    +1.,        2, "CP flg");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ENWZBases::Userout()
{
  cout << "End of ENWZBases----------------------------------- "  << endl
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
void ENWZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNf = 4;           // e-, e+, e-, nb, fu,fdb,  f, fb
   static const Int_t kHelComb[kNf][8] = {{-1, +1, -1, +1, -1, +1, -1, +1},
                                          {-1, +1, -1, +1, -1, +1, +1, -1},
                                          {+1, +1, +1, +1, -1, +1, -1, +1},
                                          {+1, +1, +1, +1, -1, +1, +1, -1}};

   // --------------------------------------------
   // Flip CP if needed
   // --------------------------------------------
   Double_t helcomb = fXComb;
   if (fXComb <= 0.5) {
      fCP = -1;
      helcomb = helcomb/0.5;
   } else {
      fCP = +1;
      helcomb = (helcomb - 0.5)/0.5;
   }
   weight = 2.;  

   // --------------------------------------------
   // Select helicities 
   // --------------------------------------------
   Double_t helm = (1. - fPolem)/2.;
   Double_t help = (1. - fPolep)/2.;
   if (fCP > 0) {
      if (helcomb < helm || helm == 1.) {
         Double_t helcombi = helcomb/helm; // (e-,e+)
         if (helcombi < help || help == 1.) {
            fJCombI = 0; // (-1,-1)
	    helcomb = helcombi/help;
         } else {
            fJCombI = 1; // (-1,+1)
	    helcomb = (helcombi-help)/(1.-help);
	 }
      } else {
         Double_t helcombi = (helcomb-helm)/(1.-helm);
         if (helcombi < help || help == 1.) {
            fJCombI = 2; // (+1,-1)
	    helcomb = helcombi/help;
	 } else {
            fJCombI = 3; // (+1,+1)
	    helcomb = (helcombi-help)/(1.-help);
	 }
      }
   } else {
      if (helcomb < helm || helm == 1.) {
         Double_t helcombi = helcomb/helm; // (e+,e-)
         if (helcombi < help || help == 1.) {
            fJCombI = 3; // (-1,-1)
	    helcomb = helcombi/help;
	 } else {
            fJCombI = 1; // (+1,-1)
	    helcomb = (helcombi-help)/(1.-help);
	 }
      } else {
         Double_t helcombi = (helcomb-helm)/(1.-helm);
         if (helcombi < help || help == 1.) {
            fJCombI = 2; // (-1,+1)
	    helcomb = helcombi/help;
	 } else {
            fJCombI = 0; // (+1,+1)
	    helcomb = (helcombi-help)/(1.-help);
	 }
      }
   }
#if 1
   if (fJCombI == 0 || fJCombI == 2) {
      weight = 0.;
      return;
   } 
#endif

   Int_t ioff = 0;
   if (fJCombI == 3) ioff = 2;

   fHelInitial[0] = kHelComb[ioff][0];
   fHelInitial[1] = kHelComb[ioff][1];

   fJCombF = (Int_t)(helcomb*(kNf/2));
   fJCombF = TMath::Min(fJCombF, (kNf/2)-1);
   fJCombF += ioff;

   fHelFinal  [0] = kHelComb[fJCombF][2];
   fHelFinal  [1] = kHelComb[fJCombF][3];
   fHelFinal  [2] = kHelComb[fJCombF][4];
   fHelFinal  [3] = kHelComb[fJCombF][5];
   fHelFinal  [4] = kHelComb[fJCombF][6];
   fHelFinal  [5] = kHelComb[fJCombF][7];
   weight *= (kNf/2);
}
