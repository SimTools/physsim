//*****************************************************************************
//* =====================
//*  ENWSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ENW generator
//*
//* (Update Record)
//*    2014/09/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ENWSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ENWSpring)
ClassImp(ENWSpringBuf)
ClassImp(ENWBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ENWSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ENWSpring::ENWSpring(const char     *name,
                     const char     *title,
                           ENWBases *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new ENWSpringBuf("ENWSpringBuf",
                               "ENWSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new ENWBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ENWSpring::~ENWSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ENWSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ENWBases *bs = static_cast<ENWBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ENWBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ENWSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ENWSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ENWBases     *bases   = static_cast<ENWBases *>(
                          static_cast<ENWSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 5;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // e-  / e+
  pv[1] = bases->fP[1]; // neb / ne
  pv[3] = bases->fP[2]; // fu  / fub
  pv[4] = bases->fP[3]; // fdb / fd
  pv[2] = pv[3] + pv[4];
  for (Int_t i=0; i<fNparton; i++) pv[i] *= ((Double_t)fCP);

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------
  fCP           = bases->GetCP();
  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fXi           = bases->GetXi();
  fEta1         = bases->GetEta1();
  fEta2         = bases->GetEta2();
  fPhi1         = bases->GetPhi1();
  fPhi2         = bases->GetPhi2();
  fQ2W          = bases->GetQ2W ();
  fCosF         = bases->GetCosThetaF();
  fPhiF         = bases->GetPhiF();
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
  Int_t    idw     = fCP*24;                          // PDG code for W+
  Int_t    ide     = fCP*11;                          // PDG code for e-
  Int_t    idn     = fCP*12;                          // PDG code for e-
  Double_t mn      = 0.;                              // electron mass
  Double_t me      = kM_e;                            // electron mass
  Double_t chrgw   = fCP;
  Double_t chrge   = -fCP;                            // electron charge

  Int_t    idf3    = fCP*bases->f3Ptr->GetPID   ();   // PDG code for fu
  Double_t chrg3   = fCP*bases->f3Ptr->GetCharge();   // fu charge
  Double_t m3      =     bases->f3Ptr->GetMass  ();   // fu mass
  Int_t    hel3    = fCP*bases->fHelFinal[2];         // fu helicity
  Int_t    idf4    = fCP*bases->f4Ptr->GetPID   ();   // PDG code for fdb
  Double_t chrg4   = fCP*bases->f4Ptr->GetCharge();   // fdb charge
  Double_t m4      =     bases->f4Ptr->GetMass  ();   // fdb mass
  Int_t    hel4    = fCP*bases->fHelFinal[3];         // fdb helicity
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 101 : 0;  	  // shower level
  Int_t    icf     = 1;                           // color flux id
  Double_t rq2w    = pv[2].Mag();
#ifdef __DEBUG__
  cerr << endl;
  ANL4DVector qcm;
  for (Int_t ip=0; ip<5; ip++) {
    cerr << "pv[" << ip << "] = (" 
         <<  pv[ip].E() << ", "
         <<  pv[ip].Px() << ", "
         <<  pv[ip].Py() << ", "
         <<  pv[ip].Pz() << ") "
	 << "m = " << pv[ip].GetMass() << endl;
    if (ip<3) qcm += pv[ip];
  }
  cerr << "qcm = ("
       <<  qcm.E() << ", "
       <<  qcm.Px() << ", "
       <<  qcm.Py() << ", "
       <<  qcm.Pz() << ") " << endl;
  cerr << "----" << endl;
#endif

  //                               No. PID   Mass  Charge   pv   Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, ide ,   me, chrge, *qp[0], 0, 0,  0,    0,   0,     0);
  new (partons[1]) JSFSpringParton(2,-idn ,   mn,    0., *qp[1], 0, 0,  0,    0,   0,     0);
  new (partons[2]) JSFSpringParton(3, idw , rq2w, chrgw, *qp[2], 2, 4,  0,    0,   0,     0);
  new (partons[3]) JSFSpringParton(4, idf3,   m3, chrg3, *qp[3], 0, 0,  3, hel3, icf, islev);
  new (partons[4]) JSFSpringParton(5, idf4,   m4, chrg4, *qp[4], 0, 0,  3, hel4, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ENWBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ENWBases::ENWBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPolem     (0.),
           fPolep     (0.),
           fWModesLo  ( 1),
           fWModesHi  (12),
	   fNCALL     (80000),
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
	   fQ2W       (0.),
           fWModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
	   fCosF      (0.),
	   fPhiF      (0.),
           fCP        ( 1),
	   fSh1       (0.),
	   fCh1       (0.),
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

  cout << "Init enwbases " << endl;
  
  using namespace std;

  stringstream ins(gJSF->Env()->GetValue("ENWBases.Ecm","500.")); // E_cm [GeV]
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.Polem","0."));       // electron polarization
  ins >> fPolem;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.Polep","0."));       // positron polarization
  ins >> fPolep;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.WModesLo","1"));      // Z decay mode lo
  ins >> fWModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.WModesHi","12"));     // Z decay mode hi
  ins >> fWModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ENWBases.NCALL","80000"));
  ins >> fNCALL;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fXCP           , 0., 1., 1, 1);
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fWDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2W          , 0., 1., 1, 1);

  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1
  //--
  DefineVariable(fXXi   , 0., 1., 1, 1);
  DefineVariable(fXEta1 , 0., 1., 1, 1);
  DefineVariable(fXEta2 , 0., 1., 1, 1);
  DefineVariable(fXPhi1 , 0., 1., 1, 1);
  DefineVariable(fXPhi21, 0., 1., 1, 1);

  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();

  DefineVariable(fCosF , fXL[0], fXU[0], 0, 1);
  DefineVariable(fPhiF , fXL[1], fXU[1], 0, 1);

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
ENWBases::~ENWBases()
{
  delete fWBosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ENWBases::Func()
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

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = kM_e; // Me
  Double_t m2   = 0.;   // Mneb
  if (fEcmIP < m1 + m2 + m3 + m4) {
    return 0.;
  }

  // --------------------------------------------
  //  Flip CP
  // --------------------------------------------
  if (fXCP <= 0.5) fCP = -1;
  else             fCP = +1;
  bsWeight *= 2.;  

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

  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - (m1 + m2);
#ifndef __ZEROWIDTH__
  fQ2W = fWBosonPtr->GetQ2BW(qmin, qmax, fXQ2W, weight);
#else
  fQ2W = TMath::Power(fWBosonPtr->GetMass(),2);
  weight = kPi*fWBosonPtr->GetMass()*fWBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t qw   = TMath::Sqrt(fQ2W);       // Mw
  Double_t s    = rs*rs;

  Double_t xilo = TMath::Log((qw-m1)*(qw+m1)/s);
  Double_t xihi = TMath::Log(1.-2.*m1/rs);
           xilo = TMath::Min(xihi,xilo);
  fXi = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;
  Double_t rxi  = TMath::Exp(fXi);

  Double_t dme  = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dpe  = 1.;

  Double_t etlo = -TMath::Log((1.+dme)/dpe)/2.;
  Double_t ethi =  TMath::Log((1.+dpe)/dme)/2.;
  fEta1     = etlo + (ethi-etlo)*fXEta1;
  bsWeight *= ethi - etlo;

  fEta2  = fXEta2;

  fPhi1  = k2Pi * fXPhi1;
  fPhi2  = k2Pi * fXPhi21 + fPhi1;
  fPhi2  = fPhi2 > k2Pi ? fPhi2 - k2Pi : fPhi2;
  bsWeight *= k2Pi*k2Pi;

  fM[0] = m1;
  fM[1] = m2;
  fM[2] = m3;
  fM[3] = m4;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX() * qed;

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP           , (bsWeight*sigma));
  Xh_fill( 2, fXi              , (bsWeight*sigma));
  Xh_fill( 3, fEta1            , (bsWeight*sigma));
  Xh_fill( 4, fEta2            , (bsWeight*sigma));
  Xh_fill( 5, fPhi1            , (bsWeight*sigma));
  Xh_fill( 6, fPhi2            , (bsWeight*sigma));
  Xh_fill( 7, k2Pi*fXPhi21     , (bsWeight*sigma));
  Xh_fill( 8, qw               , (bsWeight*sigma));
  Xh_fill( 9, fCosF            , (bsWeight*sigma));
  Xh_fill(10, fPhiF            , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill(12, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill(13, (Double_t)fWMode , (bsWeight*sigma));
  Xh_fill(14, (Double_t)fCP    , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ENWBases::DSigmaDX()
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs  = fEcmIP;
  Double_t s   = rs*rs; 
  Double_t eb  = rs/2.;
  Double_t qw2 = fQ2W;
  Double_t qw  = TMath::Sqrt(qw2); 
  Double_t m1  = fM[0]; // Me
  Double_t m2  = fM[1]; // Mneb

  Double_t rxi  = TMath::Exp(fXi);
  Double_t qwn2 = s*rxi + m1*m1;
  Double_t qwn  = TMath::Sqrt(qwn2);
  weight *= s*rxi;

  Double_t dm1 = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1 = 1.;

  Double_t exp2et1 = TMath::Exp(2.*fEta1);
  Double_t sh1     = TMath::Min(TMath::Sqrt((1.+dp1+dm1)/(1.+exp2et1)   - dm1), 1.);
  Double_t ch1     = TMath::Sqrt((1-sh1)*(1+sh1));
  static Double_t kSqrt2 = TMath::Sqrt(2.);
  Double_t cs1     = (1.-kSqrt2*sh1)*(1.+kSqrt2*sh1);
  Double_t sn1     = 2.*sh1*ch1;
  Double_t fi1     = fPhi1;
  weight *= 4.* (1.+dp1+dm1)/((1.+exp2et1)*(1.+1./exp2et1));

  Double_t dlt = 10.;
  if (qw > 0. && qwn > qw) {
    dlt = 2.*qw2/((qwn-qw)*(qwn+qw));
  }

  Double_t opdzt = TMath::Power(1. + 1./dlt, fEta2);
  Double_t sh2   = TMath::Min(TMath::Sqrt(dlt*(opdzt - 1.)), 1.);
  Double_t cs2   = (1.-kSqrt2*sh2)*(1.+kSqrt2*sh2);
  Double_t fi2   = fPhi2;
  weight *= 2.*dlt*opdzt*TMath::Log(1.+1./dlt);

  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

  GENFrame    cmframe;
  ANL4DVector qcm(rs, 0.,0.,0.);
  Double_t    m12 = m1*m1;
  GENPhase2   phaseCM(qcm, m12, qwn2, cmframe, cs1, fi1, 0);
              fP[0] = phaseCM.GetFourMomentum(0);
  ANL4DVector pwn   = phaseCM.GetFourMomentum(1);
  Double_t    betae = phaseCM.GetBetaBar();
  if (betae <= 0.) return 0.;

  Double_t ap1 = fP[0].P();
  fP[0].SetX(ap1*sn1*TMath::Cos(fi1));
  fP[0].SetY(ap1*sn1*TMath::Sin(fi1));

  Double_t    m22 = m2*m2;
  GENPhase2   phaseWN(pwn, m22, qw2, cmframe, cs2, fi2, 1);
              fP[1] = phaseWN.GetFourMomentum(0);
  ANL4DVector pw    = phaseWN.GetFourMomentum(1);
  Double_t    betan = phaseWN.GetBetaBar();
  if (betan <= 0.) return 0.;

  Double_t m3  = fM[2]; // Mfu
  Double_t m4  = fM[3]; // Mfdb
  Double_t m32 = m3*m3;
  Double_t m42 = m4*m4;
  GENFrame wnframe = phaseWN.GetFrame();
  GENPhase2 phaseW(pw, m32, m42, wnframe, fCosF, fPhiF, 1);
  fP[2] = phaseW.GetFourMomentum(0);
  fP[3] = phaseW.GetFourMomentum(1);
  Double_t betaf = phaseW.GetBetaBar();
  if (betaf <= 0.) return 0.;

  fSh1 = sh1;
  fCh1 = ch1;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
#ifdef __DEBUG__
  for (Int_t i=0; i<4; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qvw = fP[2] + fP[3];
  cerr << " pw  = (" 
       << qvw.E () << ","
       << qvw.Px() << ","
       << qvw.Py() << ","
       << qvw.Pz() << ")" << endl;
  ANL4DVector pcm = qvw + fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
  cerr << " wmass = " << qvw.GetMass() << endl;
#endif
  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();
  
  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betae * betan * betaf 
                          * weight;                // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor

  Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                       // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ENWBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor();
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
Complex_t ENWBases::FullAmplitude()
{
   Double_t mw     = fWBosonPtr->GetMass();
   Double_t gamw   = fWBosonPtr->GetWidth();
   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em (fK[0], kM_e,  fHelInitial[0], +1, kIsIncoming);
   HELFermion ep (fK[1], kM_e,  fHelInitial[1], -1, kIsOutgoing);

   HELFermion e  (fP[0], fM[0], fHelFinal  [0], +1, kIsOutgoing);
   HELFermion neb(fP[1], fM[1], fHelFinal  [1], -1, kIsIncoming);

   HELFermion fu (fP[2], fM[2], fHelFinal  [2], +1, kIsOutgoing);
   HELFermion fdb(fP[3], fM[3], fHelFinal  [3], -1, kIsIncoming);
   HELVector  wf(fdb, fu, glw, grw, mw, gamw);

   Complex_t amp = AmpEEtoENW(em, ep, e, neb, wf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoENW()
// --------------------------
Complex_t ENWBases::AmpEEtoENW(const HELFermion &em,
                               const HELFermion &ep,
                               const HELFermion &emf,
                               const HELFermion &nbf,
                               const HELVector  &wf)
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
   // Z Production Amplitude
   //---------------------------
   Double_t  mn    = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamn  = 0.;
   Double_t  game  = 0.;

   Double_t  mw    = fWBosonPtr->GetMass();
   Double_t  gamw  = fWBosonPtr->GetWidth();
   Double_t  mz    = fZBosonPtr->GetMass();
   Double_t  gamz  = fZBosonPtr->GetWidth();

   Double_t  glwf  = -kGw*kSqh;
   Double_t  grwf  = 0.;

   Double_t  qe    = -1.;
   Double_t  glae  = -kGe*qe;
   Double_t  grae  = -kGe*qe;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  t3n   = +1./2.;
   Double_t  glzn  = -kGz*t3n;
   Double_t  grzn  = 0.;

   Double_t  gwwa  = kGe;     
   Double_t  gwwz  = kGw*kCosW;

   Double_t  ebm   = TMath::Abs(em .GetFourMomentum()(0));
   Double_t  e1    = TMath::Abs(emf.GetFourMomentum()(0));

   HELVector  wrk01 (ebm, e1, fSh1, fCh1, fPhi1, 
                    fHelInitial[0], fHelFinal[0],+1, kGe, me);
   HELVector  wrk02 (em , emf, glze, grze, mz, gamz);
   HELVector  wrk03 (nbf, ep , glwf, grwf, mw, gamw);

   HELFermion wrk04 (ep , wf, glwf, grwf, mn, gamn);
   HELFermion wrk05 (nbf, wf, glwf, grwf, me, game);
   HELFermion wrk06 (emf, wf, glwf, grwf, mn, gamn);

   HELVertex  amp01 (wrk05, ep, wrk01, glae, grae);
   HELVertex  amp02 (wrk05, ep, wrk02, glze, grze);

   HELVertex  amp03 (wrk03, wf, wrk01, gwwa);
   HELVertex  amp04 (wrk03, wf, wrk02, gwwz);

   HELVertex  amp05 (nbf, wrk04, wrk02, glzn, grzn);

   HELVertex  amp06 (em, wrk06, wrk03, glwf, grwf);

   Complex_t amp = amp01 + amp02 + amp03 + amp04 + amp05 + amp06;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ENWBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ENWBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ENWBases.BeamstrahlungFilename","trc500"));
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
  if (! fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWBosonPtr->GetMode(m); 
     if (mp && (m<fWModesLo || m>fWModesHi)) {
        mp->Lock();
     }
  }
  fWBosonPtr->DebugPrint();

  if (! fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  fZBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs   = fEcmInit;
  Double_t s    = rs*rs;

  Double_t qw   = kM_w;
  Double_t m1   = kM_e;

  Double_t xilo = TMath::Log((qw-m1)*(qw+m1)/s);
  Double_t xihi = TMath::Log(1.-2.*m1/rs);

  Double_t rxi  = TMath::Exp(xilo);
  Double_t dm1 = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1 = 1.;

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,     0.,     1.,       50, "zeta"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,    55.,   105.,       50, "mw"    );
  Xh_init( 9,    -1.,    +1.,       50, "cosf"  );
  Xh_init(10,     0.,   k2Pi,       50, "phif"  );
  Xh_init(11,     0.,     4.,        4, "Helin ");
  Xh_init(12,     0.,     2.,        2, "Helot ");
  Xh_init(13,     0.,    12.,       12, "w mode");
  Xh_init(14,    -1.,    +1.,        2, "CP flg");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ENWBases::Userout()
{
  cout << "End of ENWBases----------------------------------- "  << endl
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
void ENWBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 4;
   static const Int_t kIHelComb[kNi][2] = {{-1, -1}, // zero
                                           {-1, +1},
                                           {+1, -1}, // zero
                                           {+1, +1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, -1, +1},
                                           {+1, +1, -1, +1}};
   Double_t helm = (1. - fPolem)/2.;
   Double_t help = (1. - fPolep)/2.;
   if (fCP > 0) {
      if (fHelCombInitial < helm || helm == 1.) {
         Double_t helcombi = fHelCombInitial/helm;
         if (helcombi < help) fJCombI = 0;
         else                 fJCombI = 1;
      } else {
         Double_t helcombi = (fHelCombInitial-helm)/(1.-helm);
         if (helcombi < help) fJCombI = 2;
         else                 fJCombI = 3;
      }
   } else {
      if (fHelCombInitial < helm || helm == 1.) {
         Double_t helcombi = fHelCombInitial/helm;
         if (helcombi < help) fJCombI = 3;
         else                 fJCombI = 2;
      } else {
         Double_t helcombi = (fHelCombInitial-helm)/(1.-helm);
         if (helcombi < help) fJCombI = 1;
         else                 fJCombI = 0;
      }
   }
   if (fJCombI == 0 || fJCombI == 2) {
      weight = 0.;
      return;
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
