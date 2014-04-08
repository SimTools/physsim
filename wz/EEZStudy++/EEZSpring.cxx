//*****************************************************************************
//* =====================
//*  EEZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> EEZ generator
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "EEZSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(EEZSpring)
ClassImp(EEZSpringBuf)
ClassImp(EEZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class EEZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
EEZSpring::EEZSpring(const char     *name,
                     const char     *title,
                           EEZBases *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new EEZSpringBuf("EEZSpringBuf",
                               "EEZSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new EEZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
EEZSpring::~EEZSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t EEZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    EEZBases *bs = static_cast<EEZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> EEZBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class EEZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t EEZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  EEZBases     *bases   = static_cast<EEZBases *>(
                          static_cast<EEZSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 5;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // e-
  pv[1] = bases->fP[1]; // e+
  pv[3] = bases->fP[2]; // f
  pv[4] = bases->fP[3]; // fb
  pv[2] = pv[3] + pv[4];

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
  fQ2Z          = bases->GetQ2Z ();
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
  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    ide     = 11;                          // PDG code for e
  Double_t me      = kM_e;                        // electron mass
  Double_t qe      = -1;                          // electron charge
  Int_t    idf     = bases->f3Ptr->GetPID   ();   // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge();   // f charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f mass
  Double_t m4      = m3;                          // f mass
  Int_t    hel1    = bases->fHelFinal[2];         // f helicity
  Int_t    hel2    = bases->fHelFinal[3];         // fb helicity
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	  // shower level
  Int_t    icf     = 2;                           // color flux id
  Double_t rq2z    = pv[2].Mag();
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

  //                              No. PID   Mass  Charge   pv    Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, ide,   me,    qe, *qp[0], 0, 0,  0,    0,   0,     0);
  new (partons[1]) JSFSpringParton(2,-ide,   me,   -qe, *qp[1], 0, 0,  0,    0,   0,     0);
  new (partons[2]) JSFSpringParton(3, idz, rq2z,    0., *qp[2], 0, 0,  0,    0,   0,     0);
  new (partons[3]) JSFSpringParton(4, idf,   m3,  chrg, *qp[3], 0, 0,  2, hel1, icf, islev);
  new (partons[4]) JSFSpringParton(5,-idf,   m4, -chrg, *qp[4], 0, 0,  2, hel2, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class EEZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
EEZBases::EEZBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPolem     (0.),
           fPolep     (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fZBosonPtr (0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2Z       (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
	   fCosF      (0.),
	   fPhiF      (0.),
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

  cout << "Init eezbases " << endl;
  
  using namespace std;

  stringstream ins(gJSF->Env()->GetValue("EEZBases.Ecm","500.")); // E_cm [GeV]
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.Polem","0."));       // electron polarization
  ins >> fPolem;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.Polep","0."));       // positron polarization
  ins >> fPolep;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("EEZBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2Z  , 0., 1., 1, 1);

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
  SetNoOfSample(80000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 20);
  SetIteration2(0.05, 50);

}
// --------------------------
//  D-tor
// --------------------------
EEZBases::~EEZBases()
{
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t EEZBases::Func()
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
  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = kM_e; // Me
  Double_t m2   = kM_e; // Meb
  if (fEcmIP < m1 + m2 + m3 + m4) {
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

  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - (m1 + m2);
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t q3   = TMath::Sqrt(fQ2Z);       // Mz
  Double_t s    = rs*rs;
#if 1
  Double_t xilo = TMath::Log(q3*(q3+m1+m2)/s);
  Double_t xihi = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);

  fXi = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t rxi  = TMath::Exp(fXi);
#else
  Double_t xilo = 1/(1.-2.*TMath::Min(m1,m2)/rs);
  Double_t xihi = s/(q3*(q3+m1+m2));

  fXi = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t rxi  = 1./fXi;
#endif
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
  Xh_fill( 8, q3               , (bsWeight*sigma));
  Xh_fill( 9, fCosF            , (bsWeight*sigma));
  Xh_fill(10, fPhiF            , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill(12, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill(13, (Double_t)fZMode,  (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t EEZBases::DSigmaDX()
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs  = fEcmIP;
  Double_t s   = rs*rs; 
  Double_t eb  = rs/2.;
  Double_t q32 = fQ2Z;
  Double_t q3  = TMath::Sqrt(q32); 
#if 1
  Double_t rxi = TMath::Exp(fXi);
  weight *= rxi;
#else
  Double_t rxi = 1/fXi;
  weight *= rxi*rxi;
#endif

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
    Double_t rximn = (q3-m1+m2)*(q3+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - q32/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - q32/s);
  } else {
    Double_t rximn = (q3+m1-m2)*(q3+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - q32/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - q32/s);
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
  ANL4DVector pz = qcm - fP[0] - fP[1];
  GENFrame    cmframe;

  Double_t m3  = fM[2]; // Mf
  Double_t m4  = fM[3]; // Mfb
  Double_t m32 = m3*m3;
  Double_t m42 = m4*m4;
  GENPhase2 phaseZ(pz, m32, m42, cmframe, fCosF, fPhiF, 1);
  fP[2] = phaseZ.GetFourMomentum(0);
  fP[3] = phaseZ.GetFourMomentum(1);
  Double_t betaf = phaseZ.GetBetaBar();
  if (betaf <= 0.) return 0.;

  fSh1 = sh1;
  fCh1 = ch1;
  fSh2 = sh2;
  fCh2 = ch2;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
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
  Double_t dPhase = kFact * betaf * weight;        // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor

  Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                       // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t EEZBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t EEZBases::FullAmplitude()
{
   Double_t mz     = fZBosonPtr->GetMass();
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf     = f3Ptr->GetCharge();
   Double_t t3f    = f3Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e,  fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e,  fHelInitial[1], -1, kIsOutgoing);

   HELFermion e (fP[0], fM[0], fHelFinal  [0], +1, kIsOutgoing);
   HELFermion eb(fP[1], fM[1], fHelFinal  [1], -1, kIsIncoming);

   HELFermion f (fP[2], fM[2], fHelFinal  [2], +1, kIsOutgoing);
   HELFermion fb(fP[3], fM[3], fHelFinal  [3], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, mz, gamz);

   Complex_t amp = AmpEEtoEEZ(em, ep, e, eb, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoEEZ()
// --------------------------
Complex_t EEZBases::AmpEEtoEEZ(const HELFermion &em,
                               const HELFermion &ep,
                               const HELFermion &e,
                               const HELFermion &eb,
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
   // Z Production Amplitude
   //---------------------------
   Double_t  me    = kM_e;
   Double_t  game  = 0.;
   Double_t  mz    = fZBosonPtr->GetMass();
   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);
   Double_t  glae  = -kGe*qe;
   Double_t  grae  = -kGe*qe;

   Double_t  ebm   = TMath::Abs(em.GetFourMomentum()(0));
   Double_t  e1    = TMath::Abs(e .GetFourMomentum()(0));
   Double_t  e2    = TMath::Abs(eb.GetFourMomentum()(0));

   HELVector a1(ebm, e1, fSh1, fCh1, fPhi1, 
		fHelInitial[0], fHelFinal[0],+1, kGe, me);
   HELVector z1(em, e , glze, grze, mz, gamz);

   HELVector a2(ebm, e2, fSh2, fCh2, fPhi2, 
		fHelInitial[1], fHelFinal[1],-1, kGe, me);
   HELVector z2(eb, ep, glze, grze, mz, gamz);

   HELFermion ebv(eb, zf, glze, grze, me, game);
   Complex_t amp1 = HELVertex(ebv, ep, a1, glae, grae);
   Complex_t amp2 = HELVertex(ebv, ep, z1, glze, grze);

   HELFermion epv(ep, zf, glze, grze, me, game);
   Complex_t amp3 = HELVertex(eb, epv, a1, glae, grae);
   Complex_t amp4 = HELVertex(eb, epv, z1, glze, grze);

   HELFermion ev(e, zf, glze, grze, me, game);
   Complex_t amp5 = HELVertex(em, ev, a2, glae, grae);
   Complex_t amp6 = HELVertex(em, ev, z2, glze, grze);

   HELFermion emv(em, zf, glze, grze, me, game);
   Complex_t amp7 = HELVertex(emv, e, a2, glae, grae);
   Complex_t amp8 = HELVertex(emv, e, z2, glze, grze);

#if 1
   Complex_t amp = amp1 + amp2 + amp3 + amp4
	         + amp5 + amp6 + amp7 + amp8; 
#else
   Complex_t amp = amp1 + amp3 + amp5 + amp7; 
#endif
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void EEZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("EEZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("EEZBases.BeamstrahlungFilename","trc500"));
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
  if (! fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
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

  Double_t q3   = kM_z;
  Double_t m1   = kM_e;
  Double_t m2   = kM_e;

  Double_t xilo = TMath::Log(q3*(q3+m1+m2)/s);
  Double_t xihi = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);

  Double_t rxi  = TMath::Exp(xilo);
  Double_t dm1  = (m1*m1/s)*rxi*rxi/(1.-rxi);
  Double_t dp1  = (m1*m1/s);

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,    65.,   115.,       50, "mz"    );
  Xh_init( 9,    -1.,    +1.,       50, "cosf"  );
  Xh_init(10,     0.,   k2Pi,       50, "phif"  );
  Xh_init(11,     0.,     4.,        4, "Helin ");
  Xh_init(12,     0.,     8.,        8, "Helot ");
  Xh_init(13,     0.,    12.,       12, "z mode ");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void EEZBases::Userout()
{
  cout << "End of EEZBases----------------------------------- "  << endl
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
void EEZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 4;
   static const Int_t kIHelComb[kNi][2] = {{-1, -1},
                                           {-1, +1},
                                           {+1, -1},
                                           {+1, +1}};
   static const Int_t kNf = 8;
   static const Int_t kFHelComb[kNf][4] = {{-1, -1, -1, +1},
                                           {-1, -1, +1, -1},
                                           {-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1},
                                           {+1, +1, -1, +1},
                                           {+1, +1, +1, -1}};
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
   weight = kNf;
}
