//*****************************************************************************
//* =====================
//*  NNHHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> NNHH generator
//*
//* (Update Record)
//*    2010/04/22  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "NNHHSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(NNHHSpring)
ClassImp(NNHHSpringBuf)
ClassImp(NNHHBases)

//-----------------------------------------------------------------------------
// ==============================
//  class NNHHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNHHSpring::NNHHSpring(const char      *name,
                       const char      *title,
                             NNHHBases *bases)
           : JSFSpring(name, title, bases)
{
  fEventBuf = new NNHHSpringBuf("NNHHSpringBuf",
                                "NNHHSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new NNHHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
NNHHSpring::~NNHHSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t NNHHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    NNHHBases *bs = static_cast<NNHHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> NNHHBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class NNHHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t NNHHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  NNHHBases    *bases   = static_cast<NNHHBases *>(
                          static_cast<NNHHSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 4;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // ne
  pv[1] = bases->fP[1]; // neb
  pv[2] = bases->fP[2]; // h1
  pv[3] = bases->fP[3]; // h2

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
  fQ2HH         = bases->GetQ2HH();
  fCosH         = bases->GetCosH();
  fPhiH         = bases->GetPhiH();
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
  Int_t    idne    = 12;                          // PDG code for nu_e
  Double_t mass    = bases->GetMass();
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
  new (partons[0]) JSFSpringParton(1, idne,   0.,    0., *qp[0], 0, 0,  0,    0,   0,     0);
  new (partons[1]) JSFSpringParton(2,-idne,   0.,    0., *qp[1], 0, 0,  0,    0,   0,     0);
  new (partons[2]) JSFSpringParton(3, idh , mass,    0., *qp[2], 0, 0,  0,    0,   0,     0);
  new (partons[3]) JSFSpringParton(4, idh , mass,    0., *qp[3], 0, 0,  0,    0,   0,     0);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class NNHHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNHHBases::NNHHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 120.),
           fWidth     ( 0.006),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fWBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2HH      (0.),
	   fCosH      (0.),
	   fPhiH      (0.),
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
  stringstream ins(gJSF->Env()->GetValue("NNHHBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.WidthH","0.006")); // M_x [GeV]
  ins >> fWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNHHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1
  //--
  DefineVariable(fXXi   , 0., 1., 1, 1);
  DefineVariable(fXEta1 , 0., 1., 1, 1);
  DefineVariable(fXEta2 , 0., 1., 1, 1);
  DefineVariable(fXPhi1 , 0., 1., 1, 1);
  DefineVariable(fXPhi21, 0., 1., 1, 1);
  DefineVariable(fXQ2HH , 0., 1., 1, 1);
  DefineVariable(fXCosH , 0., 1., 0, 1);
  DefineVariable(fXPhiH , 0., 1., 0, 1);

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
  SetIteration2(0.05, 40);

}
// --------------------------
//  D-tor
// --------------------------
NNHHBases::~NNHHBases()
{
  delete fWBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t NNHHBases::Func()
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
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = 0.; // Mne
  Double_t m2   = 0.; // Mneb
  if (fEcmIP < fMass + fMass + m1 + m2) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average for e- is taken
  //  care of here
  Double_t weight = 1;
  SelectHelicities(weight);
  if (weight == 0.) return 0.;
  bsWeight *= weight;

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t rs     = fEcmIP;
  Double_t s      = rs*rs;
  Double_t qhh2mn = TMath::Power(fMass + fMass,2);
  Double_t qhh2mx = TMath::Power(rs - (m1+m2),2);
  fQ2HH     = qhh2mn + (qhh2mx-qhh2mn)*fXQ2HH;
  bsWeight *= qhh2mx - qhh2mn;

  Double_t qhh   = TMath::Sqrt(fQ2HH);
  Double_t xilo  = TMath::Log(qhh*(qhh+m1+m2)/s);
  Double_t xihi  = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);
  fXi       = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t mw   = fWBosonPtr->GetMass();
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

  fCosH  = -1. + 2.*fXCosH;
  fPhiH  = k2Pi * fXPhiH;
  bsWeight *= 2.*k2Pi;

  fM[0] = m1;
  fM[1] = m2;
  fM[2] = fMass;
  fM[3] = fMass;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX();

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
  Xh_fill( 8, TMath::Sqrt(fQ2HH), (bsWeight*sigma));
  Xh_fill( 9, fCosH             , (bsWeight*sigma));
  Xh_fill(10, fPhiH             , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(12, (Double_t)fJCombF , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t NNHHBases::DSigmaDX()
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs   = fEcmIP;
  Double_t s    = rs*rs; 
  Double_t eb   = rs/2.;
  Double_t qhh2 = fQ2HH; // q^2_{HH}
  Double_t m3   = fM[2]; // Mh
  Double_t m4   = fM[3]; // Mh
  Double_t qhh  = TMath::Sqrt(qhh2);
  Double_t rxi  = TMath::Exp(fXi);
  weight *= rxi;

  Double_t m1  = fM[0]; // Mne
  Double_t m2  = fM[1]; // Mneb
  Double_t mw  = fWBosonPtr->GetMass();
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
    Double_t rximn = (qhh-m1+m2)*(qhh+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - qhh2/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - qhh2/s);
  } else {
    Double_t rximn = (qhh+m1-m2)*(qhh+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - qhh2/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - qhh2/s);
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
  ANL4DVector qvhh = qcm - fP[0] - fP[1];
  GENFrame  cmframe;
  GENPhase2 phaseHH(qvhh, m3*m3, m4*m4, cmframe, fCosH, fPhiH, 1);
  fP[2] = phaseHH.GetFourMomentum(0);
  fP[3] = phaseHH.GetFourMomentum(1);
  Double_t betah = phaseHH.GetBetaBar();
  if (betah <= 0.) return 0.;
  weight *= betah;

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

  Double_t identp = 1./2.;                         // identical particle factor
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
Double_t NNHHBases::AmpSquared()
{
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2);

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t NNHHBases::FullAmplitude()
{
   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  h1(fP[2]);
   HELScalar  h2(fP[3]);

   HELFermion ne (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing);
   HELFermion neb(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);

   Complex_t amp = AmpEEtoNNHH(em, ep, ne, neb, h1, h2);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNHH()
// --------------------------
Complex_t NNHHBases::AmpEEtoNNHH(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &ne,
                                 const HELFermion &neb,
                                 const HELScalar  &h1,
                                 const HELScalar  &h2)
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
   Double_t mw     = fWBosonPtr->GetMass();
   Double_t gamw   = fWBosonPtr->GetWidth();
   Double_t glwf   = -kGw*kSqh;
   Double_t grwf   = 0.;

   HELVector w1(em , ne, glwf, grwf, mw, gamw);
   HELVector w2(neb, ep, glwf, grwf, mw, gamw);

   Double_t mh     = fMass;
   Double_t gamh   = fWidth;
   Double_t gwwh   = kGw*mw;
   Double_t gwwhh  = kGw*kGw/2.;
   Double_t v      = 2*mw/kGw;
   Double_t ghhh   = -3.*mh*mh/v;

   // (1) self-coupling
   HELScalar hh(h1, h2, ghhh, mh, gamh);
   Complex_t amp1 = HELVertex(w1, w2, hh, gwwh);

   // (2) 4-point
   Complex_t amp2 = HELVertex(w1, w2, h1, h2, gwwhh);

   // (3) h1 and h2 from W
   HELVector w1h1(w1, h1, gwwh, mw, gamw);
   Complex_t amp3 = HELVertex(w1h1, w2, h2, gwwh);
   HELVector w1h2(w1, h2, gwwh, mw, gamw);
             amp3 += HELVertex(w1h2, w2, h1, gwwh);
   Complex_t amp = amp1 + amp2 + amp3;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void NNHHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("NNHHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("NNHHBases.BeamstrahlungFilename","trc500"));
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
  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  fWBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs    = fEcmIP;
  Double_t s     = rs*rs;
  Double_t xilo  = TMath::Log(fMass*fMass/s);
  Double_t xihi  = 0.;
  Double_t mw    = fWBosonPtr->GetMass();
  Double_t dm1   = mw*mw/s;
  Double_t dp1   = 1.;
  Double_t etlo  = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi  =  TMath::Log((1.+dp1)/dm1)/2.;
  Double_t qhhlo = 2*fMass - 40.;
  Double_t qhhhi = rs;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,  qhhlo,  qhhhi,       50, "qhh"   );
  Xh_init( 9,    -1.,    +1.,       50, "cosh"  );
  Xh_init(10,     0.,   k2Pi,       50, "phih"  );
  Xh_init(11,     0.,     2.,        2, "Helin ");
  Xh_init(12,     0.,     2.,        2, "Helot ");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void NNHHBases::Userout()
{
  cout << "End of NNHHBases----------------------------------- "  << endl
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
void NNHHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, 0, 0}};
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
