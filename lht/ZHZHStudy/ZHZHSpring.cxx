//*****************************************************************************
//* =====================
//*  ZHZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZH ZH generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/12/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZHZHSpring.h"

#include <sstream>
#include <iomanip>
//#define __NODECAY__
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PAHSESPACE__
#ifdef __NODECAY__
#ifndef __ZEROWIDTH__
#define __ZEROWIDTH__
#endif
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZHZHSpring)
ClassImp(ZHZHSpringBuf)
ClassImp(ZHZHBases)
ClassImp(ZHBoson)
ClassImp(AHBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class ZHZHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHZHSpring::ZHZHSpring(const char      *name,
                       const char      *title,
                             ZHZHBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new ZHZHSpringBuf("ZHZHSpringBuf",
                                "ZHZHSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new ZHZHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZHZHSpring::~ZHZHSpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZHZHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZHZHBases *bs = static_cast<ZHZHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZHZHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHZHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZHZHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZHZHBases  *bases     = static_cast<ZHZHBases *>
                         (static_cast<ZHZHSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 6;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];    // h  from ZH1
  pv[3] = bases->fP[1];    // AH from ZH1
  pv[4] = bases->fP[2];    // h  from ZH2
  pv[5] = bases->fP[3];    // AH from ZH2
  pv[0] = pv[2] + pv[3];   // ZH1
  pv[1] = pv[4] + pv[5];   // ZH2

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2XX         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X1         = bases->GetQ2X1();
  fCosThetaH1   = bases->GetCosThetaH1();
  fPhiH1        = bases->GetPhiH1();
  fQ2X2         = bases->GetQ2X2();
  fCosThetaH2   = bases->GetCosThetaH2();
  fPhiH2        = bases->GetPhiH2();
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);
  for (Int_t i=0; i<kNparton; i++) pv[i].Boost(boostv);

  // ----------------------------------------------
  //  Set final state parton infomation
  // ----------------------------------------------

  fNparton = kNparton;

  static TVector **qp = 0;
  if (!qp) {
    qp = new TVector* [kNparton];
    for (Int_t i=0; i<kNparton; i++) qp[i] = new TVector(4);
  }
  for (Int_t i=0; i<kNparton; i++) {
    TVector &q = *qp[i];
    q(0) = pv[i].E ();
    q(1) = pv[i].Px();
    q(2) = pv[i].Py();
    q(3) = pv[i].Pz();
  }
  Int_t    idx     = bases->fXBosonPtr ->GetPID(); // PDG code for ZH? 
  Double_t mass    = bases->GetMass();

  Int_t    iddm    = bases->fDMBosonPtr->GetPID(); // PDG code for AH? 
  Double_t msdm    = bases->GetMassDM();
  Double_t mh      = bases->GetMassHiggs();

  Int_t    idh     = 25;                           // PDG code for H

  //                               No. PID   Mass Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton( 1, idx , mass,    0., *qp[0], 2, 3,  0, 0,   0,     0);
  new (partons[1]) JSFSpringParton( 2, idx , mass,    0., *qp[1], 2, 5,  0, 0,   0,     0);
  new (partons[2]) JSFSpringParton( 3, idh , mh  ,    0., *qp[2], 0, 0,  1, 0,   0,     0);
  new (partons[3]) JSFSpringParton( 4, iddm, msdm,    0., *qp[3], 0, 0,  1, 0,   0,     0);
  new (partons[4]) JSFSpringParton( 5, idh , mh  ,    0., *qp[4], 0, 0,  2, 0,   0,     0);
  new (partons[5]) JSFSpringParton( 6, iddm, msdm,    0., *qp[5], 0, 0,  2, 0,   0,     0);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHZHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHZHBases::ZHZHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
	   fF         ( 580.),
	   fKappaL    (  0.5),
           fMass      ( 368.),
           fMassDM    ( 81.9),
           fMassT     ( 410.),
           fMassHiggs ( 134.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fPole      (0.),
           fXBosonPtr ( 0),
           fDMBosonPtr( 0),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fPhotonPtr ( 0),
           fGluonPtr  ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2XX      (0.),
           fQ2X1      (0.),
           fQ2X2      (0.),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaH1(0.),
           fPhiH1     (0.),
           fXQ2X1     (0.),
           fCosThetaH2(0.),
           fPhiH2     (0.),
           fXQ2X2     (0.),
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

  cout << "Init ZHZHBases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("ZHZHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.CosthH1Range","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.PhiH1OverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.CosthH2Range","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.PhiH2OverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.F","580.")); 	 // F [GeV]
  ins >> fF;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.KappaL","0.5")); 	 // kappa_l
  ins >> fKappaL;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.MassHiggs","134."));  // M_h [GeV]
  ins >> fMassHiggs;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHZHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fXQ2X1         , 0., 1., 0, 1);
  DefineVariable(fXQ2X2         , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();

  DefineVariable(fCosTheta  , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi       , fXL[1], fXU[1], 0, 0);
  DefineVariable(fCosThetaH1, fXL[2], fXU[2], 0, 0);
  DefineVariable(fPhiH1     , fXL[3], fXU[3], 0, 0);
  DefineVariable(fCosThetaH2, fXL[4], fXU[4], 0, 0);
  DefineVariable(fPhiH2     , fXL[5], fXU[5], 0, 0);

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
  SetNoOfSample(10000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 20);
  SetIteration2(0.05,100);

}
// --------------------------
//  D-tor
// --------------------------
ZHZHBases::~ZHZHBases()
{
  delete fXBosonPtr;
  delete fDMBosonPtr;
  delete fWBosonPtr;
  delete fZBosonPtr;
  delete fPhotonPtr;
  delete fGluonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZHZHBases::Func()
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

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < fMassHiggs + 2*fMassDM) {
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
  fQ2XX = fEcmIP*fEcmIP;

  Double_t rs   = fEcmIP;
  Double_t qmin = fMassHiggs + fMassDM;
  Double_t qmax = rs - fMassHiggs - fMassDM;
#ifndef __ZEROWIDTH__
  fQ2X1 = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2X1, weight);
  Double_t qx1 = TMath::Sqrt(fQ2X1);
  bsWeight *= weight;
  qmax = rs - qx1;
  fQ2X2 = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2X2, weight);
  bsWeight *= weight;
#else
  fQ2X1 = TMath::Power(fXBosonPtr->GetMass(),2);
  fQ2X2 = TMath::Power(fXBosonPtr->GetMass(),2);
  weight = TMath::Power(kPi*fXBosonPtr->GetMass()*fXBosonPtr->GetWidth(),2);
#ifdef __NODECAY__
  weight = 1.;
#endif
  bsWeight *= weight;
#endif

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch x1branch(fQ2X1, fCosThetaH1, fPhiH1, fMassHiggs*fMassHiggs, fMassDM*fMassDM);
  GENBranch x2branch(fQ2X2, fCosThetaH2, fPhiH2, fMassHiggs*fMassHiggs, fMassDM*fMassDM);
  GENBranch cmbranch(fQ2XX, fCosTheta  , fPhi  , &x1branch            , &x2branch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X1), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaH1       , (bsWeight*sigma));
  Xh_fill( 6, fPhiH1            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2X2), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaH2       , (bsWeight*sigma));
  Xh_fill( 9, fPhiH2            , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZHZHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  // CM -> X X
  Double_t q2x1  = cmbranch.GetM12();
  Double_t q2x2  = cmbranch.GetM22();
  Double_t cosx  = cmbranch.GetCosTheta();
  Double_t phix  = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x1, q2x2, cmframe, cosx, phix, 0);
  ANL4DVector px1 = phaseCM.GetFourMomentum(0);
  ANL4DVector px2 = phaseCM.GetFourMomentum(1);
  Double_t betax  = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  // X1 -> H1 DM1
  GENBranch &x1branch = * cmbranch.GetBranchPtr(0);
  Double_t cosh1 = x1branch.GetCosTheta();
  Double_t phih1 = x1branch.GetPhi     ();
  Double_t mh21  = x1branch.GetM12();
  Double_t mdm21 = x1branch.GetM22();
  GENPhase2 phaseX1(px1, mh21, mdm21, cmframe, cosh1, phih1, 1);
  fP[0]          = phaseX1.GetFourMomentum(0);
  fM[0]          = TMath::Sqrt(mh21);
  fP[1]          = phaseX1.GetFourMomentum(1);
  fM[1]          = TMath::Sqrt(mdm21);
  Double_t betah1 = phaseX1.GetBetaBar();
  if (betah1 <= 0.) return 0.;

  // X2 -> H2 DM2
  GENBranch &x2branch = * cmbranch.GetBranchPtr(1);
  Double_t cosh2 = x2branch.GetCosTheta();
  Double_t phih2 = x2branch.GetPhi     ();
  Double_t mh22  = x2branch.GetM12();
  Double_t mdm22 = x2branch.GetM22();
  GENPhase2 phaseX2(px2, mh22, mdm22, cmframe, cosh2, phih2, 1);
  fP[2]          = phaseX2.GetFourMomentum(0);
  fM[2]          = TMath::Sqrt(mh22);
  fP[3]          = phaseX2.GetFourMomentum(1);
  fM[3]          = TMath::Sqrt(mdm22);
  Double_t betah2 = phaseX2.GetBetaBar();
  if (betah2 <= 0.) return 0.;

  Double_t eb     = rs/2.;
  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);
#ifdef __DEBUG__
  cerr << "---------------------" << endl;
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
  ANL4DVector qx1 = fP[0] + fP[1];
  ANL4DVector qx2 = fP[2] + fP[3];
  cerr << " ---- " << endl;
  cerr << " rs = " << rs << endl;
  cerr << " mh1 = " << fP[0].GetMass() << endl;
  cerr << " md1 = " << fP[1].GetMass() << endl;
  cerr << " mh2 = " << fP[3].GetMass() << endl;
  cerr << " md2 = " << fP[4].GetMass() << endl;
  cerr << " mx1 = " << qx1.GetMass() << endl;
  cerr << " mx2 = " << qx2.GetMass() << endl;

  ANL4DVector pcm = qx1 + qx2;
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << qcm.Py() << ","
       << qcm.Pz() << ")" << endl;
#endif

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
  Double_t s      = rs * rs;

  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared(cmbranch);

  // -------------------
  //  Put them together
  // -------------------
#ifndef __NODECAY__
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
#else
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3+2));
  betah1 = 1.;
  betah2 = 1.;
#endif

  Double_t identp = 1./2.;                             // identical particle factor
  Double_t dPhase = kFact * betax * betah1 * betah2;   // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);               // beam flux factor
  Double_t spin   = 1./2.;                             // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ZHZHBases::AmpSquared(GENBranch &cmbranch)
{
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2);

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ZHZHBases::FullAmplitude()
{
   Double_t mx     = fXBosonPtr ->GetMass();
   Double_t gamx   = fXBosonPtr ->GetWidth();
   Double_t mdm    = fDMBosonPtr->GetMass();

   Double_t Ch     = fXBosonPtr->GetCh();
   Double_t Sh     = fXBosonPtr->GetSh();
   Double_t Cf     = fXBosonPtr->GetCf();
   Double_t Sf     = fXBosonPtr->GetSf();
   Double_t gxhdm  = -kGw*kGw*Sf*Cf*fF/(kCos2W*2.*TMath::Sqrt(2.))
	                  *(Ch*kCosW+Sh*kSinW)*(Sh*kCosW-Ch*kSinW);
#if 0
   cerr << " Ch = " << Ch << " Sh = " << Sh << " Cf = " << Cf << " Sf = " << Sf << endl;
   cerr << " gxhdm = " << gxhdm << " kGw = " << kGw << endl;
#endif

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

#ifndef __NODECAY__
   HELScalar  h1 (fP[0], +1);                       // H1  from ZH1
   HELVector  dm1(fP[1], mdm  , fHelFinal[1], +1);  // AH1 from ZH1
   HELVector  x1(dm1, h1, gxhdm, mx, gamx);         // ZH1

   HELScalar  h2 (fP[2], +1);                       // H2  from ZH2
   HELVector  dm2(fP[3], mdm  , fHelFinal[3], +1);  // AH2 from ZH2
   HELVector  x2(dm2, h2, gxhdm, mx, gamx);         // ZH2
#else
   ANL4DVector px1 = fP[0] + fP[1];
   ANL4DVector px2 = fP[2] + fP[3];
   HELVector   x1(px1 , mx , fHelFinal[1], +1);
   HELVector   x2(px2 , mx , fHelFinal[3], +1);
#endif
   Complex_t amp = AmpEEtoXX(em, ep, x1, x2);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoXX()
// --------------------------
Complex_t ZHZHBases::AmpEEtoXX(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &x1,
                               const HELVector  &x2)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t Ch     = fXBosonPtr->GetCh();
   Double_t Sh     = fXBosonPtr->GetSh();
   Double_t  glehzh = kGw*(Ch-Sh*kSinW/(5.*kCosW))/2.;
   Double_t  grehzh = 0.;

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
   // XX Production Amplitude
   //---------------------------
   //--
   // T-channel eH
   //--
   Double_t dummy = 5.; // dummy width for eH
   HELFermion eh1(em, x1, glehzh, grehzh, fMassT, dummy);
   HELVertex  amptee1(eh1,ep,x2,glehzh,grehzh);

   HELFermion eh2(em, x2, glehzh, grehzh, fMassT, dummy);
   HELVertex  amptee2(eh2,ep,x1,glehzh,grehzh);

   Complex_t amp = amptee1 + amptee2;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZHZHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZHZHBases.BeamstrahlungFilepath",
                                           "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZHZHBases.BeamstrahlungFilename","trc500"));
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
    fBM->SetIBParameters(0.0);

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
  fWBosonPtr->DebugPrint();
  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  fZBosonPtr->DebugPrint();
  if (!fPhotonPtr) fPhotonPtr = new GENPDTPhoton();
  //fPhotonPtr->DebugPrint();
  if (!fGluonPtr)  fGluonPtr  = new GENPDTGluon();
  //fGluonPtr->DebugPrint();

  cerr << " ZH Boson "      << endl;
  if (!fXBosonPtr)  fXBosonPtr  = new ZHBoson(fF, fMassHiggs);
  fXBosonPtr->DebugPrint();
  fMass   = fXBosonPtr->GetMass();
  fMassDM = fXBosonPtr->GetMassDM();

  fMassT  = TMath::Sqrt(2.)*fKappaL*fF;
  cerr << " M_ZH = " << fMass   << " [GeV]" << endl
       << " M_AH = " << fMassDM << " [GeV]" << endl
       << " M_eH = " << fMassT  << " [GeV]" << endl;

  cerr << " AH Boson " << endl;
  if (!fDMBosonPtr) fDMBosonPtr = new AHBoson(fMassDM);
  fDMBosonPtr->DebugPrint();

  Double_t mx = fMass;

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4, mx-10., mx+10.,       50, "Mx1"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthH1");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiH1"  );
  Xh_init( 7, mx-10., mx+10.,       50, "Mx2"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthH2");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiH2"  );
  Xh_init(10,     0.,     2.,        2, "Helin " );
  Xh_init(11,     0.,     9.,        9, "Helot " );
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZHZHBases::Userout()
{
  cout << "End of ZHZHBases----------------------------------- "  << endl
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
void ZHZHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 9;
   //                                       h  dm  h  dm  
   static const Int_t kFHelComb[kNf][4] = {{0, -1, 0, -1},
                                           {0, -1, 0,  0},
                                           {0, -1, 0, +1},
                                           {0,  0, 0, -1},
                                           {0,  0, 0,  0},
                                           {0,  0, 0, +1},
                                           {0, +1, 0, -1},
                                           {0, +1, 0,  0},
                                           {0, +1, 0, +1}};
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

//-----------------------------------------------------------------------------
// ==============================
//  class ZHBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHBoson::ZHBoson(Double_t f,
                 Double_t mh)
	: fF(f), fMassHiggs(mh)
{
   fCf          = 1.-4.*kCos2W*kSin2W*kM_z*kM_z/TMath::Power(kGe*f,2);
   fSf          = TMath::Sqrt((1.-fCf)*(1.+fCf));
   Double_t g2  = kGw;
   Double_t gp  = kGe/kCosW;
   Double_t maa = g2*g2*f*f*(fCf*fCf+7.)/8.;
   Double_t mab = g2*gp*f*f*(1.-fCf)*(1.+fCf)/8.;
   Double_t mbb = gp*gp*f*f*(5.*fCf*fCf+3.)/40.;

   fName    = TString("ZH");
   fPID     = 200000003;
   fCharge  =  0.0;
   fSpin    =  1.0;
   fMass    = TMath::Sqrt(0.5*(maa+mbb+TMath::Sqrt(TMath::Power(maa-mbb,2)+4*mab*mab))); 
   fGen     =    0;
   fIsoSpin =  0.0;
   fColor   =  1.0;

   fSh      = -2.*mab/TMath::Sqrt(
                  TMath::Power(maa-mbb+TMath::Sqrt(TMath::Power(maa-mbb,2)+4*mab*mab),2)
                  +4.*mab*mab);
   fCh      = TMath::Sqrt((1.-fSh)*(1.+fSh));
   fMassDM  = TMath::Sqrt(0.5*(maa+mbb-TMath::Sqrt(TMath::Power(maa-mbb,2)+4.*mab*mab)));

   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void ZHBoson::Initialize()
{
   //--
   // Couplings
   //--
   Double_t  a  = -kGw*kGw*fSf*fCf*fF/(kCos2W*2.*TMath::Sqrt(2.))
	                  *(fCh*kCosW+fSh*kSinW)*(fSh*kCosW-fCh*kSinW);
   //--
   // ZH --> H + AH
   //--
   GENDecayMode *dmp;
   GENPDTEntry  *d1p, *d2p;
   Int_t    idh   = 25;
   Double_t chgh  = 0.0;
   Double_t spinh = 0.0;
   d1p = new GENPDTEntry("H",idh,chgh,spinh,fMassHiggs);
   d2p = new AHBoson(fMassDM);
   Double_t m1    = d1p->GetMass();
   Double_t m2    = d2p->GetMass();
   Double_t ident = 1.;
   Double_t gam = GamToSV(m1, m2, a)/ident;
   if (gam > 0.) {
      dmp = new GENDecayMode(gam);
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
}

//_____________________________________________________________________________
// --------------------------
//  GamToSV
// --------------------------
Double_t ZHBoson::GamToSV(Double_t m1, // 1st daughter mass
                          Double_t m2, // 2nd daughter mass
                          Double_t a)  // coupling
{
   Double_t x1   = TMath::Power(m1/fMass,2);
   Double_t x2   = TMath::Power(m2/fMass,2);
   Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

   if (beta <= 0.) return 0.;
   beta = TMath::Sqrt(beta);

   ANL4DVector px(fMass,0.,0.,0.);
   double p = fMass*beta/2;
   ANL4DVector ph (sqrt(p*p+m1*m1),0.,0., p);
   ANL4DVector pdm(sqrt(p*p+m2*m2),0.,0.,-p);
   HELScalar h(ph, +1);
   HELVector *xPtr[3], *dPtr[3];
   for (Int_t i=0; i<3; i++) {
     xPtr[i] = new HELVector(px , fMass, i-1, +1);
     dPtr[i] = new HELVector(pdm, m2   , i-1, +1);
   } 
   Double_t amp2 = 0.;
   for (Int_t jx=0; jx<3; jx++) {
      for (Int_t jd=0; jd<3; jd++) {
         HELVertex amp(*xPtr[jx], *dPtr[jd], h, a);
         amp2 += pow(abs(amp),2);
      }
   }
   Double_t spin = 2*fSpin+1;
   Double_t fac  = 1./(16.*kPi)/fMass/spin;
   Double_t gam  = fac*amp2*beta;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class AHBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
AHBoson::AHBoson(Double_t m)
{
   fName    = TString("AH");
#if 0
   fPID     = 200000001;
#else
   fPID     = 220000; // LSP code for JSFHadronizer.
#endif
   fCharge  =  0.0;
   fSpin    =  1.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin =  0.0;
   fColor   =  1.0;
}
