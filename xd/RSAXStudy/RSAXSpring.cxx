//*****************************************************************************
//* =====================
//*  RSAXSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- --> AX generator
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "RSAXSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(RSAXSpring)
ClassImp(RSAXSpringBuf)
ClassImp(RSAXBases)
ClassImp(RSXBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class RSAXSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSAXSpring::RSAXSpring(const char      *name,
                       const char      *title,
                             RSAXBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new RSAXSpringBuf("RSAXSpringBuf",
                                "RSAXSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new RSAXBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
RSAXSpring::~RSAXSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class RSAXSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t RSAXSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  RSAXBases    *bases   = (RSAXBases*)((RSAXSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 4;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];
  pv[3] = bases->fP[1];
  pv[1] = bases->fP[2];
  pv[0] = pv[2] + pv[3];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2AX         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X          = bases->GetQ2X();
  fCosThetaA    = bases->GetCosThetaA();
  fPhiA         = bases->GetPhiA();
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
  Int_t    idx     = bases->fXBosonPtr->GetPID(); // PDG code for X??
  Int_t    ida     = 22;                          // PDG code for gamma
  Int_t    id1     = bases->f1Ptr->GetPID();      // PDG code for 1st daughter
  Int_t    id2     = bases->f2Ptr->GetPID();      // PDG code for 1st daughter
  Double_t m1      = bases->f1Ptr->GetMass();     // 1st daughter mass
  Double_t m2      = bases->f2Ptr->GetMass();     // 1st daughter mass
  Double_t color   = bases->f1Ptr->GetColor();    // color factor
  Double_t chg1    = bases->f1Ptr->GetCharge();   // 1st daughter charge
  Double_t chg2    = bases->f2Ptr->GetCharge();   // 2nd daughter charge
  Int_t    islev   = color > 1. ? 101 : 0;        // shower level
  Int_t    icf     = 1;                           // color flux id

  Double_t mass    = bases->GetMass();
#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idx << " m=" << mass  << " Q=" << 0.
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 2 pid=" << idz << " m=" << rq2z << " Q=" << 0.
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 3 pid=" << ida << " m=" << 0. << " Q=" << 0.
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 4 pid=" << ida << " m=" << 0. << " Q=" << 0.
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  cerr << " 5 pid=" << idf << " m=" << m3 << " Q=" << chrg
       << " pv=(" << (*qp[4])(0) << ", "
                  << (*qp[4])(1) << ", "
                  << (*qp[4])(2) << ", "
                  << (*qp[4])(3) << ") " << endl;
  TVector qcm(4), qx(4), qa(4);
  qx  = *qp[2] + *qp[3];
  qa  = *qp[1];
  qcm = qx +qa;
  cerr << " px=(" << qx[0] << ", "
                  << qx[1] << ", "
                  << qx[2] << ", "
                  << qx[3] << ") " << endl;
  cerr << " pa=(" << qa[0] << ", "
                  << qa[1] << ", "
                  << qa[2] << ", "
                  << qa[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                              No. PID  Mass  Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton(1, idx, mass,    0., *qp[0], 2, 3,  0, 0,   0,     0);
  new (partons[1]) JSFSpringParton(2, ida,   0.,    0., *qp[1], 0, 0,  0, 0,   0,     0);
  new (partons[2]) JSFSpringParton(3, id1,   m1,  chg1, *qp[2], 0, 0,  1, 0, icf, islev);
  new (partons[3]) JSFSpringParton(4, id2,   m2,  chg2, *qp[3], 0, 0,  1, 0, icf, islev);
  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class RSAXBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSAXBases::RSAXBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fLambda    (1000.),
           fC0        (0.),
           fC1        (1.),
           fC2        (1.),
           fC3        (1.),
           fC4        (1.),
           fMass      ( 120.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fPole      (0.),
           fXBosonPtr ( 0),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fPhotonPtr ( 0),
           fGluonPtr  ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2AX      (0.),
           fQ2X       (0.),
           fXModePtr  (0),
           f1Ptr      (0),
           f2Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaA (0.),
           fPhiA      (0.),
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

  cout << "Init rszxbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("RSAXBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.CosthARange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.PhiAOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Lambda","1000."));    // Lambda [GeV]
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C0","0.")); 	  	 // C_0
  ins >> fC0;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C1","1.")); 		 // C_1
  ins >> fC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C2","1.")); 		 // C_2
  ins >> fC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C3","1.")); 		 // C_3
  ins >> fC3;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C4","1.")); 		 // C_4
  ins >> fC4;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.MassX","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Bremstrahlung","1")); // ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSAXBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSAXBases.BeamstrahlungFilename","trc500"));
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
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);

  DefineVariable(fXDecayMode    , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaA, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiA     , fXL[3], fXU[3], 0, 1);

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
  SetNoOfSample(20000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 20);
  SetIteration2(0.05, 50);

}
// --------------------------
//  D-tor
// --------------------------
RSAXBases::~RSAXBases()
{
  delete fXBosonPtr;
  delete fWBosonPtr;
  delete fZBosonPtr;
  delete fPhotonPtr;
  delete fGluonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t RSAXBases::Func()
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
  } // if Bremstrahlung


  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung

  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  Double_t weight = 1;
  fXModePtr = fXBosonPtr->PickMode(fXDecayMode, weight, fXMode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fXModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fXModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();
  Double_t m3   = 0.;

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < fMass + m3) {
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
  fQ2AX = fEcmIP*fEcmIP;
  fQ2X  = fMass*fMass; // narrow width approx. for X

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch xbranch (fQ2X , fCosThetaA, fPhiA, m1, m2);
  GENBranch cmbranch(fQ2AX, fCosTheta , fPhi , &xbranch, m3);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP           , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta        , (bsWeight*sigma));
  Xh_fill( 3, fPhi             , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaA       , (bsWeight*sigma));
  Xh_fill( 6, fPhiA            , (bsWeight*sigma));
  Xh_fill( 7, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill( 8, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill( 9, (Double_t)fXMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t RSAXBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2x  = cmbranch.GetM12();
  Double_t q2z  = cmbranch.GetM22();
  Double_t cosx = cmbranch.GetCosTheta();
  Double_t phix = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x, q2z, cmframe, cosx, phix, 0);
  ANL4DVector px = phaseCM.GetFourMomentum(0);
           fP[2] = phaseCM.GetFourMomentum(1);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  GENBranch &xbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosa = xbranch.GetCosTheta();
  Double_t phia = xbranch.GetPhi     ();
  Double_t m12  = xbranch.GetM12();
  Double_t m22  = xbranch.GetM22();
  GENPhase2 phaseX(px, m12, m22, cmframe, cosa, phia, 1);
  fP[0] = phaseX.GetFourMomentum(0);
  fP[1] = phaseX.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betaa = phaseX.GetBetaBar();
  if (betaa <= 0.) return 0.;

  Double_t eb     = rs/2.;
  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

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
  static const Int_t    kNbr  = 2;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax * betaa;         // phase space factor
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
Double_t RSAXBases::AmpSquared(GENBranch &cmbranch)
{
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2);
#if 1
  amp2 *= k2Pi * k8Pi; // branch factor correction: no phase space for X
  amp2 *= fXModePtr->GetBR();
#endif

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t RSAXBases::FullAmplitude()
{
   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  xf(fP[0]+fP[1], kFALSE);
   HELVector  af(fP[2], 0., fHelFinal[2], +1);

   Complex_t amp = AmpEEtoAX(em, ep, xf, af);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoAX()
// --------------------------
Complex_t RSAXBases::AmpEEtoAX(const HELFermion &em,
                               const HELFermion &ep,
                               const HELScalar  &xf,
                               const HELVector  &af)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t  a1    = fC1*kCos2W + fC2*kSin2W;
   Double_t  a2    = (fC1 - fC2) * kSinCosW;
   Double_t  a3    = fC1*kSin2W + fC2*kCos2W;
   Double_t  gaax  = a1*kSqh/fLambda;
   Double_t  gazx  = a2*kSqh/fLambda;
   Double_t  gzzx  = a3*kSqh/fLambda;

   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

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
   // AX Production Amplitude
   //---------------------------
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel   Z: backward-going
   ANL4DVector k1   = zs.GetFourMomentum();      // s-channel   Z: backward-going
   ANL4DVector k2   = af.GetFourMomentum();      // final state A: forward-going 

   Double_t    k1k2 = k1*k2;
   Complex_t   k1af = k1(0)*af[0] - k1(1)*af[1] - k1(2)*af[2] - k1(3)*af[3];
   Complex_t   zsaf = zs[0]*af[0] - zs[1]*af[1] - zs[2]*af[2] - zs[3]*af[3];
   Complex_t   k2zs = zs[0]*k2(0) - zs[1]*k2(1) - zs[2]*k2(2) - zs[3]*k2(3);
   Complex_t   ampzax = gazx * (k1k2*zsaf - k1af*k2zs);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma

   Complex_t   asaf = as[0]*af[0] - as[1]*af[1] - as[2]*af[2] - as[3]*af[3];
   Complex_t   k2as = as[0]*k2(0) - as[1]*k2(1) - as[2]*k2(2) - as[3]*k2(3);
   Complex_t   ampaax = gaax * (k1k2*asaf - k1af*k2as);

   //--
   // Sum of the two
   //--
   Complex_t amp = ampzax + ampaax;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void RSAXBases::Userin()
{
  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
   if (!fXBosonPtr) fXBosonPtr = new RSXBoson(fMass, 
                                              fLambda,
                                              fC0,
                                              fC1,
                                              fC2,
                                              fC3,
                                              fC4);
   fXBosonPtr->DebugPrint();
   if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
   fWBosonPtr->DebugPrint();
   if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
   fZBosonPtr->DebugPrint();
   if (!fPhotonPtr) fPhotonPtr = new GENPDTPhoton();
   //fPhotonPtr->DebugPrint();
   if (!fGluonPtr)  fGluonPtr  = new GENPDTGluon();
   //fGluonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh" );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init( 4,   100.,   140.,       50, "Mx"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthA");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiA"  );
  Xh_init( 7,     0.,     2.,        2, "Helin ");
  Xh_init( 8,     0.,     2.,        2, "Helot ");
  Xh_init( 9,     0.,     2.,        2, "X mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void RSAXBases::Userout()
{
  cout << "End of RSAXBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremstrahlung        = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void RSAXBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][3] = {{0, 0, +1},
                                           {0, 0, -1}};
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
   weight = kNf;
}

//-----------------------------------------------------------------------------
// ==============================
//  class RSXBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSXBoson::RSXBoson(Double_t m,
                   Double_t lambda,
                   Double_t c0,
                   Double_t c1,
                   Double_t c2,
                   Double_t c3,
                   Double_t c4)
        : fLambda(lambda),
          fC0(c0),
          fC1(c1),
          fC2(c2),
          fC3(c3),
          fC4(c4)
{
   fName    = TString("X");
   fPID     = 200000000;
   fCharge  = 0.;
   fSpin    = 0.;
   fMass    = m;
   fGen     = 0;
   fIsoSpin = 0.;
   fColor   = 1.;

   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void RSXBoson::Initialize()
{
   //--
   // Couplings
   //--
   Double_t  a1    = fC1*kCos2W + fC2*kSin2W;
   Double_t  a2    = (fC1 - fC2) * kSinCosW;
   Double_t  a3    = fC1*kSin2W + fC2*kCos2W;
   //--
   // X --> gluon gluon
   //--
   GENDecayMode *dmp;
   GENPDTEntry  *d1p, *d2p;
   d1p = new GENPDTGluon();
   d2p = new GENPDTGluon();
   Double_t cf    = 8.;
   Double_t m1    = 0.;
   Double_t m2    = 0.;
   Double_t ident = 2.;
   Double_t gam = GamToVV(m1, m2, fC4, cf)/ident;
   if (gam > 0.) {
      dmp = new GENDecayMode(gam);
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
   //--
   // X --> gamma gamma 
   //--
   d1p = new GENPDTPhoton();
   d2p = new GENPDTPhoton();
   cf    = 1.;
   m1    = 0.;
   m2    = 0.;
   ident = 2.;
   gam = GamToVV(m1, m2, a1, cf)/ident;
   if (gam > 0.) {
      dmp = new GENDecayMode(gam);
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
   //--
   // X --> Z gamma 
   //--
   d1p = new GENPDTZBoson();
   d2p = new GENPDTPhoton();
   cf    = 1.;
   m1    = d1p->GetMass();
   m2    = d2p->GetMass();
   ident = 1.;
   gam = GamToVV(m1, m2, a2, cf)/ident;
   if (gam > 0.) {
      dmp = new GENDecayMode(gam);
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
}

//_____________________________________________________________________________
// --------------------------
//  GamToVV
// --------------------------
Double_t RSXBoson::GamToVV(Double_t m1, // 1st daughter mass
                           Double_t m2, // 2nd daughter mass
                           Double_t a,  // coupling
                           Double_t cf) // color factor
{
   Double_t x1   = TMath::Power(m1/fMass,2);
   Double_t x2   = TMath::Power(m2/fMass,2);
   Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

   if (beta <= 0.) return 0.;
   beta = TMath::Sqrt(beta);

   Double_t tta = a*a*fMass*TMath::Power(fMass/fLambda,2)
	          * (3. + 2.*TMath::Power(beta,2) + 3.*TMath::Power(beta,4))/8.;
   Double_t fac = 1./(64.*kPi);
   Double_t gam = fac*tta*beta*cf;

   return gam;
}
