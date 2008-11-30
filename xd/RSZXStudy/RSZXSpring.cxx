//*****************************************************************************
//* =====================
//*  RSZXSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- --> ZX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "RSZXSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __HIGGS__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(RSZXSpring)
ClassImp(RSZXSpringBuf)
ClassImp(RSZXBases)
ClassImp(RSXBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class RSZXSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSZXSpring::RSZXSpring(const char      *name,
                       const char      *title,
                             RSZXBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new RSZXSpringBuf("RSZXSpringBuf",
                                "RSZXSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new RSZXBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
RSZXSpring::~RSZXSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class RSZXSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t RSZXSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  RSZXBases    *bases   = (RSZXBases*)((RSZXSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 6;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];
  pv[3] = bases->fP[1];
  pv[4] = bases->fP[2];
  pv[5] = bases->fP[3];
  pv[0] = pv[2] + pv[3];
  pv[1] = pv[4] + pv[5];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2ZX         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X          = bases->GetQ2X();
  fCosThetaA    = bases->GetCosThetaA();
  fPhiA         = bases->GetPhiA();
  fQ2Z          = bases->GetQ2Z();
  fCosThetaF    = bases->GetCosThetaF();
  fPhiF         = bases->GetPhiF();
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
  Int_t    id1     = bases->f1Ptr->GetPID();      // PDG code for 1st daughter
  Int_t    id2     = bases->f2Ptr->GetPID();      // PDG code for 1st daughter
  Double_t m1      = bases->f1Ptr->GetMass();     // 1st daughter mass
  Double_t m2      = bases->f2Ptr->GetMass();     // 1st daughter mass
  Double_t colorx  = bases->f1Ptr->GetColor();    // color factor
  Double_t chg1    = bases->f1Ptr->GetCharge();   // 1st daughter charge
  Double_t chg2    = bases->f2Ptr->GetCharge();   // 2nd daughter charge
  Int_t    islevx  = colorx > 1. ? 101 : 0;       // shower level
  Int_t    icfx    = 1;                           // color flux id

  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    idf     = bases->f3Ptr->GetPID   ();   // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge();   // F charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // F mass
  Double_t m4      = m3;                          // F mass
  Int_t    hel3    = bases->fHelFinal[2];         // f1 helicity
  Int_t    hel4    = bases->fHelFinal[3];         // f2 helicity
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	  // shower level
  Int_t    icf     = 2;                           // color flux id
  Double_t rq2z    = pv[1].Mag();

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
  cerr << " 6 pid=" <<-idf << " m=" << m4 << " Q=" << -chrg
       << " pv=(" << (*qp[5])(0) << ", "
                  << (*qp[5])(1) << ", "
                  << (*qp[5])(2) << ", "
                  << (*qp[5])(3) << ") " << endl;
  TVector qcm(4), qx(4), qz(4);
  qx  = *qp[2] + *qp[3];
  qz  = *qp[4] + *qp[5];
  qcm = qx +qz;
  cerr << " px=(" << qx[0] << ", "
                  << qx[1] << ", "
                  << qx[2] << ", "
                  << qx[3] << ") " << endl;
  cerr << " pz=(" << qz[0] << ", "
                  << qz[1] << ", "
                  << qz[2] << ", "
                  << qz[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                              No. PID  Mass  Charge   pv   Nd 1st Mom  hel  col shower
  new (partons[0]) JSFSpringParton(1, idx, mass,    0., *qp[0], 2, 3,  0,    0,   0,     0);
  new (partons[1]) JSFSpringParton(2, idz, rq2z,    0., *qp[1], 2, 5,  0,    0,   0,     0);
  new (partons[2]) JSFSpringParton(3, id1,   m1,  chg1, *qp[2], 0, 0,  1,    0,icfx,islevx);
  new (partons[3]) JSFSpringParton(4, id2,   m2,  chg2, *qp[3], 0, 0,  1,    0,icfx,islevx);
  new (partons[4]) JSFSpringParton(5, idf,   m3,  chrg, *qp[4], 0, 0,  2, hel3, icf, islev);
  new (partons[5]) JSFSpringParton(6,-idf,   m4, -chrg, *qp[5], 0, 0,  2, hel4, icf, islev);
  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class RSZXBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSZXBases::RSZXBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fLambda    (1000.),
           fC0        (0.),
           fC1        (1.),
           fC2        (1.),
           fC3        (1.),
           fMass      ( 120.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fPole      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fXBosonPtr ( 0),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fPhotonPtr ( 0),
           fGluonPtr  ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2ZX      (0.),
           fQ2X       (0.),
           fXModePtr  (0),
           f1Ptr      (0),
           f2Ptr      (0),
           fQ2Z       (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaA (0.),
           fPhiA      (0.),
           fXQ2Z      (0.),
           fCosThetaF (0.),
           fPhiF      (0.),
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
  stringstream ins(gJSF->Env()->GetValue("RSZXBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.CosthARange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiAOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Lambda","1000."));    // Lambda [GeV]
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C0","0.")); 	  	 // C_0
  ins >> fC0;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C1","1.")); 		 // C_1
  ins >> fC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C2","1.")); 		 // C_2
  ins >> fC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C3","1.")); 		 // C_3
  ins >> fC3;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.MassX","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fXDecayMode    , 0., 1., 0, 1);
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
  DefineVariable(fCosThetaA, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiA     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaF, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiF     , fXL[5], fXU[5], 0, 1);

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
RSZXBases::~RSZXBases()
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
Double_t RSZXBases::Func()
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

  fXModePtr = fXBosonPtr->PickMode(fXDecayMode, weight, fXMode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fXModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fXModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < fMass + m3 + m4) {
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
  fQ2ZX = fEcmIP*fEcmIP;
  fQ2X  = fMass*fMass; // narrow width approx. for X

  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - fMass;
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch xbranch (fQ2X , fCosThetaA, fPhiA, m1*m1, m2*m2);
  GENBranch zbranch (fQ2Z , fCosThetaF, fPhiF, m3*m3, m4*m4);
  GENBranch cmbranch(fQ2ZX, fCosTheta , fPhi , &xbranch, &zbranch);

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
  Xh_fill( 7, TMath::Sqrt(fQ2Z), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaF       , (bsWeight*sigma));
  Xh_fill( 9, fPhiF            , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill(12, (Double_t)fXMode , (bsWeight*sigma));
  Xh_fill(13, (Double_t)fZMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t RSZXBases::DSigmaDX(GENBranch &cmbranch)
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
  ANL4DVector pz = phaseCM.GetFourMomentum(1);
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

  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared(cmbranch);

  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax * betaa * betaf; // phase space factor
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
Double_t RSZXBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;
#if 1
  amp2 *= k2Pi * k8Pi; // branch factor correction: no phase space for X
  amp2 *= fXModePtr->GetBR();
#endif

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
Complex_t RSZXBases::FullAmplitude()
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

   HELScalar  xf(fP[0]+fP[1]);

   HELFermion f (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
   HELFermion fb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoZX(em, ep, xf, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZX()
// --------------------------
Complex_t RSZXBases::AmpEEtoZX(const HELFermion &em,
                               const HELFermion &ep,
                               const HELScalar  &xf,
                               const HELVector  &zf)
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
#if 1
   static Int_t ncall = 0;
   if (!ncall) {
      ncall = 1;
      cerr << " a1 = " << a1
           << " a2 = " << a2
           << " a3 = " << a3 << endl;
      cerr << " kSqh = "    << kSqh 
           << " fLambda = " << fLambda << endl;
   }
#endif

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
#ifdef __HIGGS__
   //---------------------------
   // Higgs Production Amplitude
   //---------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Double_t gzzh  = kGz*kM_z;
   Complex_t amp = HELVertex(zs, zf, xf, gzzh);
#else
   //---------------------------
   // ZX Production Amplitude
   //---------------------------
#if 0
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   ANL4DVector k1b  = zs.GetFourMomentum();      // s-channel   Z: backward-going
   ANL4DVector k1   = ANL4DVector(-k1b.E(),      // s-channel   Z: forward-going
                                  -k1b.Px(), 
                                  -k1b.Py(),
                                  -k1b.Pz());
   ANL4DVector k2   = zf.GetFourMomentum();      // final state Z: forward-going 

   Double_t    k1k2 = k1*k2;
   Complex_t   k1zf = k1(0)*zf[0] - k1(1)*zf[1] - k1(2)*zf[2] - k1(3)*zf[3];
   Complex_t   zszf = conj(zs[0])*zf[0] - conj(zs[1])*zf[1]
                                        - conj(zs[2])*zf[2]
                                        - conj(zs[3])*zf[3];
   Complex_t   k2zs = conj(zs[0])*k2(0) - conj(zs[1])*k2(1)
                                        - conj(zs[2])*k2(2)
                                        - conj(zs[3])*k2(3);
   Complex_t   ampzzx = gzzx * (k1k2*zszf - k1zf*k2zs);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma

   Complex_t   aszf = conj(as[0])*zf[0] - conj(as[1])*zf[1]
                                        - conj(as[2])*zf[2]
                                        - conj(as[3])*zf[3];
   Complex_t   k2as = conj(as[0])*k2(0) - conj(as[1])*k2(1)
                                        - conj(as[2])*k2(2)
                                        - conj(as[3])*k2(3);
   Complex_t   ampazx = gazx * (k1k2*aszf - k1zf*k2as);
#else
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel   Z: backward-going
   ANL4DVector k1   = zs.GetFourMomentum();      // s-channel   Z: backward-going
   ANL4DVector k2   = zf.GetFourMomentum();      // final state Z: forward-going 

   Double_t    k1k2 = k1*k2;
   Complex_t   k1zf = k1(0)*zf[0] - k1(1)*zf[1] - k1(2)*zf[2] - k1(3)*zf[3];
   Complex_t   zszf = zs[0]*zf[0] - zs[1]*zf[1] - zs[2]*zf[2] - zs[3]*zf[3];
   Complex_t   k2zs = zs[0]*k2(0) - zs[1]*k2(1) - zs[2]*k2(2) - zs[3]*k2(3);
   Complex_t   ampzzx = gzzx * (k1k2*zszf - k1zf*k2zs);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma

   Complex_t   aszf = as[0]*zf[0] - as[1]*zf[1] - as[2]*zf[2] - as[3]*zf[3];
   Complex_t   k2as = as[0]*k2(0) - as[1]*k2(1) - as[2]*k2(2) - as[3]*k2(3);
   Complex_t   ampazx = gazx * (k1k2*aszf - k1zf*k2as);
#endif

   //--
   // Sum of the two
   //--
   Complex_t amp = ampzzx + ampazx;
#endif /* end __HIGGS__ */
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void RSZXBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("RSZXBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSZXBases.BeamstrahlungFilename","trc500"));
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
  if (!fXBosonPtr) fXBosonPtr = new RSXBoson(fMass,
                                             fLambda,
                                             fC0,
                                             fC1,
                                             fC2,
                                             fC3);
  fXBosonPtr->DebugPrint();
  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  fWBosonPtr->DebugPrint();
  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZBosonPtr->GetMode(m); 
     if (mp && (m<fZModesLo || m>fZModesHi)) {
        mp->Lock();
     }
  }
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
  Xh_init( 7,    70.,   110.,       50, "Mz"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthF");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiF"  );
  Xh_init(10,     0.,     2.,        2, "Helin ");
  Xh_init(11,     0.,     2.,        2, "Helot ");
  Xh_init(12,     0.,     2.,        2, "X mode");
  Xh_init(13,     0.,    12.,       12, "Z mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void RSZXBases::Userout()
{
  cout << "End of RSZXBases----------------------------------- "  << endl
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
void RSZXBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][2] = {{-1, +1},
                                           {+1, -1}};
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
   fHelFinal  [0] = 0;
   fHelFinal  [1] = 0;
   fHelFinal  [2] = kFHelComb[fJCombF][0];
   fHelFinal  [3] = kFHelComb[fJCombF][1];
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
                   Double_t c3)
        : fLambda(lambda),
          fC0(c0),
          fC1(c1),
          fC2(c2),
          fC3(c3)
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
   Double_t gam = GamToVV(m1, m2, fC3, cf)/ident;
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
