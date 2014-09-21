//*****************************************************************************
//* =====================
//*  ZZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZ generator
//*
//* (Update Record)
//*    2014/09/21  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZZSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZZSpring)
ClassImp(ZZSpringBuf)
ClassImp(ZZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ZZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZSpring::ZZSpring(const char *name,
                   const char *title,
                   ZZBases    *bases)
        : JSFSpring(name, title, bases)
{
  fEventBuf = new ZZSpringBuf("ZZSpringBuf",
                              "ZZSpring event buffer",
                              this);
  if (!bases) { 
    SetBases(new ZZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZZSpring::~ZZSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZZBases *bs = static_cast<ZZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZZBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ZZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZZBases      *bases   = static_cast<ZZBases *>(
                          static_cast<ZZSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 6;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[1] = bases->fP[1];
  pv[2] = bases->fP[2];
  pv[3] = bases->fP[3];
  pv[4] = pv[0] + pv[1];
  pv[5] = pv[2] + pv[3];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2Z1         = bases->GetQ2Z1();
  fCosThetaZ1F  = bases->GetCosThetaZ1F();
  fPhiZ1F       = bases->GetPhiZ1F();
  fQ2Z2         = bases->GetQ2Z2();
  fCosThetaZ2F  = bases->GetCosThetaZ2F();
  fPhiZ2F       = bases->GetPhiZ2F();
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
  // Z1
  Int_t    idf1    = bases->f1Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg1   = bases->f1Ptr->GetCharge();   // f1 charge
  Double_t m1      = bases->f1Ptr->GetMass  ();   // f1 mass
  Int_t    hel1    = bases->fHelFinal[0];         // f1 helicity
  Double_t color1  = bases->f1Ptr->GetColor();    // color factor for f1

  Int_t    idf2    = bases->f2Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg2   = bases->f2Ptr->GetCharge();   // f2 charge
  Double_t m2      = bases->f2Ptr->GetMass  ();   // f2 mass
  Int_t    hel2    = bases->fHelFinal[1];         // f2 helicity

  Int_t    islevz1 = color1 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfz1   = 2;                           // color flux id
  Double_t rq2z1   = pv[4].Mag();

  // Z2
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f1 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f1 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f1 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f1

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f2 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f2 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f2 helicity

  Int_t    islevz2 = color3 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfz2   = 3;                           // color flux id
  Double_t rq2z2   = pv[5].Mag();

#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 2 pid=" << idf1 << " m=" << m1 << " Q=" << chrg1 
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 3 pid=" << idf2 << " m=" << m2 << " Q=" << chrg2 
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 6 pid=" << idf3 << " m=" << m3 << " Q=" << chrg3 
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 7 pid=" << idf4 << " m=" << m4 << " Q=" << chrg4
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  TVector qcm(4), qz1(4), qz2(4);
  qz1 = *qp[0] + *qp[1];
  qz2 = *qp[2] + *qp[3];
  qcm = qz1 + qz2;

  cerr << " pz1=(" << qz1[0] << ", "
                   << qz1[1] << ", "
                   << qz1[2] << ", "
                   << qz1[3] << ") " << endl;
  cerr << " pz2=(" << qz2[0] << ", "
                   << qz2[1] << ", "
                   << qz2[2] << ", "
                   << qz2[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                                No.  PID  Mass  Charge   pv   Nd 1st Mom  hel  col  shower
  new (partons[0]) JSFSpringParton( 1,  idz,rq2z1,   -1., *qp[4], 2, 3,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2,  idz,rq2z2,   +1., *qp[5], 2, 5,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3, idf1,   m1, chrg1, *qp[0], 0, 0,  1, hel1,icfz1,islevz1);
  new (partons[3]) JSFSpringParton( 4, idf2,   m2, chrg2, *qp[1], 0, 0,  1, hel2,icfz1,islevz1);
  new (partons[4]) JSFSpringParton( 5, idf3,   m3, chrg3, *qp[2], 0, 0,  2, hel3,icfz2,islevz2);
  new (partons[5]) JSFSpringParton( 6, idf4,   m4, chrg4, *qp[3], 0, 0,  2, hel4,icfz2,islevz2);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZBases::ZZBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fEcmInit    (1000.),
           fISR        ( 1),
           fBeamStr    ( 1),
           fBeamWidth  (0.002),
           fPole       (0.),
           fZ1ModesLo  ( 1),
           fZ1ModesHi  (12),
           fZ2ModesLo  ( 1),
           fZ2ModesHi  (12),
	   fNCALL      (80000),
	   fACC1       (0.05),
	   fACC2       (0.05),
	   fITMX1      (20),
	   fITMX2      (40),
           fZ1BosonPtr ( 0),
           fZ2BosonPtr ( 0),
           fZBoost     (0.),
           fEcmIP      (fEcmInit),
           fQ2Z1       (0.),
           fQ2Z2       (0.),
           fZ1ModePtr  (0),
           f1Ptr       (0),
           f2Ptr       (0),
           fZ2ModePtr  (0),
           f3Ptr       (0),
           f4Ptr       (0),
           fCosTheta   (0.),
           fPhi        (0.),
           fXQ2Z1      (0.),
           fCosThetaZ1F(0.),
           fPhiZ1F     (0.),
           fXQ2Z2      (0.),
           fCosThetaZ2F(0.),
           fPhiZ2F     (0.),
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
  stringstream ins(gJSF->Env()->GetValue("ZZBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.CosthZ1FRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.PhiZ1FOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.CosthZ2FRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.PhiZ2FOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.ACC1","0.05"));
  ins >> fACC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.ACC2","0.05"));
  ins >> fACC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.ITMX1","20"));
  ins >> fITMX1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.ITMX2","40"));
  ins >> fITMX2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZBases.NCALL","80000"));
  ins >> fNCALL;

  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fZ1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZ2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2Z1         , 0., 1., 0, 1);
  DefineVariable(fXQ2Z2         , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] *= TMath::Pi();
  fXU[1] *= TMath::Pi();
  fXL[3] *= TMath::Pi();
  fXU[3] *= TMath::Pi();
  fXL[5] *= TMath::Pi();
  fXU[5] *= TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaZ1F, fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiZ1F     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaZ2F, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiZ2F     , fXL[5], fXU[5], 0, 1);

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
ZZBases::~ZZBases()
{
  delete fZ1BosonPtr;
  delete fZ2BosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZZBases::Func()
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

  GENDecayMode *fZ1ModePtr = fZ1BosonPtr->PickMode(fZ1DecayMode, weight, fZ1Mode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fZ2ModePtr = fZ2BosonPtr->PickMode(fZ2DecayMode, weight, fZ2Mode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < m1 + m2 + m3 + m4) {
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
  Double_t qz1min = m1 + m2;
  Double_t qz1max = rs - (m3 + m4);
#ifndef __ZEROWIDTH__
  fQ2Z1 = fZ1BosonPtr->GetQ2BW(qz1min, qz1max, fXQ2Z1, weight);
#else
  fQ2Z1 = TMath::Power(fZ1BosonPtr->GetMass(),2);
  weight = kPi*fZ1BosonPtr->GetMass()*fZ1BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2z1  = TMath::Sqrt(fQ2Z1);
  Double_t qz2min = m3 + m4;
  Double_t qz2max = rs - rq2z1;
#ifndef __ZEROWIDTH__
  fQ2Z2 = fZ2BosonPtr->GetQ2BW(qz2min, qz2max, fXQ2Z2, weight);
#else
  fQ2Z2 = TMath::Power(fZ2BosonPtr->GetMass(),2);
  weight = kPi*fZ2BosonPtr->GetMass()*fZ2BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch z1branch(fQ2Z1, fCosThetaZ1F, fPhiZ1F, m1*m1    , m2*m2);
  GENBranch z2branch(fQ2Z2, fCosThetaZ2F, fPhiZ2F, m3*m3    , m4*m4);
  GENBranch cmbranch(s    , fCosTheta   , fPhi   , &z1branch, &z2branch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2Z1), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaZ1F      , (bsWeight*sigma));
  Xh_fill( 6, fPhiZ1F           , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Z2), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaZ2F      , (bsWeight*sigma));
  Xh_fill( 9, fPhiZ2F           , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(12, (Double_t)fZ1Mode , (bsWeight*sigma));
  Xh_fill(13, (Double_t)fZ2Mode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZZBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2z1 = cmbranch.GetM12();
  Double_t q2z2 = cmbranch.GetM22();
  Double_t cosz = cmbranch.GetCosTheta();
  Double_t phiz = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2z1, q2z2, cmframe, cosz, phiz, 0);
  ANL4DVector pz1 = phaseCM.GetFourMomentum(0);
  ANL4DVector pz2 = phaseCM.GetFourMomentum(1);
  Double_t betaz  = phaseCM.GetBetaBar();
  if (betaz <= 0.) return 0.;

  GENBranch &z1branch = *cmbranch.GetBranchPtr(0);
  Double_t cosz1f = z1branch.GetCosTheta();
  Double_t phiz1f = z1branch.GetPhi     ();
  Double_t m12    = z1branch.GetM12();
  Double_t m22    = z1branch.GetM22();
  GENPhase2 phaseZ1(pz1, m12, m22, cmframe, cosz1f, phiz1f, 1);
  fP[0] = phaseZ1.GetFourMomentum(0);
  fP[1] = phaseZ1.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betaz1f = phaseZ1.GetBetaBar();
  if (betaz1f <= 0.) return 0.;

  GENBranch &z2branch = *cmbranch.GetBranchPtr(1);
  Double_t cosz2f = z2branch.GetCosTheta();
  Double_t phiz2f = z2branch.GetPhi     ();
  Double_t m32    = z2branch.GetM12();
  Double_t m42    = z2branch.GetM22();
  GENPhase2 phaseZ2(pz2, m32, m42, cmframe, cosz2f, phiz2f, 1);
  fP[2] = phaseZ2.GetFourMomentum(0);
  fP[3] = phaseZ2.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betaz2f = phaseZ2.GetBetaBar();
  if (betaz2f <= 0.) return 0.;

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
  for (Int_t i=0; i<4; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qz1 = fP[0] + fP[1];
  ANL4DVector qz2 = fP[2] + fP[3];
  cerr << " qz1 = (" 
       << qz1.E () << ","
       << qz1.Px() << ","
       << qz1.Py() << ","
       << qz1.Pz() << ")" << endl;
  cerr << " qz2 = (" 
       << qz2.E () << ","
       << qz2.Px() << ","
       << qz2.Py() << ","
       << qz2.Pz() << ")" << endl;
  ANL4DVector pcm = qz1 + qz2;
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
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1./2.;                             // identical particle factor
  Double_t dPhase = kFact * betaz * betaz1f * betaz2f; // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);               // beam flux factor
  Double_t spin   = 1./2.;                             // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase;      // in [1/GeV^2]
           sigma *= kGeV2fb;                                   // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ZZBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor();

  Complex_t amp    = FullAmplitude();
  Double_t  amp2   = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ZZBases::FullAmplitude()
{
   Double_t mz     = fZ1BosonPtr->GetMass();
   Double_t gamz   = fZ1BosonPtr->GetWidth();

   Double_t qf1    = f1Ptr->GetCharge();
   Double_t t3f1   = f1Ptr->GetISpin();
   Double_t glzf1  = -kGz*(t3f1 - qf1*kSin2W);
   Double_t grzf1  = -kGz*(     - qf1*kSin2W);

   Double_t qf3    = f3Ptr->GetCharge();
   Double_t t3f3   = f3Ptr->GetISpin();
   Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
   Double_t grzf3  = -kGz*(     - qf3*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

   HELFermion f1 (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // f1
   HELFermion f2b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // f2b
   HELVector  z1(f2b, f1, glzf1, grzf1, mz, gamz);               // Z1

   HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fu
   HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // fdbar
   HELVector  z2(f4b, f3, glzf3, grzf3, mz, gamz);               // Z2

   Complex_t amp = AmpEEtoZZ(em, ep, z1, z2);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZZ()
// --------------------------
Complex_t ZZBases::AmpEEtoZZ(const HELFermion &em,
                             const HELFermion &ep,
                             const HELVector  &z1,
                             const HELVector  &z2)
{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamz  = fZ1BosonPtr->GetWidth();

   Double_t  me    = kMass[0][1][0];
   Double_t  game  = 0.;

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
   // ZZ Production Amplitude
   //---------------------------
   //--
   // T-channel
   //--
   HELFermion emz1 (em, z1, glze, grze, me, game);
   HELVertex  amp01(emz1, ep, z2, glze, grze);

   //--
   // U-channel ne
   //--
   HELFermion emz2 (em, z2, glze, grze, me, game);
   HELVertex  amp02(emz2, ep, z1, glze, grze);

   //--
   // Sum
   //--
   Complex_t amp = amp01 + amp02;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZZBases.BeamstrahlungFilename","trc500"));
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
  if (!fZ1BosonPtr) fZ1BosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZ1BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZ1BosonPtr->GetMode(m); 
     if (mp && (m<fZ1ModesLo || m>fZ1ModesHi)) {
        mp->Lock();
     }
  }
  fZ1BosonPtr->DebugPrint();

  if (!fZ2BosonPtr) fZ2BosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZ2BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZ2BosonPtr->GetMode(m); 
     if (mp && (m<fZ2ModesLo || m>fZ2ModesHi)) {
        mp->Lock();
     }
  }
  fZ2BosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4,    70.,   110.,       50, "Mz1"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosF1"  );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF1"  );
  Xh_init( 7,    70.,   110.,       50, "Mz2"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosF3"  );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiF3"  );
  Xh_init(10,     0.,     2.,        2, "Helin " );
  Xh_init(11,     0.,     4.,        4, "Helot " );
  Xh_init(12,     0.,    12.,       12, "Z1 mode");
  Xh_init(13,     0.,    12.,       12, "Z2 mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZZBases::Userout()
{
  cout << "End of ZZBases----------------------------------- "  << endl
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
void ZZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1}};
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
