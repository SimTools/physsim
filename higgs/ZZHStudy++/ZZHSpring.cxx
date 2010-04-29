//*****************************************************************************
//* =====================
//*  ZZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZH generator
//*
//* (Update Record)
//*    2010/04/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZZHSpring.h"

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

ClassImp(ZZHSpring)
ClassImp(ZZHSpringBuf)
ClassImp(ZZHBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ZZHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZHSpring::ZZHSpring(const char *name,
                     const char *title,
                     ZZHBases   *bases)
         : JSFSpring(name, title, bases)
{
  fEventBuf = new ZZHSpringBuf("ZZHSpringBuf",
                               "ZZHSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new ZZHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZZHSpring::~ZZHSpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZZHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZZHBases *bs = static_cast<ZZHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZZHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ZZHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZZHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZZHBases     *bases   = static_cast<ZZHBases *>(
		          static_cast<ZZHSpring *>(Module())->GetBases());

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
  fQ2ZZ         = bases->GetQ2ZZ();
  fCosThetaZ    = bases->GetCosThetaZ();
  fPhiZ         = bases->GetPhiZ();
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
  Int_t    idh     = 25;                          // PDG code for H
  Double_t mass    = bases->GetMass();

  Int_t    idz     = 23;                          // PDG code for Z
  // Z1
  Int_t    idf1    = bases->f1Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg1   = bases->f1Ptr->GetCharge();   // f1 charge
  Double_t m1      = bases->f1Ptr->GetMass  ();   // f1 mass
  Int_t    hel1    = bases->fHelFinal[1];         // f1 helicity
  Double_t color1  = bases->f1Ptr->GetColor();    // color factor for f1

  Int_t    idf2    = bases->f2Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg2   = bases->f2Ptr->GetCharge();   // f2 charge
  Double_t m2      = bases->f2Ptr->GetMass  ();   // f2 mass
  Int_t    hel2    = bases->fHelFinal[2];         // f2 helicity

  Int_t    islevz1 = color1 > 1. ? 201 : 0; 	  // shower level
  Int_t    icfz1   = 2;                           // color flux id
  Double_t rq2z1   = pv[1].Mag();

  // Z2
  Int_t    idf3    = bases->f1Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg3   = bases->f1Ptr->GetCharge();   // f1 charge
  Double_t m3      = bases->f1Ptr->GetMass  ();   // f1 mass
  Int_t    hel3    = bases->fHelFinal[1];         // f1 helicity
  Double_t color3  = bases->f1Ptr->GetColor();    // color factor for f1

  Int_t    idf4    = bases->f2Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg4   = bases->f2Ptr->GetCharge();   // f2 charge
  Double_t m4      = bases->f2Ptr->GetMass  ();   // f2 mass
  Int_t    hel4    = bases->fHelFinal[2];         // f2 helicity

  Int_t    islevz2 = color3 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfz2   = 3;                           // color flux id
  Double_t rq2z2   = pv[2].Mag();

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
  TVector qcm(4), qh(4), qz1(4), qz2(4);
  qh  = *qp[0];
  qz1 = *qp[1] + *qp[2];
  qz2 = *qp[3] + *qp[4];
  qcm = qh + qz1 + qz2;

  cerr << " ph=(" << qh[0] << ", "
                  << qh[1] << ", "
                  << qh[2] << ", "
                  << qh[3] << ") " << endl;
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
  new (partons[0]) JSFSpringParton( 1,  idh, mass,    0., *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2,  idz,rq2z1,    0., *qp[5], 2, 4,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3,  idz,rq2z2,    0., *qp[6], 2, 6,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton( 4, idf1,   m1, chrg1, *qp[1], 0, 0,  2, hel1,icfz1,islevz1);
  new (partons[4]) JSFSpringParton( 5, idf2,   m2, chrg2, *qp[2], 0, 0,  2, hel2,icfz1,islevz1);
  new (partons[5]) JSFSpringParton( 6, idf3,   m3, chrg3, *qp[3], 0, 0,  3, hel3,icfz2,islevz2);
  new (partons[6]) JSFSpringParton( 7, idf4,   m4, chrg4, *qp[4], 0, 0,  3, hel4,icfz2,islevz2);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZZHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZHBases::ZZHBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fMass       ( 120.),
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
           fQ2ZZ       (0.),
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
           fXQ2ZZ      (0.),
           fCosThetaZ  (0.),
           fPhiZ       (0.),
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
  stringstream ins(gJSF->Env()->GetValue("ZZHBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.CosthWRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.PhiZOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.CosthZ1FRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.PhiZ1FOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.CosthZ2FRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.PhiZ2FOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.ACC1","0.05"));
  ins >> fACC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.ACC2","0.05"));
  ins >> fACC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.ITMX1","20"));
  ins >> fITMX1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.ITMX2","40"));
  ins >> fITMX2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZHBases.NCALL","80000"));
  ins >> fNCALL;

  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fZ1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZ2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2ZZ         , 0., 1., 1, 1);
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
  fXL[7] *= TMath::Pi();
  fXU[7] *= TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaZ  , fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiZ       , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaZ1F, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiZ1F     , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaZ2F, fXL[6], fXU[6], 0, 1);
  DefineVariable(fPhiZ2F     , fXL[7], fXU[7], 0, 1);

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
ZZHBases::~ZZHBases()
{
  delete fZ1BosonPtr;
  delete fZ2BosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZZHBases::Func()
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
  Double_t qz1min = m1 + m2;
  Double_t qz1max = rs - (m3 + m4 + mh);
#ifndef __ZEROWIDTH__
  fQ2Z1 = fZ1BosonPtr->GetQ2BW(qz1min, qz1max, fXQ2Z1, weight);
#else
  fQ2Z1 = TMath::Power(fZ1BosonPtr->GetMass(),2);
  weight = kPi*fZ1BosonPtr->GetMass()*fZ1BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2z1  = TMath::Sqrt(fQ2Z1);
  Double_t qz2min = m3 + m4;
  Double_t qz2max = rs - (rq2z1 + mh);
#ifndef __ZEROWIDTH__
  fQ2Z2 = fZ2BosonPtr->GetQ2BW(qz2min, qz2max, fXQ2Z2, weight);
#else
  fQ2Z2 = TMath::Power(fZ2BosonPtr->GetMass(),2);
  weight = kPi*fZ2BosonPtr->GetMass()*fZ2BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // ZZ system
  Double_t rq2z2   = TMath::Sqrt(fQ2Z2);
  Double_t qzz2min = TMath::Power(rq2z1 + rq2z2,2);
  Double_t qzz2max = TMath::Power(rs - mh,2);
           fQ2ZZ   = qzz2min + (qzz2max - qzz2min)*fXQ2ZZ;
  bsWeight *= qzz2max - qzz2min;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch z1branch(fQ2Z1, fCosThetaZ1F, fPhiZ1F, m1*m1    , m2*m2);
  GENBranch z2branch(fQ2Z2, fCosThetaZ2F, fPhiZ2F, m3*m3    , m4*m4);
  GENBranch zzbranch(fQ2ZZ, fCosThetaZ  , fPhiZ  , &z1branch, &z2branch);
  GENBranch cmbranch(s    , fCosTheta   , fPhi   , mh*mh    , &zzbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2ZZ), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaZ        , (bsWeight*sigma));
  Xh_fill( 6, fPhiZ             , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Z1), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaZ1F      , (bsWeight*sigma));
  Xh_fill( 9, fPhiZ1F           , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2Z2), (bsWeight*sigma));
  Xh_fill(11, fCosThetaZ2F      , (bsWeight*sigma));
  Xh_fill(12, fPhiZ2F           , (bsWeight*sigma));
  Xh_fill(13, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(14, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(15, (Double_t)fZ1Mode , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fZ2Mode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZZHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2h  = cmbranch.GetM12();
  Double_t q2zz = cmbranch.GetM22();
  Double_t cosh = cmbranch.GetCosTheta();
  Double_t phih = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2h, q2zz, cmframe, cosh, phih, 0);
  ANL4DVector ph  = phaseCM.GetFourMomentum(0);
  ANL4DVector pzz = phaseCM.GetFourMomentum(1);
  Double_t betah  = phaseCM.GetBetaBar();
  if (betah <= 0.) return 0.;
  fP[0] = ph;
  fM[0] = fMass;

  GENBranch &zzbranch = * cmbranch.GetBranchPtr(1);
  Double_t cosz = zzbranch.GetCosTheta();
  Double_t phiz = zzbranch.GetPhi     ();
  Double_t qz12 = zzbranch.GetM12();
  Double_t qz22 = zzbranch.GetM22();
  GENPhase2 phaseZZ(pzz, qz12, qz22, cmframe, cosz, phiz, 1);
  ANL4DVector pz1 = phaseZZ.GetFourMomentum(0);
  ANL4DVector pz2 = phaseZZ.GetFourMomentum(1);
  Double_t betaz  = phaseZZ.GetBetaBar();
  if (betaz <= 0.) return 0.;

  GENBranch &z1branch = * zzbranch.GetBranchPtr(0);
  Double_t cosz1f = z1branch.GetCosTheta();
  Double_t phiz1f = z1branch.GetPhi     ();
  Double_t m12    = z1branch.GetM12();
  Double_t m22    = z1branch.GetM22();
  GENPhase2 phaseZ1(pz1, m12, m22, cmframe, cosz1f, phiz1f, 1);
  fP[1] = phaseZ1.GetFourMomentum(0);
  fP[2] = phaseZ1.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m12);
  fM[2] = TMath::Sqrt(m22);
  Double_t betaz1f = phaseZ1.GetBetaBar();
  if (betaz1f <= 0.) return 0.;

  GENBranch &z2branch = * zzbranch.GetBranchPtr(1);
  Double_t cosz2f = z2branch.GetCosTheta();
  Double_t phiz2f = z2branch.GetPhi     ();
  Double_t m32    = z2branch.GetM12();
  Double_t m42    = z2branch.GetM22();
  GENPhase2 phaseZ2(pz2, m32, m42, cmframe, cosz2f, phiz2f, 1);
  fP[3] = phaseZ2.GetFourMomentum(0);
  fP[4] = phaseZ2.GetFourMomentum(1);
  fM[3] = TMath::Sqrt(m32);
  fM[4] = TMath::Sqrt(m42);
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
  for (Int_t i=0; i<5; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qz1 = fP[1] + fP[2];
  ANL4DVector qz2 = fP[3] + fP[4];
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
  ANL4DVector pcm = fP[0] + qz1 + qz2;
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

  Double_t identp = 1./2.;                                     // identical particle factor
  Double_t dPhase = kFact * betah * betaz * betaz1f * betaz2f; // phase space factor
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
Double_t ZZHBases::AmpSquared(GENBranch &cmbranch)
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
Complex_t ZZHBases::FullAmplitude()
{
   Double_t gamz   = fZ1BosonPtr->GetWidth();

   Double_t qf1    = f1Ptr->GetCharge();
   Double_t t3f1   = f1Ptr->GetISpin();
   Double_t glfz1  = -kGz*(t3f1 - qf1*kSin2W);
   Double_t grfz1  = -kGz*(     - qf1*kSin2W);

   Double_t qf3    = f3Ptr->GetCharge();
   Double_t t3f3   = f3Ptr->GetISpin();
   Double_t glfz3  = -kGz*(t3f3 - qf3*kSin2W);
   Double_t grfz3  = -kGz*(     - qf3*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

   HELScalar  hs(fP[0]); // higgs

   HELFermion f1 (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing); // f1
   HELFermion f2b(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming); // f2b
   HELVector  z1 (f2b, f1, glfz1, grfz1, kM_z, gamz);            // Z1

   HELFermion f3 (fP[3], fM[3], fHelFinal [3], +1, kIsOutgoing); // f3
   HELFermion f4b(fP[4], fM[4], fHelFinal [4], -1, kIsIncoming); // f4d
   HELVector  z2 (f4b, f3, glfz3, grfz3, kM_z, gamz);            // Z2

   Complex_t amp = AmpEEtoZZH(em, ep, z1, z2, hs);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZZH()
// --------------------------
Complex_t ZZHBases::AmpEEtoZZH(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &z1,
                               const HELVector  &z2,
                               const HELScalar  &hs)
{
   Double_t  qe    = -1.;
   Double_t  game  = 0.;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamz  = fZ1BosonPtr->GetWidth();

   Double_t  gzzh  = kGz*kM_z;

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
   // ZZH Production Amplitude
   //---------------------------
   HELFermion emz1  (em, z1, glze, grze, kM_e, game);
   HELVector  z2star(z2, hs, gzzh, kM_z, gamz);
   HELVertex  ampz1z2s(emz1, ep, z2star, glze, grze);

   HELFermion emz2(em, z2, glze, grze, kM_e, game);
   HELVector  z1star(z1, hs, gzzh, kM_z, gamz);
   HELVertex  ampz2z1s(emz2, ep, z1star, glze, grze);

   HELFermion emz1s (em, z1star, glze, grze, kM_e, game);
   HELVertex  ampz1sz2(emz1s, ep, z2, glze, grze);

   HELFermion emz2s (em, z2star, glze, grze, kM_e, game);
   HELVertex  ampz2sz1(emz2s, ep, z1, glze, grze);

   Complex_t amp = ampz1z2s + ampz2z1s + ampz1sz2 + ampz2sz1;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZZHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZZHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZZHBases.BeamstrahlungFilename","trc500"));
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
  Xh_init( 4,     0., fEcmInit    , 50, "Qzz"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosW"   );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiZ"   );
  Xh_init( 7,    70.,   110.,       50, "Mz1"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosZ1F" );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiZ1F" );
  Xh_init(10,    70.,   110.,       50, "Mz2"    );
  Xh_init(11, fXL[6], fXU[6],       50, "CosZ2F" );
  Xh_init(12, fXL[7], fXU[7],       50, "PhiZ2F" );
  Xh_init(13,     0.,     2.,        2, "Helin " );
  Xh_init(14,     0.,     1.,        2, "Helot " );
  Xh_init(15,     0.,    12.,       12, "Z1 mode");
  Xh_init(16,     0.,    12.,       12, "Z2 mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZZHBases::Userout()
{
  cout << "End of ZZHBases----------------------------------- "  << endl
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
void ZZHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][5] = {{0, -1, +1, -1, +1},
                                           {0, -1, +1, +1, -1},
                                           {0, +1, -1, -1, +1},
                                           {0, +1, -1, +1, -1}};
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
