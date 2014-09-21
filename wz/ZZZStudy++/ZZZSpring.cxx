//*****************************************************************************
//* =====================
//*  ZZZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZZ generator
//*
//* (Update Record)
//*    2014/09/21  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZZZSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZZZSpring)
ClassImp(ZZZSpringBuf)
ClassImp(ZZZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ZZZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZZSpring::ZZZSpring(const char *name,
                     const char *title,
                     ZZZBases   *bases)
        : JSFSpring(name, title, bases)
{
  fEventBuf = new ZZZSpringBuf("ZZZSpringBuf",
                               "ZZZSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new ZZZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZZZSpring::~ZZZSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZZZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZZZBases *bs = static_cast<ZZZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZZZBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ZZZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZZZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZZZBases     *bases   = static_cast<ZZZBases *>(
                          static_cast<ZZZSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 9;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[1] = bases->fP[1];
  pv[2] = bases->fP[2];
  pv[3] = bases->fP[3];
  pv[4] = bases->fP[4];
  pv[5] = bases->fP[5];
  pv[6] = pv[0] + pv[1];
  pv[7] = pv[2] + pv[3];
  pv[8] = pv[4] + pv[5];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2ZZ         = bases->GetQ2ZZ();
  fCosThetaZ1   = bases->GetCosThetaZ1();
  fPhiZ1        = bases->GetPhiZ1();
  fQ2Z1         = bases->GetQ2Z1();
  fCosThetaZ1F  = bases->GetCosThetaZ1F();
  fPhiZ1F       = bases->GetPhiZ1F();
  fQ2Z2         = bases->GetQ2Z2();
  fCosThetaZ2F  = bases->GetCosThetaZ2F();
  fPhiZ2F       = bases->GetPhiZ2F();
  fQ2Z          = bases->GetQ2Z();
  fCosThetaZF   = bases->GetCosThetaZF();
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
  Double_t rq2z1   = pv[6].Mag();

  // Z2
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f3 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f4 helicity

  Int_t    islevz2 = color3 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfz2   = 3;                           // color flux id
  Double_t rq2z2   = pv[7].Mag();

  // Z
  Int_t    idf5    = bases->f5Ptr->GetPID   ();   // PDG code for f5
  Double_t chrg5   = bases->f5Ptr->GetCharge();   // f5 charge
  Double_t m5      = bases->f5Ptr->GetMass  ();   // f5 mass
  Int_t    hel5    = bases->fHelFinal[4];         // f5 helicity
  Double_t color5  = bases->f5Ptr->GetColor();    // color factor for f5

  Int_t    idf6    = bases->f6Ptr->GetPID   ();   // PDG code for f6
  Double_t chrg6   = bases->f6Ptr->GetCharge();   // f6 charge
  Double_t m6      = bases->f6Ptr->GetMass  ();   // f6 mass
  Int_t    hel6    = bases->fHelFinal[5];         // f6 helicity

  Int_t    islevz  = color5 > 1. ? 301 : 0;  	  // shower level
  Int_t    icfz    = 4;                           // color flux id
  Double_t rq2z    = pv[8].Mag();

#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idf1 << " m=" << m1 << " Q=" << chrg1 
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 2 pid=" << idf2 << " m=" << m2 << " Q=" << chrg2 
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 3 pid=" << idf3 << " m=" << m3 << " Q=" << chrg3 
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 4 pid=" << idf4 << " m=" << m4 << " Q=" << chrg4
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  cerr << " 5 pid=" << idf5 << " m=" << m5 << " Q=" << chrg5 
       << " pv=(" << (*qp[4])(0) << ", "
                  << (*qp[4])(1) << ", "
                  << (*qp[4])(2) << ", "
                  << (*qp[4])(3) << ") " << endl;
  cerr << " 6 pid=" << idf6 << " m=" << m6 << " Q=" << chrg6
       << " pv=(" << (*qp[5])(0) << ", "
                  << (*qp[5])(1) << ", "
                  << (*qp[5])(2) << ", "
                  << (*qp[5])(3) << ") " << endl;
  TVector qcm(4), qz1(4), qz2(4), qz(4);
  qz1 = *qp[0] + *qp[1];
  qz2 = *qp[2] + *qp[3];
  qz  = *qp[4] + *qp[5];
  qcm = qz1 + qz2 + qz;

  cerr << " pz1=(" << qz1[0] << ", "
                   << qz1[1] << ", "
                   << qz1[2] << ", "
                   << qz1[3] << ") " << endl;
  cerr << " pz2=(" << qz2[0] << ", "
                   << qz2[1] << ", "
                   << qz2[2] << ", "
                   << qz2[3] << ") " << endl;
  cerr << " pz =(" << qz [0] << ", "
                   << qz [1] << ", "
                   << qz [2] << ", "
                   << qz [3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                                No.  PID  Mass  Charge   pv   Nd 1st Mom  hel  col  shower
  new (partons[0]) JSFSpringParton( 1,  idz,rq2z1,    0., *qp[6], 2, 4,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2,  idz,rq2z2,    0., *qp[7], 2, 6,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3,  idz, rq2z,    0., *qp[8], 2, 8,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton( 4, idf1,   m1, chrg1, *qp[0], 0, 0,  1, hel1,icfz1,islevz1);
  new (partons[4]) JSFSpringParton( 5, idf2,   m2, chrg2, *qp[1], 0, 0,  1, hel2,icfz1,islevz1);
  new (partons[5]) JSFSpringParton( 6, idf3,   m3, chrg3, *qp[2], 0, 0,  2, hel3,icfz2,islevz2);
  new (partons[6]) JSFSpringParton( 7, idf4,   m4, chrg4, *qp[3], 0, 0,  2, hel4,icfz2,islevz2);
  new (partons[7]) JSFSpringParton( 8, idf5,   m5, chrg5, *qp[4], 0, 0,  3, hel5, icfz, islevz);
  new (partons[8]) JSFSpringParton( 9, idf6,   m6, chrg6, *qp[5], 0, 0,  3, hel6, icfz, islevz);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZZZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZZZBases::ZZZBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fEcmInit    (500.),
           fISR        ( 1),
           fBeamStr    ( 1),
           fBeamWidth  (0.002),
           fPole       (0.),
           fZ1ModesLo  ( 1),
           fZ1ModesHi  (12),
           fZ2ModesLo  ( 1),
           fZ2ModesHi  (12),
           fZModesLo   ( 1),
           fZModesHi   (12),
	   fNCALL      (80000),
	   fACC1       (0.05),
	   fACC2       (0.05),
	   fITMX1      (20),
	   fITMX2      (40),
           fZ1BosonPtr ( 0),
           fZ2BosonPtr ( 0),
           fZBosonPtr  ( 0),
           fZBoost     (0.),
           fEcmIP      (fEcmInit),
           fQ2ZZ       (0.),
           fQ2Z1       (0.),
           fQ2Z2       (0.),
           fQ2Z        (0.),
           fZ1ModePtr  (0),
           f1Ptr       (0),
           f2Ptr       (0),
           fZ2ModePtr  (0),
           f3Ptr       (0),
           f4Ptr       (0),
           fZModePtr   (0),
           f5Ptr       (0),
           f6Ptr       (0),
           fCosTheta   (0.),
           fPhi        (0.),
           fXQ2ZZ      (0.),
           fCosThetaZ1 (0.),
           fPhiZ1      (0.),
           fXQ2Z1      (0.),
           fCosThetaZ1F(0.),
           fPhiZ1F     (0.),
           fXQ2Z2      (0.),
           fCosThetaZ2F(0.),
           fPhiZ2F     (0.),
           fXQ2Z       (0.),
           fCosThetaZF (0.),
           fPhiZF      (0.),
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
  stringstream ins(gJSF->Env()->GetValue("ZZZBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.CosthZ1Range","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.PhiZ1OverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.CosthZ1FRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.PhiZ1FOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.CosthZ2FRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.PhiZ2FOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.CosthZFRange","-1.0 1.0"));
  ins >> fXL[8] >> fXU[8];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.PhiZFOverPiRange","0.0 2.0"));
  ins >> fXL[9] >> fXU[9];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Z1ModesLo","1"));      // Z1 decay mode lo
  ins >> fZ1ModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Z1ModesHi","12"));     // Z1 decay mode hi
  ins >> fZ1ModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Z2ModesLo","1"));      // Z2 decay mode lo
  ins >> fZ2ModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.Z2ModesHi","12"));     // Z2 decay mode hi
  ins >> fZ2ModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZZZBases.NCALL","80000"));
  ins >> fNCALL;

  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fZ1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZ2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2ZZ         , 0., 1., 1, 1);
  DefineVariable(fXQ2Z1         , 0., 1., 0, 1);
  DefineVariable(fXQ2Z2         , 0., 1., 0, 1);
  DefineVariable(fXQ2Z          , 0., 1., 0, 1);
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
  fXL[9] *= TMath::Pi();
  fXU[9] *= TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaZ1 , fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiZ1      , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaZ1F, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiZ1F     , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaZ2F, fXL[6], fXU[6], 1, 1);
  DefineVariable(fPhiZ2F     , fXL[7], fXU[7], 0, 1);
  DefineVariable(fCosThetaZF , fXL[8], fXU[8], 1, 1);
  DefineVariable(fPhiZF      , fXL[9], fXU[9], 0, 1);

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
ZZZBases::~ZZZBases()
{
  delete fZ1BosonPtr;
  delete fZ2BosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZZZBases::Func()
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

  GENDecayMode *fZModePtr  = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f5Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f6Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m5   = f5Ptr->GetMass();
  Double_t m6   = f6Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < m1 + m2 + m3 + m4 + m5 + m6) {
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
  Double_t qz1max = rs - (m3 + m4 + m5 + m6);
#ifndef __ZEROWIDTH__
  fQ2Z1 = fZ1BosonPtr->GetQ2BW(qz1min, qz1max, fXQ2Z1, weight);
#else
  fQ2Z1 = TMath::Power(fZ1BosonPtr->GetMass(),2);
  weight = kPi*fZ1BosonPtr->GetMass()*fZ1BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2z1  = TMath::Sqrt(fQ2Z1);
  Double_t qz2min = m3 + m4;
  Double_t qz2max = rs - (rq2z1 + m5 + m6);
#ifndef __ZEROWIDTH__
  fQ2Z2 = fZ2BosonPtr->GetQ2BW(qz2min, qz2max, fXQ2Z2, weight);
#else
  fQ2Z2 = TMath::Power(fZ2BosonPtr->GetMass(),2);
  weight = kPi*fZ2BosonPtr->GetMass()*fZ2BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2z2  = TMath::Sqrt(fQ2Z2);
  Double_t qzmin  = m5 + m6;
  Double_t qzmax  = rs - (rq2z1 + rq2z2);
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qzmin, qzmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2z    = TMath::Sqrt(fQ2Z);
  Double_t qzzmin  = rq2z1 + rq2z2;
  Double_t qzzmax  = rs - rq2z;
  Double_t q2zzmin = qzzmin*qzzmin;
  Double_t q2zzmax = qzzmax*qzzmax;
  fQ2ZZ  = q2zzmin + (q2zzmax - q2zzmin)*fXQ2ZZ;
  weight = q2zzmax - q2zzmin;
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch z1branch(fQ2Z1, fCosThetaZ1F, fPhiZ1F, m1*m1    , m2*m2);
  GENBranch z2branch(fQ2Z2, fCosThetaZ2F, fPhiZ2F, m3*m3    , m4*m4);
  GENBranch zbranch (fQ2Z , fCosThetaZF , fPhiZF , m5*m5    , m6*m6);
  GENBranch zzbranch(fQ2ZZ, fCosThetaZ1 , fPhiZ1 , &z1branch, &z2branch);
  GENBranch cmbranch(s    , fCosTheta   , fPhi   , &zzbranch, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2ZZ), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaZ1       , (bsWeight*sigma));
  Xh_fill( 6, fPhiZ1            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Z1), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaZ1F      , (bsWeight*sigma));
  Xh_fill( 9, fPhiZ1F           , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2Z2), (bsWeight*sigma));
  Xh_fill(11, fCosThetaZ2F      , (bsWeight*sigma));
  Xh_fill(12, fPhiZ2F           , (bsWeight*sigma));
  Xh_fill(13, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill(14, fCosThetaZF       , (bsWeight*sigma));
  Xh_fill(15, fPhiZF            , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fZ1Mode , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fZ2Mode , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fZMode  , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZZZBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2zz  = cmbranch.GetM12();
  Double_t q2z   = cmbranch.GetM22();
  Double_t coszz = cmbranch.GetCosTheta();
  Double_t phizz = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2zz, q2z, cmframe, coszz, phizz, 0);
  ANL4DVector pzz = phaseCM.GetFourMomentum(0);
  ANL4DVector pz  = phaseCM.GetFourMomentum(1);
  Double_t betazz = phaseCM.GetBetaBar();
  if (betazz <= 0.) return 0.;

  GENBranch &zzbranch = *cmbranch.GetBranchPtr(0);
  Double_t q2z1  = zzbranch.GetM12();
  Double_t q2z2  = zzbranch.GetM22();
  Double_t cosz1 = zzbranch.GetCosTheta();
  Double_t phiz1 = zzbranch.GetPhi     ();
  GENPhase2 phaseZZ(pzz, q2z1, q2z2, cmframe, cosz1, phiz1, 1);
  ANL4DVector pz1 = phaseZZ.GetFourMomentum(0);
  ANL4DVector pz2 = phaseZZ.GetFourMomentum(1);
  Double_t betaz1  = phaseZZ.GetBetaBar();
  if (betaz1 <= 0.) return 0.;

  GENBranch &z1branch = *zzbranch.GetBranchPtr(0);
  Double_t cosz1f  = z1branch.GetCosTheta();
  Double_t phiz1f  = z1branch.GetPhi     ();
  Double_t m12     = z1branch.GetM12();
  Double_t m22     = z1branch.GetM22();
  GENFrame zzframe = phaseZZ.GetFrame();
  GENPhase2 phaseZ1(pz1, m12, m22, zzframe, cosz1f, phiz1f, 1);
  fP[0] = phaseZ1.GetFourMomentum(0);
  fP[1] = phaseZ1.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betaz1f = phaseZ1.GetBetaBar();
  if (betaz1f <= 0.) return 0.;

  GENBranch &z2branch = *zzbranch.GetBranchPtr(1);
  Double_t cosz2f = z2branch.GetCosTheta();
  Double_t phiz2f = z2branch.GetPhi     ();
  Double_t m32    = z2branch.GetM12();
  Double_t m42    = z2branch.GetM22();
  GENPhase2 phaseZ2(pz2, m32, m42, zzframe, cosz2f, phiz2f, 1);
  fP[2] = phaseZ2.GetFourMomentum(0);
  fP[3] = phaseZ2.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betaz2f = phaseZ2.GetBetaBar();
  if (betaz2f <= 0.) return 0.;

  GENBranch &zbranch = *cmbranch.GetBranchPtr(1);
  Double_t coszf  = zbranch.GetCosTheta();
  Double_t phizf  = zbranch.GetPhi     ();
  Double_t m52    = zbranch.GetM12();
  Double_t m62    = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m52, m62, cmframe, coszf, phizf, 1);
  fP[4] = phaseZ.GetFourMomentum(0);
  fP[5] = phaseZ.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(m52);
  fM[5] = TMath::Sqrt(m62);
  Double_t betazf = phaseZ.GetBetaBar();
  if (betazf <= 0.) return 0.;

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
  for (Int_t i=0; i<6; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qz1 = fP[0] + fP[1];
  ANL4DVector qz2 = fP[2] + fP[3];
  ANL4DVector qz  = fP[4] + fP[5];
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
  cerr << " qz  = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector pcm = qz1 + qz2 + qz;
  cerr << " mz1 = " << qz1.GetMass() << endl;
  cerr << " mz2 = " << qz2.GetMass() << endl;
  cerr << " mz  = " << qz .GetMass() << endl;
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
  static const Int_t    kNbr  = 5;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1./(1.*2.*3);                      // identical particle factor
  Double_t dPhase = kFact * betazz * betaz1 
	                  * betaz1f * betaz2f * betazf;// phase space factor
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
Double_t ZZZBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor() * f5Ptr->GetColor();

  Complex_t amp    = FullAmplitude();
  Double_t  amp2   = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ZZZBases::FullAmplitude()
{
   Double_t mz     = fZBosonPtr->GetMass();
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf1    = f1Ptr->GetCharge();
   Double_t t3f1   = f1Ptr->GetISpin();
   Double_t glzf1  = -kGz*(t3f1 - qf1*kSin2W);
   Double_t grzf1  = -kGz*(     - qf1*kSin2W);

   Double_t qf3    = f3Ptr->GetCharge();
   Double_t t3f3   = f3Ptr->GetISpin();
   Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
   Double_t grzf3  = -kGz*(     - qf3*kSin2W);

   Double_t qf5    = f5Ptr->GetCharge();
   Double_t t3f5   = f5Ptr->GetISpin();
   Double_t glzf5  = -kGz*(t3f5 - qf5*kSin2W);
   Double_t grzf5  = -kGz*(     - qf5*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

   HELFermion f1 (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // f1
   HELFermion f2b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // f2b
   HELVector  z1(f2b, f1, glzf1, grzf1, mz, gamz);               // Z1

   HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // f3
   HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // f4b
   HELVector  z2(f4b, f3, glzf3, grzf3, mz, gamz);               // Z2

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // f5
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // f6b
   HELVector  zf(f6b, f5, glzf5, grzf5, mz, gamz);               // Z

   Complex_t amp = AmpEEtoZZZ(em, ep, z1, z2, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZZZ()
// --------------------------
Complex_t ZZZBases::AmpEEtoZZZ(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &z1,
                               const HELVector  &z2,
                               const HELVector  &zf)
{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  mz    = fZBosonPtr->GetMass();
   Double_t  gamz  = fZBosonPtr->GetWidth();

   Double_t  me    = kMass[0][1][0];
   Double_t  game  = 0.;

   Double_t  gzzh  = kGz*kM_z;
   Double_t  mh    = 9999.;//125.;
   Double_t  gamh  = 4.e-3;

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
   // ZZZ Production Amplitude
   //---------------------------
#if 1
#if 1
   // (1)
   HELFermion emz1 (em, z1, glze, grze, me, game);
   HELFermion epz2 (ep, z2, glze, grze, me, game);
   HELVertex  amp01(emz1, epz2, zf, glze, grze);

   // (2)
   HELFermion epzf (ep, zf, glze, grze, me, game);
   HELVertex  amp02(emz1, epzf, z2, glze, grze);

   // (3)
   HELFermion emz2 (em, z2, glze, grze, me, game);
   HELFermion epz1 (ep, z1, glze, grze, me, game);
   HELVertex  amp03(emz2, epz1, zf, glze, grze);

   // (4)
   HELVertex  amp04(emz2, epzf, z1, glze, grze);

   // (5)
   HELFermion emzf (em, zf, glze, grze, me, game);
   HELVertex  amp05(emzf, epz1, z2, glze, grze);

   // (6)
   HELVertex  amp06(emzf, epz2, z1, glze, grze);
#else
   // (1)
   HELFermion emz1   (em  , z1, glze, grze, me, game);
   HELFermion emz1zf (emz1, zf, glze, grze, me, game);
   HELVertex  amp01(emz1zf, ep, z2, glze, grze);

   // (2)
   HELFermion emz1z2 (emz1, z2, glze, grze, me, game);
   HELVertex  amp02(emz1z2, ep, zf, glze, grze);

   // (3)
   HELFermion emz2   (em  , z2, glze, grze, me, game);
   HELFermion emz2zf (emz2, zf, glze, grze, me, game);
   HELVertex  amp03(emz2zf, ep, z1, glze, grze);

   // (4)
   HELFermion emz2z1 (emz2, z1, glze, grze, me, game);
   HELVertex  amp04(emz2z1, ep, zf, glze, grze);

   // (5)
   HELFermion emzf   (em  , zf, glze, grze, me, game);
   HELFermion emzfz2 (emzf, z2, glze, grze, me, game);
   HELVertex  amp05(emzfz2, ep, z1, glze, grze);

   // (6)
   HELFermion emzfz1 (emzf, z1, glze, grze, me, game);
   HELVertex  amp06(emzfz1, ep, z2, glze, grze);
#endif
#else
   HELFermion emz1 (em, z1, glze, grze, me, game);
   HELFermion emz2 (em, z2, glze, grze, me, game);
   HELFermion emzf (em, zf, glze, grze, me, game);

   Complex_t  amp01 = AmpEEtoZZ(emz1, ep, z2, zf);
   Complex_t  amp02 = AmpEEtoZZ(emz2, ep, z1, zf);
   Complex_t  amp03 = AmpEEtoZZ(emzf, ep, z1, z2);
   Complex_t  amp04 = 0.;
   Complex_t  amp05 = 0.;
   Complex_t  amp06 = 0.;
#endif

   //--
   // Higgs
   //--
   HELVector  zs (em, ep, glze, grze, mz, gamz);
   // (7)
   HELScalar  h12(z1, z2, gzzh, mh, gamh);
   HELVertex  amp07(zf, zs, h12, gzzh);

   // (8)
   HELScalar  h1f(z1, zf, gzzh, mh, gamh);
   HELVertex  amp08(z2, zs, h1f, gzzh);

   // (9)
   HELScalar  h2f(z2, zf, gzzh, mh, gamh);
   HELVertex  amp09(z1, zs, h2f, gzzh);

   //--
   // Sum
   //--
   Complex_t amp  = amp01 + amp02 + amp03 + amp04 + amp05 + amp06
                  + amp07 + amp08 + amp09;

#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZZ()
// --------------------------
Complex_t ZZZBases::AmpEEtoZZ(const HELFermion &em,
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
void ZZZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZZZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZZZBases.BeamstrahlungFilename","trc500"));
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
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4,   100., fEcmInit,     50, "Mzz"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosZ1"  );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiZ1"  );
  Xh_init( 7,    70.,   110.,       50, "Mz1"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosFd"  );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiFd"  );
  Xh_init(10,    70.,   110.,       50, "Mz2"    );
  Xh_init(11, fXL[6], fXU[6],       50, "CosFdb" );
  Xh_init(12, fXL[7], fXU[7],       50, "PhiFdb" );
  Xh_init(13,    70.,   110.,       50, "Mz"     );
  Xh_init(14, fXL[8], fXU[8],       50, "CosF"   );
  Xh_init(15, fXL[9], fXU[9],       50, "PhiF"   );
  Xh_init(16,     0.,     2.,        2, "Helin " );
  Xh_init(17,     0.,     8.,        8, "Helot " );
  Xh_init(18,     0.,    12.,       12, "Z1 mode");
  Xh_init(19,     0.,    12.,       12, "Z2 mode");
  Xh_init(20,     0.,    12.,       12, "Z  mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZZZBases::Userout()
{
  cout << "End of ZZZBases----------------------------------- "  << endl
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
void ZZZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 8;
   static const Int_t kFHelComb[kNf][6] = {{-1, +1, -1, +1, -1, +1},
                                           {-1, +1, -1, +1, +1, -1},
                                           {-1, +1, +1, -1, -1, +1},
                                           {-1, +1, +1, -1, +1, -1},
                                           {+1, -1, -1, +1, -1, +1},
                                           {+1, -1, -1, +1, +1, -1},
                                           {+1, -1, +1, -1, -1, +1},
                                           {+1, -1, +1, -1, +1, -1}};
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
   fHelFinal  [5] = kFHelComb[fJCombF][5];
   weight = kNf;
}
