//*****************************************************************************
//* =====================
//*  WWZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> WWZ generator
//*
//* (Update Record)
//*    2010/04/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "WWZSpring.h"
#include "HBoson.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(WWZSpring)
ClassImp(WWZSpringBuf)
ClassImp(WWZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class WWZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WWZSpring::WWZSpring(const char *name,
                     const char *title,
                     WWZBases   *bases)
        : JSFSpring(name, title, bases)
{
  fEventBuf = new WWZSpringBuf("WWZSpringBuf",
                               "WWZSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new WWZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
WWZSpring::~WWZSpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t WWZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    WWZBases *bs = static_cast<WWZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> WWZBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class WWZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t WWZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  WWZBases     *bases   = static_cast<WWZBases *>(
                          static_cast<WWZSpring *>(Module())->GetBases());

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
  fQ2WW         = bases->GetQ2WW();
  fCosThetaWm   = bases->GetCosThetaWm();
  fPhiWm        = bases->GetPhiWm();
  fQ2Wm         = bases->GetQ2Wm();
  fCosThetaWmF  = bases->GetCosThetaWmF();
  fPhiWmF       = bases->GetPhiWmF();
  fQ2Wp         = bases->GetQ2Wp();
  fCosThetaWpF  = bases->GetCosThetaWpF();
  fPhiWpF       = bases->GetPhiWpF();
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
  Int_t    idw     = 24;                          // PDG code for W
  Int_t    idz     = 23;                          // PDG code for Z
  // W-
  Int_t    idf1    = bases->f1Ptr->GetPID   ();   // PDG code for f1
  Double_t chrg1   = bases->f1Ptr->GetCharge();   // f1 charge
  Double_t m1      = bases->f1Ptr->GetMass  ();   // f1 mass
  Int_t    hel1    = bases->fHelFinal[0];         // f1 helicity
  Double_t color1  = bases->f1Ptr->GetColor();    // color factor for f1

  Int_t    idf2    = bases->f2Ptr->GetPID   ();   // PDG code for f2
  Double_t chrg2   = bases->f2Ptr->GetCharge();   // f2 charge
  Double_t m2      = bases->f2Ptr->GetMass  ();   // f2 mass
  Int_t    hel2    = bases->fHelFinal[1];         // f2 helicity

  Int_t    islevwm = color1 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfwm   = 2;                           // color flux id
  Double_t rq2wm   = pv[6].Mag();

  // W+
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f3 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f4 helicity

  Int_t    islevwp = color3 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfwp   = 3;                           // color flux id
  Double_t rq2wp   = pv[7].Mag();

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
  TVector qcm(4), qwm(4), qwp(4), qz(4);
  qwm = *qp[0] + *qp[1];
  qwp = *qp[2] + *qp[3];
  qz  = *qp[4] + *qp[5];
  qcm = qwm + qwp + qz;

  cerr << " pwm=(" << qwm[0] << ", "
                   << qwm[1] << ", "
                   << qwm[2] << ", "
                   << qwm[3] << ") " << endl;
  cerr << " pwp=(" << qwp[0] << ", "
                   << qwp[1] << ", "
                   << qwp[2] << ", "
                   << qwp[3] << ") " << endl;
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
  new (partons[0]) JSFSpringParton( 1, -idw,rq2wm,   -1., *qp[6], 2, 4,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2,  idw,rq2wp,   +1., *qp[7], 2, 6,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3,  idz, rq2z,    0., *qp[8], 2, 8,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton( 4,-idf1,   m1,-chrg1, *qp[0], 0, 0,  1, hel1,icfwm,islevwm);
  new (partons[4]) JSFSpringParton( 5,-idf2,   m2,-chrg2, *qp[1], 0, 0,  1, hel2,icfwm,islevwm);
  new (partons[5]) JSFSpringParton( 6, idf3,   m3, chrg3, *qp[2], 0, 0,  2, hel3,icfwp,islevwp);
  new (partons[6]) JSFSpringParton( 7, idf4,   m4, chrg4, *qp[3], 0, 0,  2, hel4,icfwp,islevwp);
  new (partons[7]) JSFSpringParton( 8, idf5,   m5, chrg5, *qp[4], 0, 0,  3, hel5, icfz, islevz);
  new (partons[8]) JSFSpringParton( 9, idf6,   m6, chrg6, *qp[5], 0, 0,  3, hel6, icfz, islevz);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class WWZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WWZBases::WWZBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fMass       (9999.),
           fLambda     (1000.),
           fA          (   0.),
           fB          (   0.),
           fBtilde     (   0.),
           fEcmInit    (500.),
           fISR        ( 1),
           fBeamStr    ( 1),
           fBeamWidth  (0.002),
           fPole       (0.),
           fWmModesLo  ( 1),
           fWmModesHi  (12),
           fWpModesLo  ( 1),
           fWpModesHi  (12),
           fZModesLo   ( 1),
           fZModesHi   (12),
	   fNCALL      (80000),
	   fACC1       (0.05),
	   fACC2       (0.05),
	   fITMX1      (20),
	   fITMX2      (40),
           fHBosonPtr  ( 0),
           fWmBosonPtr ( 0),
           fWpBosonPtr ( 0),
           fZBosonPtr  ( 0),
           fZBoost     (0.),
           fEcmIP      (fEcmInit),
           fQ2WW       (0.),
           fQ2Wm       (0.),
           fQ2Wp       (0.),
           fQ2Z        (0.),
           fWmModePtr  (0),
           f1Ptr       (0),
           f2Ptr       (0),
           fWpModePtr  (0),
           f3Ptr       (0),
           f4Ptr       (0),
           fZModePtr   (0),
           f5Ptr       (0),
           f6Ptr       (0),
           fCosTheta   (0.),
           fPhi        (0.),
           fXQ2WW      (0.),
           fCosThetaWm (0.),
           fPhiWm      (0.),
           fXQ2Wm      (0.),
           fCosThetaWmF(0.),
           fPhiWmF     (0.),
           fXQ2Wp      (0.),
           fCosThetaWpF(0.),
           fPhiWpF     (0.),
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
  stringstream ins(gJSF->Env()->GetValue("WWZBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.CosthWmRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.PhiWmOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.CosthWmFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.PhiWmFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.CosthWpFRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.PhiWpFOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.CosthZFRange","-1.0 1.0"));
  ins >> fXL[8] >> fXU[8];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.PhiZFOverPiRange","0.0 2.0"));
  ins >> fXL[9] >> fXU[9];

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.WmModesLo","1"));      // W- decay mode lo
  ins >> fWmModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.WmModesHi","12"));     // W- decay mode hi
  ins >> fWmModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.WpModesLo","1"));      // W+ decay mode lo
  ins >> fWpModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.WpModesHi","12"));     // W+ decay mode hi
  ins >> fWpModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.HiggsMass","9999."));  // m_H [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.Lambda","1000.")); 	 // Lambda [GeV]
  ins >> fLambda;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.A","0.")); 	 // a
  ins >> fA;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.B","0.")); 	 // b
  ins >> fB;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.Btilde","0.")); 	 // btilde
  ins >> fBtilde;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("WWZBases.NCALL","80000"));
  ins >> fNCALL;

  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 1, 1);
  DefineVariable(fWmDecayMode   , 0., 1., 0, 1);
  DefineVariable(fWpDecayMode   , 0., 1., 0, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2WW         , 0., 1., 1, 1);
  DefineVariable(fXQ2Wm         , 0., 1., 1, 1);
  DefineVariable(fXQ2Wp         , 0., 1., 1, 1);
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
  DefineVariable(fCosThetaWm , fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiWm      , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaWmF, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiWmF     , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaWpF, fXL[6], fXU[6], 1, 1);
  DefineVariable(fPhiWpF     , fXL[7], fXU[7], 0, 1);
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
WWZBases::~WWZBases()
{
  delete fHBosonPtr;
  delete fWmBosonPtr;
  delete fWpBosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t WWZBases::Func()
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

  GENDecayMode *fWmModePtr = fWmBosonPtr->PickMode(fWmDecayMode, weight, fWmMode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fWpModePtr = fWpBosonPtr->PickMode(fWpDecayMode, weight, fWpMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(1));
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

  Double_t qwwmin  = m1 + m2 + m3 + m4;
  Double_t qwwmax  = rs - (m5 + m6);
#ifndef HIGGS_ONLY
  Double_t q2wwmin = qwwmin*qwwmin;
  Double_t q2wwmax = qwwmax*qwwmax;
  fQ2WW  = q2wwmin + (q2wwmax - q2wwmin)*fXQ2WW;
  weight = q2wwmax - q2wwmin;
#else
#if 0
  fQ2WW  = fHBosonPtr->GetQ2BW(qwwmin, qwwmax, fXQ2WW, weight);
#else
  // Zero width approx.
  fQ2WW  = TMath::Power(fHBosonPtr->GetMass(),2);
  weight = kPi*fHBosonPtr->GetMass()*fHBosonPtr->GetWidth();
#endif
#endif
  bsWeight *= weight;

  Double_t rq2ww  = TMath::Sqrt(fQ2WW);
  Double_t qzmin  = m5 + m6;
  Double_t qzmax  = rs - rq2ww;
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qzmin, qzmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t qwmmin = m1 + m2;
  Double_t qwmmax = rq2ww - (m3 + m4);
#ifndef __ZEROWIDTH__
  fQ2Wm = fWmBosonPtr->GetQ2BW(qwmmin, qwmmax, fXQ2Wm, weight);
#else
  fQ2Wm = TMath::Power(fWmBosonPtr->GetMass(),2);
  weight = kPi*fWmBosonPtr->GetMass()*fWmBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  Double_t rq2wm  = TMath::Sqrt(fQ2Wm);
  Double_t qwpmin = m3 + m4;
  Double_t qwpmax = rq2ww - rq2wm;
#ifndef __ZEROWIDTH__
  fQ2Wp = fWpBosonPtr->GetQ2BW(qwpmin, qwpmax, fXQ2Wp, weight);
#else
  fQ2Wp = TMath::Power(fWpBosonPtr->GetMass(),2);
  weight = kPi*fWpBosonPtr->GetMass()*fWpBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch wmbranch(fQ2Wm, fCosThetaWmF, fPhiWmF, m1*m1    , m2*m2);
  GENBranch wpbranch(fQ2Wp, fCosThetaWpF, fPhiWpF, m3*m3    , m4*m4);
  GENBranch zbranch (fQ2Z , fCosThetaZF , fPhiZF , m5*m5    , m6*m6);
  GENBranch wwbranch(fQ2WW, fCosThetaWm , fPhiWm , &wmbranch, &wpbranch);
  GENBranch cmbranch(s    , fCosTheta   , fPhi   , &wwbranch, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  ANL4DVector pwm = fP[0] + fP[1];
  ANL4DVector pwp = fP[2] + fP[3];
  ANL4DVector pz  = fP[4] + fP[5];
  ANL4DVector ph  = pwm + pwp;
  ANL4DVector phb(ph.E(), -ph.Px(), -ph.Py(), -ph.Pz());
  TVector3 boostv = phb.BoostVector();
  ANL4DVector qf[4];
  for (Int_t i=0; i<4; i++) {
	  qf[i] = fP[i];
	  qf[i].Boost(boostv);
  }
  ANL4DVector qwm = pwm;
  ANL4DVector qwp = pwp;
  qwm.Boost(boostv);
  qwp.Boost(boostv);
  TVector3 vwm = qwm.Vect().Unit();
  TVector3 v0  = qf[0].Vect().Unit();
  TVector3 v2  = qf[2].Vect().Unit();
  TVector3 vn1 = v0.Cross(vwm);
  TVector3 vn2 = v2.Cross(vwm);
  TVector3 vn3 = vn1.Cross(vn2);
  Double_t sn = vwm*vn3;
  Double_t cs = vn1*vn2;
  Double_t phidiff = TMath::ATan2(sn,cs);
  phidiff = phidiff < 0. ? phidiff + k2Pi : phidiff;
  Double_t aqw = qwm.Vect().Mag();

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2WW), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaWm       , (bsWeight*sigma));
  Xh_fill( 6, fPhiWm            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Wm), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaWmF      , (bsWeight*sigma));
  Xh_fill( 9, fPhiWmF           , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2Wp), (bsWeight*sigma));
  Xh_fill(11, fCosThetaWpF      , (bsWeight*sigma));
  Xh_fill(12, fPhiWpF           , (bsWeight*sigma));
  Xh_fill(13, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill(14, fCosThetaZF       , (bsWeight*sigma));
  Xh_fill(15, fPhiZF            , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fWmMode , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fWpMode , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fZMode  , (bsWeight*sigma));
  Xh_fill(21, phidiff           , (bsWeight*sigma));
  Xh_fill(22, aqw               , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t WWZBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2ww  = cmbranch.GetM12();
  Double_t q2z   = cmbranch.GetM22();
  Double_t cosww = cmbranch.GetCosTheta();
  Double_t phiww = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2ww, q2z, cmframe, cosww, phiww, 0);
  ANL4DVector pww = phaseCM.GetFourMomentum(0);
  ANL4DVector pz  = phaseCM.GetFourMomentum(1);
  Double_t betaww = phaseCM.GetBetaBar();
  if (betaww <= 0.) return 0.;

  GENBranch &wwbranch = *cmbranch.GetBranchPtr(0);
  Double_t q2wm = wwbranch.GetM12();
  Double_t q2wp = wwbranch.GetM22();
  Double_t cosw = wwbranch.GetCosTheta();
  Double_t phiw = wwbranch.GetPhi     ();
  GENPhase2 phaseWW(pww, q2wm, q2wp, cmframe, cosw, phiw, 1);
  ANL4DVector pwm = phaseWW.GetFourMomentum(0);
  ANL4DVector pwp = phaseWW.GetFourMomentum(1);
  Double_t betaw  = phaseWW.GetBetaBar();
  if (betaw <= 0.) return 0.;

  GENBranch &wmbranch = *wwbranch.GetBranchPtr(0);
  Double_t coswmf  = wmbranch.GetCosTheta();
  Double_t phiwmf  = wmbranch.GetPhi     ();
  Double_t m12     = wmbranch.GetM12();
  Double_t m22     = wmbranch.GetM22();
  GENFrame wwframe = phaseWW.GetFrame();
  GENPhase2 phaseWm(pwm, m12, m22, wwframe, coswmf, phiwmf, 1);
  fP[0] = phaseWm.GetFourMomentum(0);
  fP[1] = phaseWm.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betawmf = phaseWm.GetBetaBar();
  if (betawmf <= 0.) return 0.;

  GENBranch &wpbranch = *wwbranch.GetBranchPtr(1);
  Double_t coswpf = wpbranch.GetCosTheta();
  Double_t phiwpf = wpbranch.GetPhi     ();
  Double_t m32    = wpbranch.GetM12();
  Double_t m42    = wpbranch.GetM22();
  GENPhase2 phaseWp(pwp, m32, m42, wwframe, coswpf, phiwpf, 1);
  fP[2] = phaseWp.GetFourMomentum(0);
  fP[3] = phaseWp.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betawpf = phaseWp.GetBetaBar();
  if (betawpf <= 0.) return 0.;

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
  ANL4DVector qwm = fP[0] + fP[1];
  ANL4DVector qwp = fP[2] + fP[3];
  ANL4DVector qz  = fP[4] + fP[5];
  cerr << " qwm = (" 
       << qwm.E () << ","
       << qwm.Px() << ","
       << qwm.Py() << ","
       << qwm.Pz() << ")" << endl;
  cerr << " qwp = (" 
       << qwp.E () << ","
       << qwp.Px() << ","
       << qwp.Py() << ","
       << qwp.Pz() << ")" << endl;
  cerr << " qz  = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector pcm = qwm + qwp + qz;
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

  Double_t identp = 1.;                                // identical particle factor
  Double_t dPhase = kFact * betaww * betaw 
	                  * betawmf * betawpf * betazf;// phase space factor
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
Double_t WWZBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor() * f5Ptr->GetColor();
  Int_t     ig1   = f1Ptr->GetGenNo() - 1;
  Int_t     ig2   = f2Ptr->GetGenNo() - 1;
  Int_t     ig3   = f3Ptr->GetGenNo() - 1;
  Int_t     ig4   = f4Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig1][ig2]*kVkm[ig3][ig4],2);

  Complex_t amp    = FullAmplitude();
  Double_t  amp2   = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t WWZBases::FullAmplitude()
{
   Double_t gamw   = fWmBosonPtr->GetWidth();
   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;

   Double_t gamz   = fZBosonPtr->GetWidth();
   Double_t qf     = f5Ptr->GetCharge();
   Double_t t3f    = f5Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

   HELFermion f1b(fP[0], fM[0], fHelFinal [0], -1, kIsIncoming); // fubar
   HELFermion f2 (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing); // fd
   HELVector  wm(f1b, f2, glw, grw, kM_w, gamw);                 // W-

   HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fu
   HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // fdbar
   HELVector  wp(f4b, f3, glw, grw, kM_w, gamw);                 // W+

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing);
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming);
   HELVector  zf(f6b, f5, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoWWZ(em, ep, wm, wp, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoWWZ()
// --------------------------
Complex_t WWZBases::AmpEEtoWWZ(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &wm,
                               const HELVector  &wp,
                               const HELVector  &zf)
{
   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  t3n   = +1./2.;
   Double_t  glzn  = -kGz*t3n;
   Double_t  grzn  = 0.;

   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  gamw  = fWmBosonPtr->GetWidth();

   Double_t  gwwz   = kGw*kCosW ;
   Double_t  gwwz3  = kGw;
   Double_t  glw   = -kGw*kSqh;
   Double_t  grw   = 0.;

   Double_t  mne   = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamne = 0.;
   Double_t  game  = 0.;

   Double_t  ghzz  = kGz*kM_z;
   Double_t  ghww  = kGw*kM_w;
   Double_t  mh    = fMass;
   Double_t  gamh  = fHBosonPtr->GetWidth();

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
   // WWZ Production Amplitude
   //---------------------------
   //--
   // S-channel
   //--
   // (1-3)
   HELVector w3(em, ep, glae, grae, glze, grze, kM_z, gamz);
   HELVertex ampwww3w3(wm, w3, wp, zf, gwwz3, gwwz, kM_w, gamw, kTRUE);
   // (4)
   HELFermion emv(em, zf, glze, grze, me, game);
   HELVector w3emv(emv, ep, glae, grae, glze, grze, kM_z, gamz);
   HELVertex ampwmwpz3emv(wm, wp, w3emv, gwwz3);
   // (5)
   HELFermion epv(ep, zf, glze, grze, me, game);
   HELVector w3epv(em, epv, glae, grae, glze, grze, kM_z, gamz);
   HELVertex ampwmwpz3epv(wm, wp, w3epv, gwwz3);

   //--
   // T-channel ne
   //--
   // (6)
   HELVector  wmv(zf, wm, gwwz, kM_w, gamw);
   HELFermion newmv(em, wmv, glw, grw, mne, gamne);
   HELVertex  ampteewmv(newmv, ep, wp, glw, grw);
   // (7)
   HELVector  wpv(wp, zf, gwwz, kM_w, gamw);
   HELFermion ne(em, wm, glw, grw, mne, gamne);
   HELVertex  ampteewpv(ne, ep, wpv, glw, grw);
   // (8)
   HELFermion nemv(emv, wm, glw, grw, mne, gamne);
   HELVertex  ampteeemv(nemv, ep, wp, glw, grw);
   // (9)
   HELVertex  ampteeepv(ne, epv, wp, glw, grw);
   // (10)
   HELFermion nev (ne, zf, glzn, grzn, mne, gamne);
   HELVertex  ampteenev(nev, ep, wp, glw, grw);

   //--
   // Higgs
   //--
   // (11)
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
#ifndef ANOM_WWH
   HELScalar hf(wm, wp, ghww, mh, gamh);
   Complex_t amph = HELVertex(zs, zf, hf, ghzz);
#else
   Double_t g1     = ghww + 2 * kM_w * kM_w * (fA/fLambda);
   Double_t g2     = -2 * (fB/fLambda);
   Double_t g3     = -4 * (fBtilde/fLambda);
   HELScalar hf(wm, wp, g1, g2, g3, mh, gamh);

   Double_t g1z    = ghzz + 2 * kM_z * kM_z * (fA/fLambda);
   Double_t g2z    = -2 * (fB/fLambda);
   Double_t g3z    = -4 * (fBtilde/fLambda);
   Complex_t amph = HELVertex(zs, zf, hf, g1z, g2z, g3z);
#endif

   //--
   // Sum
   //--
   Complex_t amp  = ampwww3w3 
                   + ampwmwpz3emv + ampwmwpz3epv
	           + ampteewmv + ampteewpv 
		   + ampteeemv + ampteeepv 
		   + ampteenev
		   + amph
		   ;

#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void WWZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("WWZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("WWZBases.BeamstrahlungFilename","trc500"));
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
  if (!fWmBosonPtr) fWmBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWmBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWmBosonPtr->GetMode(m); 
     if (mp && (m<fWmModesLo || m>fWmModesHi)) {
        mp->Lock();
     }
  }
  fWmBosonPtr->DebugPrint();

  if (!fWpBosonPtr) fWpBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWpBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWpBosonPtr->GetMode(m); 
     if (mp && (m<fWpModesLo || m>fWpModesHi)) {
        mp->Lock();
     }
  }
  fWpBosonPtr->DebugPrint();

  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZBosonPtr->GetMode(m); 
     if (mp && (m<fZModesLo || m>fZModesHi)) {
        mp->Lock();
     }
  }
  fZBosonPtr->DebugPrint();

  if (!fHBosonPtr) fHBosonPtr = new HBoson(fMass,fLambda,fA,fB,fBtilde);
  fHBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4,   100.,   300.,       50, "Mww"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosW-"  );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiW-"  );
  Xh_init( 7,    60.,   100.,       50, "Mw-"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosFd"  );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiFd"  );
  Xh_init(10,    60.,   100.,       50, "Mw+"    );
  Xh_init(11, fXL[6], fXU[6],       50, "CosFdb" );
  Xh_init(12, fXL[7], fXU[7],       50, "PhiFdb" );
  Xh_init(13,    70.,   110.,       50, "Mz"     );
  Xh_init(14, fXL[8], fXU[8],       50, "CosF"   );
  Xh_init(15, fXL[9], fXU[9],       50, "PhiF"   );
  Xh_init(16,     0.,     2.,        2, "Helin " );
  Xh_init(17,     0.,     1.,        2, "Helot " );
  Xh_init(18,     0.,    12.,       12, "W- mode");
  Xh_init(19,     0.,    12.,       12, "W+ mode");
  Xh_init(20,     0.,    12.,       12, "Z  mode");
  Xh_init(21, fXL[5], fXU[5],       50, "DPhi"   );
  Xh_init(22,     0.,   150.,       50, "|qw|"   );
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void WWZBases::Userout()
{
  cout << "End of WWZBases----------------------------------- "  << endl
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
void WWZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][6] = {{+1, -1, -1, +1, -1, +1},
                                           {+1, -1, -1, +1, +1, -1}};
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
