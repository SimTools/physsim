//*****************************************************************************
//* =====================
//*  WHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> WH generator
//*
//* (Update Record)
//*    2014/01/30  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "WHSpring.h"
#include "XBoson.h"

#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

#ifdef __NODECAY__
#define __ZEROWIDTH__
#endif
//#define TEMP_H
#ifdef TEMP_H
static TH1F *hMh       = 0;
static TH1F *hRSH      = 0;
static TH1F *hEsum     = 0;
static TH1F *hBSweight = 0;
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(WHSpring)
ClassImp(WHSpringBuf)
ClassImp(WHBases)

//-----------------------------------------------------------------------------
// ==============================
//  class WHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WHSpring::WHSpring(const char      *name,
                   const char      *title,
                           WHBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new WHSpringBuf("WHSpringBuf",
                              "WHSpring event buffer",
                              this);
  if (!bases) { 
    SetBases(new WHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
WHSpring::~WHSpring()
{
  //delete fEventBuf;   // JSFSpring takes care of deleting these
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t WHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    WHBases *bs = static_cast<WHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> WHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class WHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t WHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  WHBases      *bases   = (WHBases*)((WHSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 10;
  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[4] = bases->fP[0]; // u-f  from W1
  pv[5] = bases->fP[1]; // d-fb from W1
  pv[6] = bases->fP[2]; // u-fb from W2
  pv[7] = bases->fP[3]; // d-f  from W2
  pv[8] = bases->fP[4]; // f    from Z
  pv[9] = bases->fP[5]; // fb   from Z

  pv[1] = pv[4] + pv[5]; // W1
  pv[2] = pv[6] + pv[7]; // W2 from H
  pv[3] = pv[8] + pv[9]; // Z  from H
  pv[0] = pv[2] + pv[3]; // H

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2WH         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2W1         = bases->GetQ2W1();
  fCosThetaW1F  = bases->GetCosThetaW1F();
  fPhiW1F       = bases->GetPhiW1F();
  fQ2H          = bases->GetQ2H();
  fCosThetaW2   = bases->GetCosThetaW2();
  fPhiW2        = bases->GetPhiW2();
  fQ2W2         = bases->GetQ2W2();
  fCosThetaW2F  = bases->GetCosThetaW2F();
  fPhiW2F       = bases->GetPhiW2F();
  fQ2Z          = bases->GetQ2H();
  fCosThetaZF   = bases->GetCosThetaZF();
  fPhiZF        = bases->GetPhiZF();
  fCP           = bases->GetCP();

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
  Int_t    cp      = bases->GetCP();               // (W+H-,W-H+)=(+1,-1)
  Int_t    idh     = -37*cp;                       // PDG code for H-
  Int_t    idw1    = 24*cp;                        // PDG code for W+
  Int_t    idw2    = -24*cp;                       // PDG code for W-
  Int_t    idz     = 23;                           // PDG code for Z

  Double_t rq2h    = pv[0].Mag();

  Int_t    idf1    = bases->f1Ptr->GetPID   ()*cp; // PDG code for u-f from W1
  Double_t chrg1   = bases->f1Ptr->GetCharge()*cp; // charge
  Double_t m1      = bases->f1Ptr->GetMass  ();    // mass
  Int_t    idf2    = bases->f2Ptr->GetPID   ()*cp; // PDG code for d-fb from W1
  Double_t chrg2   = bases->f2Ptr->GetCharge()*cp; // charge
  Double_t m2      = bases->f2Ptr->GetMass  ();    // mass
  Int_t    hel1    = bases->fHelFinal[0];          // f1 helicity
  Int_t    hel2    = bases->fHelFinal[1];          // f2 helicity
  Double_t color1  = bases->f1Ptr->GetColor();     // color factor
  Int_t    islev1  = color1 > 1. ? 201 : 0;  	   // shower level
  Int_t    icf1    = 2;                            // color flux id
  Double_t rq2w1   = pv[1].Mag();

  Int_t    idf3    = -bases->f3Ptr->GetPID   ()*cp; // PDG code for u-fb from W2
  Double_t chrg3   = -bases->f3Ptr->GetCharge()*cp; // charge
  Double_t m3      = bases->f3Ptr->GetMass  ();    // mass
  Int_t    idf4    = -bases->f4Ptr->GetPID   ()*cp; // PDG code for d-f from W2
  Double_t chrg4   = -bases->f4Ptr->GetCharge()*cp; // charge
  Double_t m4      = bases->f4Ptr->GetMass  ();    // mass
  Int_t    hel3    = bases->fHelFinal[2];          // f3 helicity
  Int_t    hel4    = bases->fHelFinal[3];          // f4 helicity
  Double_t color2  = bases->f3Ptr->GetColor();     // color factor
  Int_t    islev2  = color2 > 1. ? 301 : 0;  	   // shower level
  Int_t    icf2    = 3;                            // color flux id
  Double_t rq2w2   = pv[2].Mag();

  Int_t    idf5    = bases->f5Ptr->GetPID   ()*cp; // PDG code for f from Z
  Double_t chrg5   = bases->f5Ptr->GetCharge()*cp; // charge
  Double_t m5      = bases->f5Ptr->GetMass  ();    // mass
  Int_t    idf6    = bases->f6Ptr->GetPID   ()*cp; // PDG code for fb from Z
  Double_t chrg6   = bases->f6Ptr->GetCharge()*cp; // charge
  Double_t m6      = bases->f6Ptr->GetMass  ();    // mass
  Int_t    hel5    = bases->fHelFinal[4];          // f5 helicity
  Int_t    hel6    = bases->fHelFinal[5];          // f6 helicity
  Double_t color3  = bases->f5Ptr->GetColor();     // color factor
  Int_t    islev3  = color3 > 1. ? 401 : 0;  	   // shower level
  Int_t    icf3    = 4;                            // color flux id
  Double_t rq2z    = pv[3].Mag();

  Double_t mass    = bases->GetMass();
#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idh << " m=" << mass  << " Q=" << 0.
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 2 pid=" << idw << " m=" << rq2w << " Q=" << 0.
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 3 pid=" << idf << " m=" << m3 << " Q=" << 0.
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 4 pid=" << idf << " m=" << m4 << " Q=" << 0.
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  TVector qcm(4), qh(4), qw(4);
  qh  = *qp[0];
  qw  = *qp[2] + *qp[3];
  qcm = qh + qw;
  cerr << " ph=(" << qh[0] << ", "
                  << qh[1] << ", "
                  << qh[2] << ", "
                  << qh[3] << ") " << endl;
  cerr << " pw=(" << qw[0] << ", "
                  << qw[1] << ", "
                  << qw[2] << ", "
                  << qw[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                                No. PID  Mass   Charge   pv    Nd 1st  Mom hel  col shower
  new (partons[0]) JSFSpringParton( 1, idh ,  rq2h,    -cp, *qp[0], 2, 5,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton( 2, idw1, rq2w1,     cp, *qp[1], 2, 3,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton( 3, idf1,    m1,  chrg1, *qp[4], 0, 0,  2, hel1, icf1, islev1);
  new (partons[3]) JSFSpringParton( 4, idf2,    m2,  chrg2, *qp[5], 0, 0,  2, hel2, icf1, islev1);
  new (partons[4]) JSFSpringParton( 5, idw2, rq2w2,    -cp, *qp[2], 2, 7,  1,    0,    0,      0);
  new (partons[5]) JSFSpringParton( 6, idz , rq2z ,      0, *qp[3], 2, 9,  1,    0,    0,      0);
  new (partons[6]) JSFSpringParton( 7, idf3,    m3,  chrg3, *qp[6], 0, 0,  5, hel3, icf2, islev2);
  new (partons[7]) JSFSpringParton( 8, idf4,    m4,  chrg4, *qp[7], 0, 0,  5, hel4, icf2, islev2);
  new (partons[8]) JSFSpringParton( 9, idf5,    m5,  chrg5, *qp[8], 0, 0,  6, hel5, icf3, islev3);
  new (partons[9]) JSFSpringParton(10, idf6,    m6,  chrg6, *qp[9], 0, 0,  6, hel6, icf3, islev3);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class WHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
WHBases::WHBases(const char *name, const char *title)
         : JSFBases    (name, title), 
           fMass       ( 150.),
           fFhwz       ( 1.),
           fFhwa       ( 0.),
           fEcmInit    (1000.),
           fISR        ( 1),
           fBeamStr    ( 1),
           fBeamWidth  (0.002),
           fPole       (0.),
           fPolp       (0.),
	   fFixCP      ( 1),
           fW1ModesLo  ( 1),
           fW1ModesHi  (12),
           fW2ModesLo  ( 1),
           fW2ModesHi  (12),
           fZModesLo   ( 1),
           fZModesHi   (12),
           fW1BosonPtr ( 0),
           fW2BosonPtr ( 0),
           fZBosonPtr  ( 0),
           fZBoost     (0.),
           fEcmIP      (fEcmInit),
           fQ2WH       (0.),
           fQ2W1       (0.),
           fQ2W2       (0.),
           fQ2Z        (0.),
	   fCP         (1),
           fW1ModePtr  (0),
           fW2ModePtr  (0),
           fZModePtr   (0),
           f1Ptr       (0),
           f2Ptr       (0),
           f3Ptr       (0),
           f4Ptr       (0),
           f5Ptr       (0),
           f6Ptr       (0),
           fCosTheta   (0.),
           fPhi        (0.),
           fXQ2W1      (0.),
           fCosThetaW1F(0.),
           fPhiW1F     (0.),
           fXQ2H       (0.),
           fCosThetaW2 (0.),
           fPhiW2      (0.),
           fXQ2W2      (0.),
           fCosThetaW2F(0.),
           fPhiW2F     (0.),
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

  cout << "Init whbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("WHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.CosthW1FRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiW1FOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.CosthW2Range","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiW2OverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.CosthW2FRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiW2FOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.CosthZFRange","-1.0 1.0"));
  ins >> fXL[8] >> fXU[8];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiZFOverPiRange","0.0 2.0"));
  ins >> fXL[9] >> fXU[9];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.MassH","150.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Fhwz","1.")); 	       // HWZ form factor
  ins >> fFhwz;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Fhwa","0.")); 	       // HWA form factor
  ins >> fFhwa;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.BeamWidth","0.002")); // Beam energy spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Polp","0."));         // positron polarization
  ins >> fPolp;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.FixCP","1"));         // CP combination (W+H-,W-H+)=(1,-1)
  ins >> fFixCP;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.W1ModesLo","1"));      // W1 decay mode lo
  ins >> fW1ModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.W1ModesHi","12"));     // W1 decay mode hi
  ins >> fW1ModesHi;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.W2ModesLo","1"));      // W2 decay mode lo
  ins >> fW2ModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.W2ModesHi","12"));     // W2 decay mode hi
  ins >> fW2ModesHi;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.ZModesLo","1"));       // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.ZModesHi","12"));      // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fCPFinal       , 0., 1., 0, 1);
  DefineVariable(fW1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fW2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2H          , 0., 1., 0, 1);
  DefineVariable(fXQ2W1         , 0., 1., 0, 1);
  DefineVariable(fXQ2W2         , 0., 1., 1, 1);
  DefineVariable(fXQ2Z          , 0., 1., 1, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();
  fXL[7] = fXL[7]*TMath::Pi();
  fXU[7] = fXU[7]*TMath::Pi();
  fXL[9] = fXL[9]*TMath::Pi();
  fXU[9] = fXU[9]*TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaW1F, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiW1F     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaW2 , fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiW2      , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaW2F, fXL[6], fXU[6], 0, 1);
  DefineVariable(fPhiW2F     , fXL[7], fXU[7], 0, 1);
  DefineVariable(fCosThetaZF , fXL[8], fXU[8], 0, 1);
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
  SetNoOfSample(80000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 10);
  SetIteration2(0.05, 20);

}
// --------------------------
//  D-tor
// --------------------------
WHBases::~WHBases()
{
  delete fW1BosonPtr;
  delete fW2BosonPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t WHBases::Func()
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
  //  Select final state CP
  // --------------------------------------------
  Double_t weight = 1.;
  if (fFixCP > 0) {
     fCP = +1;
  } else if (fFixCP < 0) {
     fCP = -1;
  } else {
     if (fCPFinal < 0.5) {
        fCP = +1;
     } else {
        fCP = -1;
     }
     weight = 2.;
  }
  bsWeight *= weight;
  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  GENDecayMode *fW1ModePtr = fW1BosonPtr->PickMode(fW1DecayMode, weight, fW1Mode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fW2ModePtr = fW2BosonPtr->PickMode(fW2DecayMode, weight, fW2Mode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fW2ModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fW2ModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f5Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f6Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m5   = f5Ptr->GetMass();
  Double_t m6   = f6Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  const Double_t kEps = 0.1;
  if (fEcmIP < m1 + m2 + m3 + m4 + m5 + m6 + 6*kEps) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average is taken
  //  care of here
  SelectHelicities(weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  fQ2WH = fEcmIP*fEcmIP;

  Double_t rs   = fEcmIP;
  Double_t qmin = m1 + m2 + 2*kEps;
  Double_t qmax = rs - (m3 + m4 + m5 + m6 + 4*kEps);
  if (qmin >= qmax) return 0.;
#ifndef __ZEROWIDTH__
  fQ2W1 = fW1BosonPtr->GetQ2BW(qmin, qmax, fXQ2W1, weight);
#else
  fQ2W1 = TMath::Power(fW1BosonPtr->GetMass(),2);
  weight = kPi*fW1BosonPtr->GetMass()*fW1BosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t qw1 = TMath::Sqrt(fQ2W1);

  qmin = m3 + m4 + m5 + m6 + 4*kEps;
  qmax = rs - qw1;
  if (qmin >= qmax) return 0.;
#ifndef __ZEROWIDTH___
  fQ2H = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2H, weight);
#else
  fQ2H = TMath::Power(fXBosonPtr->GetMass(),2);
  weight = kPi*fXBosonPtr->GetMass()*fXBosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t qh = TMath::Sqrt(fQ2H);

  qmin = m3 + m4 + 2*kEps;
  qmax = qh - (m5 + m6 + 2*kEps);
  if (qmin >= qmax) return 0.;
#ifndef __ZEROWIDTH__
  fQ2W2 = fW2BosonPtr->GetQ2BW(qmin, qmax, fXQ2W2, weight);
#else
  fQ2W2 = TMath::Power(fW2BosonPtr->GetMass(),2);
  weight = kPi*fW2BosonPtr->GetMass()*fW2BosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t qw2 = TMath::Sqrt(fQ2W2);

  qmin = m5 + m6 + 2*kEps;
  qmax = qh - qw2;
  if (qmin >= qmax) return 0.;
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

  GENBranch w1branch (fQ2W1, fCosThetaW1F, fPhiW1F, m1*m1, m2*m2);
  GENBranch w2branch (fQ2W2, fCosThetaW2F, fPhiW2F, m3*m3, m4*m4);
  GENBranch zbranch  (fQ2Z , fCosThetaZF , fPhiZF , m5*m5, m6*m6);

  GENBranch hbranch  (fQ2H , fCosThetaW2 , fPhiW2 , &w2branch, &zbranch );
  GENBranch cmbranch (fQ2WH, fCosTheta   , fPhi   , &hbranch , &w1branch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

#ifdef TEMP_H
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);
  ANL4DVector qw1 = fP[0] + fP[1];
  ANL4DVector qw2 = fP[2] + fP[3];
  ANL4DVector qz  = fP[4] + fP[5];
  qw.Boost(boostv);
  ANL4DVector qcm(fEcmInit,0.,0.,0.);
  ANL4DVector qmm = qcm - qw1;

  hMh  ->Fill(qmm.Mag()   , (bsWeight*sigma));
  hRSH ->Fill(fEcmIP      , (bsWeight*sigma));
  hEsum->Fill(eplus+eminus, (bsWeight*sigma));
  hBSweight->Fill(eplus+eminus, (bsWeight));
#endif
  Xh_fill( 1, fEcmIP             , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta          , (bsWeight*sigma));
  Xh_fill( 3, fPhi               , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2W1) , (bsWeight*sigma));
  Xh_fill( 5, fCosThetaW1F       , (bsWeight*sigma));
  Xh_fill( 6, fPhiW1F            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2H)  , (bsWeight*sigma));
  Xh_fill( 8, fCosThetaW2        , (bsWeight*sigma));
  Xh_fill( 9, fPhiW2             , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2W2) , (bsWeight*sigma));
  Xh_fill(11, fCosThetaW2F       , (bsWeight*sigma));
  Xh_fill(12, fPhiW2F            , (bsWeight*sigma));
  Xh_fill(13, TMath::Sqrt(fQ2Z)  , (bsWeight*sigma));
  Xh_fill(14, fCosThetaZF        , (bsWeight*sigma));
  Xh_fill(15, fPhiZF             , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fJCombI  , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombF  , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fW1Mode  , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fW2Mode  , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fZMode   , (bsWeight*sigma));
  Xh_fill(21, (Double_t)(fCP+1)/2, (bsWeight*sigma));
  Dh_fill(22, TMath::Sqrt(fQ2W2) ,
              TMath::Sqrt(fQ2Z)  , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t WHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2x  = cmbranch.GetM12();
  Double_t q2w  = cmbranch.GetM22();
  Double_t cosx = cmbranch.GetCosTheta();
  Double_t phix = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x, q2w, cmframe, cosx, phix, 0);
  ANL4DVector px  = phaseCM.GetFourMomentum(0);
  ANL4DVector pw1 = phaseCM.GetFourMomentum(1);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  GENBranch &w1branch = *cmbranch.GetBranchPtr(1);
  Double_t cosw1f = w1branch.GetCosTheta();
  Double_t phiw1f = w1branch.GetPhi     ();
  Double_t m12    = w1branch.GetM12();
  Double_t m22    = w1branch.GetM22();
  GENPhase2 phaseW1(pw1, m12, m22, cmframe, cosw1f, phiw1f, 1);
  fP[0] = phaseW1.GetFourMomentum(0);
  fP[1] = phaseW1.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betaw1f = phaseW1.GetBetaBar();
  if (betaw1f <= 0.) return 0.;

  GENBranch &hbranch = *cmbranch.GetBranchPtr(0);
  Double_t cosw2 = hbranch.GetCosTheta();
  Double_t phiw2 = hbranch.GetPhi     ();
  Double_t mw22  = hbranch.GetM12();
  Double_t mz2   = hbranch.GetM22();
  GENPhase2 phaseH(px, mw22, mz2, cmframe, cosw2, phiw2, 1);
  ANL4DVector pw2 = phaseH.GetFourMomentum(0);
  ANL4DVector pz  = phaseH.GetFourMomentum(1);
  Double_t betaw2 = phaseH.GetBetaBar();
  if (betaw2 <= 0.) return 0.;

  GENBranch &w2branch = *hbranch.GetBranchPtr(0);
  Double_t cosw2f = w2branch.GetCosTheta();
  Double_t phiw2f = w2branch.GetPhi     ();
  Double_t m32    = w2branch.GetM12();
  Double_t m42    = w2branch.GetM22();
  GENFrame hframe = phaseH.GetFrame();
  GENPhase2 phaseW2(pw2, m32, m42, hframe, cosw2f, phiw2f, 1);
  fP[2] = phaseW2.GetFourMomentum(0);
  fP[3] = phaseW2.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betaw2f = phaseW2.GetBetaBar();
  if (betaw2f <= 0.) return 0.;

  GENBranch &zbranch = *hbranch.GetBranchPtr(1);
  Double_t coszf = zbranch.GetCosTheta();
  Double_t phizf = zbranch.GetPhi     ();
  Double_t m52   = zbranch.GetM12();
  Double_t m62   = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m52, m62, hframe, coszf, phizf, 1);
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
  cerr << " fP[4] = (" 
       << fP[4].E () << ","
       << fP[4].Px() << ","
       << fP[4].Py() << ","
       << fP[4].Pz() << ")" << endl;
  cerr << " fP[5] = (" 
       << fP[5].E () << ","
       << fP[5].Px() << ","
       << fP[5].Py() << ","
       << fP[5].Pz() << ")" << endl;
  ANL4DVector qw1 = fP[0] + fP[1];
  cerr << " qw1= (" 
       << qw1.E () << ","
       << qw1.Px() << ","
       << qw1.Py() << ","
       << qw1.Pz() << ")" << endl;
  ANL4DVector qw2 = fP[2] + fP[3];
  cerr << " qw2= (" 
       << qw2.E () << ","
       << qw2.Px() << ","
       << qw2.Py() << ","
       << qw2.Pz() << ")" << endl;
  ANL4DVector qz = fP[4] + fP[5];
  cerr << " qz= (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector qh = qw2 + qz;
  ANL4DVector pcm = qh + qw1;
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

  Double_t identp = 1.;                                                  // identical particle factor
  Double_t dPhase = kFact * betax * betaw1f * betaw2 * betaw2f * betazf; // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);                                 // beam flux factor

  Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                       // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t WHBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor() * f5Ptr->GetColor();
  Int_t     ig11  = f1Ptr->GetGenNo() - 1;
  Int_t     ig12  = f2Ptr->GetGenNo() - 1;
  Int_t     ig21  = f3Ptr->GetGenNo() - 1;
  Int_t     ig22  = f4Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig11][ig12]*kVkm[ig21][ig22],2);
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t WHBases::FullAmplitude()
{
   Double_t gamw   = fW1BosonPtr->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0;

   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf     = f5Ptr->GetCharge();
   Double_t t3f    = f5Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   Double_t mx     = fXBosonPtr->GetMass();
   Double_t gamx   = fXBosonPtr->GetWidth();
   Double_t gzwh   = kGw*kM_w*fFhwz;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   Int_t i1 = 0;
   Int_t i2 = 1;
   Int_t i3 = 2;
   Int_t i4 = 3;
   Int_t i5 = 4;
   Int_t i6 = 5;
   if (fCP < 0) {
      i1 = 1;
      i2 = 0;
      i3 = 3;
      i4 = 2;
      i5 = 5;
      i6 = 4;
   }
   HELFermion f1 (fP[i1], fM[i1], fHelFinal [i1], +1, kIsOutgoing);
   HELFermion f2b(fP[i2], fM[i2], fHelFinal [i2], -1, kIsIncoming);
   HELVector  w1f(f2b, f1, glw, grw, kM_w, gamw);

   HELFermion f3b(fP[i3], fM[i3], fHelFinal [i3], -1, kIsIncoming);
   HELFermion f4 (fP[i4], fM[i4], fHelFinal [i4], +1, kIsOutgoing);
   HELVector  w2f(f3b, f4, glw, grw, kM_w, gamw);

   HELFermion f5 (fP[i5], fM[i5], fHelFinal [i5], +1, kIsOutgoing);
   HELFermion f6b(fP[i6], fM[i6], fHelFinal [i6], -1, kIsIncoming);
   HELVector  zf (f6b, f5, glz, grz, kM_z, gamz);

   HELScalar  xf(w2f, zf, gzwh, mx, gamx);

   Complex_t amp = AmpEEtoWH(em, ep, xf, w1f);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoWH()
// --------------------------
Complex_t WHBases::AmpEEtoWH(const HELFermion &em,
                             const HELFermion &ep,
                             const HELScalar  &xf,
                             const HELVector  &wf)
{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);
   Double_t  glae  = -kGe*qe;
   Double_t  grae  = glae;

   Double_t  gamz  = fZBosonPtr->GetWidth();

   Double_t  ma    = 0.; // photon mass
   Double_t  gama  = 0.; // photon width

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
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Double_t gzwh   = kGw*kM_w*fFhwz;

   Complex_t amp = HELVertex(zs, wf, xf, gzwh);

   HELVector as(em, ep, glae, grae, ma, gama);
   Double_t gawh   = kGw*kM_w*fFhwa;
             amp += HELVertex(as, wf, xf, gawh);
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void WHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("WHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("WHBases.BeamstrahlungFilename","trc500"));
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
  //  Initialize X decay table
  // --------------------------------------------
  if (!fXBosonPtr) fXBosonPtr = new XBoson(fMass,fFhwz,fFhwa);
  fXBosonPtr->DebugPrint();
  
  // --------------------------------------------
  //  Initialize W decay table
  // --------------------------------------------
  if (!fW1BosonPtr) fW1BosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fW1BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fW1BosonPtr->GetMode(m); 
     if (mp && (m<fW1ModesLo || m>fW1ModesHi)) {
        mp->Lock();
     }
  }
  fW1BosonPtr->DebugPrint();

  if (!fW2BosonPtr) fW2BosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fW2BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fW2BosonPtr->GetMode(m); 
     if (mp && (m<fW2ModesLo || m>fW2ModesHi)) {
        mp->Lock();
     }
  }
  fW2BosonPtr->DebugPrint();

  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
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
  Xh_init( 2, fXL[0], fXU[0],       50, "CosthH" );
  Xh_init( 3, fXL[1], fXU[1],       50, "PhiH"   );
  Xh_init( 4,    50.,   100.,       50, "Mw"     );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthF" );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF"   );
  Xh_init( 7, fMass-10., fMass+10., 50, "Mh"     );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthW" );
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiW"   );
  Xh_init(10,    10.,   110.,       50, "Mw2"    );
  Xh_init(11, fXL[6], fXU[6],       50, "CosthF" );
  Xh_init(12, fXL[7], fXU[7],       50, "PhiF"   );
  Xh_init(13,    10.,   110.,       50, "Mz"     );
  Xh_init(14, fXL[8], fXU[8],       50, "CosthF" );
  Xh_init(15, fXL[9], fXU[9],       50, "PhiF"   );
  Xh_init(16,     0.,     2.,        2, "Helin  ");
  Xh_init(17,     0.,     2.,        2, "Helot  ");
  Xh_init(18,     0.,    12.,       12, "W1 mode");
  Xh_init(19,     0.,    12.,       12, "W2 mode");
  Xh_init(20,     0.,    12.,       12, "Z  mode");
  Xh_init(21,     0.,     2.,        2, "CP     ");
  Dh_init(22,    10.,   110.,       50,
                 10.,   110.,       50, "(Mw,Mz)");
#ifdef TEMP_H
  if (!hMh)   hMh   = new TH1F("hMh"  ,"", 200,115.,135.);
  if (!hRSH)  hRSH  = new TH1F("hRSH" ,"",1100,  0.,fEcmInit*1.1);
  if (!hEsum) hEsum = new TH1F("hEsum","",1100,  0.,fEcmInit*1.1);
  if (!hBSweight) hBSweight = new TH1F("hBSweight","",1100,0,fEcmInit*1.1);
#endif
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void WHBases::Userout()
{
  cout << "End of WHBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremsstrahlung       = " << (fISR     ? "on" : "off")  << endl
       << "e- Polarization      = " << GetPole()                  << endl
       << "e+ Polarization      = " << GetPolp()                  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
#ifdef TEMP_H
  hMh  ->Write();
  hRSH ->Write();
  hEsum->Write();
  hBSweight->Write();
#endif
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void WHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][6] = {{-1, +1, +1, -1, -1, +1},
                                           {-1, +1, +1, -1, +1, -1}};
   Double_t helm = (1. - fPole)/2.;
   if (fHelCombInitial < helm) {
      fJCombI = 0;
      weight  = (1. + fPolp)/2.;
   } else {
      fJCombI = 1;
      weight  = (1. - fPolp)/2.;
   }
   fHelInitial[0] = kIHelComb[fJCombI][0];
   fHelInitial[1] = kIHelComb[fJCombI][1];
   fJCombF = (Int_t)(fHelCombFinal*kNf);
   fJCombF = TMath::Min(fJCombF, kNf-1);
   fHelFinal  [0] = kFHelComb[fJCombF][0]*fCP;
   fHelFinal  [1] = kFHelComb[fJCombF][1]*fCP;
   fHelFinal  [2] = kFHelComb[fJCombF][2]*fCP;
   fHelFinal  [3] = kFHelComb[fJCombF][3]*fCP;
   fHelFinal  [4] = kFHelComb[fJCombF][4]*fCP;
   fHelFinal  [5] = kFHelComb[fJCombF][5]*fCP;
   weight *= kNf;
}
