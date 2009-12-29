//*****************************************************************************
//* =====================
//*  NnSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> Nn  generator
//*
//* (Update Record)
//*    2009/08/11  T.Saito	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "NnSpring.h"

#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __NOWDECAY__
//#define __DEBUG__
//#define __ZEROWIDTH__
#ifdef __NOWDECAY__
#ifndef __ZEROWIDTH__
#define __ZEROWIDTH__
#endif
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"
//                                       N_e, N_mu, N_tau
static const Double_t kMNS [3][3] = {{  0.83,  0.56, -0.07},  // e
                                     { -0.35,  0.61,  0.71},  // mu
                                     {  0.44, -0.56,  0.70}}; // tau
static const Double_t kMnu [3] = {0., 9.0e-12, 5.9e-11};

static TH1D *hEcm   = 0;
static TH1D *hCosNR = 0;
static TH1D *hPhiNR = 0;
static TH1D *hMNR   = 0;
static TH1D *hCosW  = 0;
static TH1D *hPhiW  = 0;
static TH1D *hCosle = 0;
static TH1D *hPhile = 0;
static TH1D *hMw    = 0;
static TH1D *hCosU  = 0;
static TH1D *hPhiU  = 0;
static TH1D *hWDK   = 0;
static TH1D *hHelIn = 0;
static TH1D *hHelOt = 0;

static Double_t g_CosThetaLepton = 100.;
static Double_t g_PhiLepton      = 100.;

ClassImp(NnSpring)
ClassImp(NnSpringBuf)
ClassImp(NnBases)
ClassImp(RNeutrino)


//-----------------------------------------------------------------------------
// ==============================
//  class NnSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NnSpring::NnSpring(const char      *name,
                   const char      *title,
                         NnBases   *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new NnSpringBuf("NnSpringBuf",
                              "NnSpring event buffer",
                              this);

  if (!bases) { 
    SetBases(new NnBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
NnSpring::~NnSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class NnSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t NnSpringBuf::SetPartons()
{
  std::cout << "NnSprintBuf::SetPartons()" << std::endl;
  TClonesArray &partons = *fPartons;
  NnBases      *bases   = (NnBases*)((NnSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 6;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // neutrino
  pv[3] = bases->fP[1]; // charged lepton
  pv[4] = bases->fP[2]; // up quark
  pv[5] = bases->fP[3]; // down quark
  pv[2] = pv[4] + pv[5]; // W
  pv[1] = pv[2] + pv[3]; // right handed neutrino


  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2XX         = fEcmIP*fEcmIP;
  fQ2X1         = bases->GetQ2X1();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fCosThetaW1   = bases->GetCosThetaW1();
  fPhiW1        = bases->GetPhiW1();
  fQ2W1         = bases->GetQ2W1();
  fCosThetaF1   = bases->GetCosThetaF1();
  fPhiF1        = bases->GetPhiF1();

  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);

  for (Int_t i=0; i<fNparton; i++) pv[i].Boost(boostv);

  // ----------------------------------------------
  //  Set final state parton infomation
  // ----------------------------------------------
  Int_t cp = bases->fHelCombInitial < 0.5 ? +1 : -1; // cp flip factor

  static TVector **qp = 0;
  if (!qp) {
    qp = new TVector* [fNparton];
    for (Int_t i=0; i<fNparton; i++) qp[i] = new TVector(4);
  }
  for (Int_t i=0; i<fNparton; i++) {
    TVector &q = *qp[i];
    q(0) =    pv[i].E ();
    q(1) = cp*pv[i].Px();
    q(2) = cp*pv[i].Py();
    q(3) = cp*pv[i].Pz();
  }

  Int_t    gennu   = bases->GetGenNu();
  Int_t    idnu    = -cp*(gennu == 1 ? 12 : (gennu == 2 ? 14 : 16));   // PDG code for f
  Double_t mnu     = kMnu[gennu-1];       // nu mass

  Int_t    idx     = cp*bases->fNRPtr->GetPID();  // N PDG : meaningless!
  Int_t    idl     = cp*(2*bases->fGenLepton + 9);

  Double_t ml      = bases->fFPtr->GetMass();

  Int_t    helnu   = cp*bases->fHelFinal[0];
  Int_t    hell    = cp*bases->fHelFinal[1];

  Double_t chgl    = -cp;
  Double_t chgw    =  cp;

  Int_t    idu     = cp*bases->f1Ptr->GetPID();
  Int_t    idd     = cp*bases->f2Ptr->GetPID();
  Double_t mu      = bases->f1Ptr->GetMass();
  Double_t md      = bases->f2Ptr->GetMass();
  Int_t    helu    = cp*bases->fHelFinal[2];
  Int_t    held    = cp*bases->fHelFinal[3];
  Double_t colorp  = bases->f1Ptr->GetColor();
  Double_t chgu    = cp*bases->f1Ptr->GetCharge();
  Double_t chgd    = cp*bases->f2Ptr->GetCharge();
  Int_t    islev   = colorp > 1. ? 201 : 0;
  Int_t    icfp    = 2;

  Double_t mrn     = bases->GetMass();

  Int_t    idw     = cp*24;
  Double_t mw      = bases->fW1BosonPtr->GetMass();

  new (partons[0]) JSFSpringParton(1, idnu, mnu,    0, *qp[0], 0, 0, 0, helnu,    0,     0);
  new (partons[1]) JSFSpringParton(2,  idx, mrn,    0, *qp[1], 2, 3, 0,     0,    0,     0);
  new (partons[2]) JSFSpringParton(3,  idw,  mw, chgw, *qp[2], 2, 5, 2,     0,    0,     0);
  new (partons[3]) JSFSpringParton(4,  idl,  ml, chgl, *qp[3], 0, 0, 2,  hell,    0,     0);
  new (partons[4]) JSFSpringParton(5,  idu,  mu, chgu, *qp[4], 0, 0, 3,  helu, icfp, islev);
  new (partons[5]) JSFSpringParton(6,  idd,  md, chgd, *qp[5], 0, 0, 3,  held, icfp, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class NnBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NnBases::NnBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fEcmInit   (500.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0),
           fMass      (100.),
           fMass4     (3.3e-9),
           fNkk       (1),
           fGenNu     (1),
           fGenLepton (2),
           fWpModesLo ( 1),
           fWpModesHi (12),
           fFPtr      (0),
           fZBosonPtr (0),
           fNRPtr     (0),
           fW1BosonPtr(0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2XX      (0.),
           fQ2X1      (0.),
           fW1ModePtr (0),
           f1Ptr      (0),
           f2Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fR_BW_m    (0),
           fR_BW_p    (0),
           fR_BS_m    (0),
           fR_BS_p    (0),
           fR_ISR_var (0),
           fR_ISR_side(0),
           fCosThetaW1(0.),
           fPhiW1     (0.),
           fCosThetaF1(0.),
           fPhiF1     (0.)
{
  //  Constructor of bases.  Default parameter should be initialized here
  // --------------------------------------------
  //  Get parameters from jsf.conf, if specified
  // --------------------------------------------

  cout << "Init ffbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("NnBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.CosthWRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.PhiOverPiWRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.PhiOverPiFRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.Ecm","500."));       // E_cm (0.5TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.BeamWidth","0.002")); // Beam spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.WpModesLo","1"));     // W- decay mode lo
  ins >> fWpModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.WpModesHi","12"));    // W- decay mode hi
  ins >> fWpModesHi;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.NRMass","100."));    // NR mass
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.NRMass4","3.3e-9")); // 4-dim NR mass
  ins >> fMass4;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.Nkk","1"));          //  NR KK mode number
  ins >> fNkk;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.HiggsMass","120.")); //  NR KK mode number
  ins >> fMassHiggs;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.GenNu","1"));        //  nu Generation
  ins >> fGenNu;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("NnBases.GenLepton","2"));    // Lepton Generation
  ins >> fGenLepton;


  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);

  DefineVariable(fXQ2X1  , 0., 1., 0, 1);
  DefineVariable(fXQ2W1  , 0., 1., 0, 1);
  DefineVariable(fW1DecayMode   , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();

  DefineVariable(fCosTheta   , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi        , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaW1 , fXL[2], fXU[2], 1, 1);
  DefineVariable(fPhiW1      , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaF1 , fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiF1      , fXL[5], fXU[5], 0, 1);

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
  SetIteration1(0.05, 5);
  SetIteration2(0.05, 10);

}
// --------------------------
//  D-tor
// --------------------------
NnBases::~NnBases()
{
  delete fZBosonPtr;
  delete fNRPtr;
  delete fW1BosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t NnBases::Func()
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

  fW1ModePtr = fW1BosonPtr->PickMode(fW1DecayMode, weight, fW1Mode);
#ifndef __NOWDECAY__
  bsWeight *= weight;
#endif
  f1Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(1));

  Double_t mu = f1Ptr->GetMass();
  Double_t md = f2Ptr->GetMass();

  Double_t mlep   = fFPtr->GetMass();
  Double_t mnu    = kMnu[fGenNu-1]; //takubo

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < mnu + mlep + mu + md) { 
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
  //  Calcuate differential cross section
  // --------------------------------------------
  fQ2XX = fEcmIP*fEcmIP;
  
#ifndef __ZEROWIDTH__
  Double_t rs   = fEcmIP;
  Double_t qmin = mlep + mu + md;
  Double_t qmax = rs - mnu;
  fQ2X1 = fNRPtr->GetQ2BW(qmin,qmax, fXQ2X1, weight);
#else
  fQ2X1  = TMath::Power(fNRPtr->GetMass(),2);
  weight = kPi*fNRPtr->GetMass()*fNRPtr->GetWidth();
#endif
  Double_t qx1 = TMath::Sqrt(fQ2X1);
  bsWeight *= weight;

#ifndef __ZEROWIDTH__
  qmin = mu + md;
  qmax = qx1 - mlep;
  fQ2W1 = fW1BosonPtr->GetQ2BW(qmin, qmax, fXQ2W1, weight);
#else
  fQ2W1 = TMath::Power(fW1BosonPtr->GetMass(),2);
  weight = kPi*fW1BosonPtr->GetMass()*fW1BosonPtr->GetWidth();
#ifdef __NOWDECAY__
  weight = 1.;
#endif
#endif

  Double_t qw1 = TMath::Sqrt(fQ2W1);
  bsWeight *= weight;


  GENBranch w1branch(fQ2W1, fCosThetaF1, fPhiF1, mu*mu, md*md);
  GENBranch x1branch(fQ2X1, fCosThetaW1, fPhiW1, &w1branch, mlep*mlep);
  GENBranch cmbranch(fQ2XX, fCosTheta, fPhi, &x1branch, mnu*mnu);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------
  //  cout << "CosThetaLepton= " << g_CosThetaLepton << endl;

  if (!hEcm  ) hEcm   = new TH1D( "hEcm"  , "Ecm"    , 50,     0., fEcmInit*1.1 );
  if (!hCosNR) hCosNR = new TH1D( "hCosNR", "Costh"  , 50, fXL[0],        fXU[0]);
  if (!hPhiNR) hPhiNR = new TH1D( "hPhiNR", "Phi"    , 50, fXL[1],        fXU[1]);
  if (!hMNR  ) hMNR   = new TH1D( "hMNR"  , "MNR"    , 50, 100-5.e-2,  100+5.e-2);
  if (!hCosW ) hCosW  = new TH1D( "hCosW" , "CosthW" , 50, fXL[2],        fXU[2]);
  if (!hPhiW ) hPhiW  = new TH1D( "hPhiW" , "PhiW"   , 50, fXL[3],        fXU[3]);
  if (!hCosle) hCosle = new TH1D( "hCosle","CosthLepton", 50, fXL[2],     fXU[2]);
  if (!hPhile) hPhile = new TH1D( "hPhile","PhiLepton", 50, fXL[3],      fXU[3]);
  if (!hMw   ) hMw    = new TH1D( "hMw"   , "Mw"     , 80,    20.,          100.);
  if (!hCosU ) hCosU  = new TH1D( "hCosU" , "CosthU" , 50, fXL[4],        fXU[4]);
  if (!hPhiU ) hPhiU  = new TH1D( "hPhiU" , "PhiU"   , 50, fXL[5],        fXU[5]);
  if (!hWDK  ) hWDK   = new TH1D( "hWDK"  , "W Mode" , 12,     0.,           12.);
  if (!hHelIn) hHelIn = new TH1D( "hHelIn", "Helin"  ,  2,     0.,            2.);
  if (!hHelOt) hHelOt = new TH1D( "hHelOt", "Helot"  ,  4,     0.,            4.);

  //test in order to get the angular distribution of Lepton
  
      //  std::cout << "NnSprintBuf::Saito()" << std::endl;
  int S_fNparton = 6;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [S_fNparton];
  }
  pv[0] = fP[0]; // neutrino
  pv[3] = fP[1]; // charged lepton
  pv[4] = fP[2]; // up quark
  pv[5] = fP[3]; // down quark
  pv[2] = pv[4] + pv[5]; // W
  pv[1] = pv[2] + pv[3]; // right handed neutrino

  /*
  Int_t cp = fHelCombInitial < 0.5 ? +1 : -1;
  static TVector **qp = 0;
  if (!qp) {
    qp = new TVector* [S_fNparton];
    for (Int_t i=0; i<S_fNparton; i++) qp[i] = new TVector(4);
  }
  for (Int_t i=0; i<S_fNparton; i++) {
    TVector &q = *qp[i];
    q(0) =    pv[i].E ();
    q(1) = cp*pv[i].Px();
    q(2) = cp*pv[i].Py();
    q(3) = cp*pv[i].Pz();
  }
  */
  
  g_CosThetaLepton = pv[3].CosTheta();
  g_PhiLepton      = pv[3].Phi();

  // test end

  hEcm  ->Fill(fEcmIP           , (bsWeight*sigma));
  hCosNR->Fill(fCosTheta        , (bsWeight*sigma));
  hPhiNR->Fill(fPhi             , (bsWeight*sigma));
  hMNR  ->Fill(qx1              , (bsWeight*sigma));
  hCosW ->Fill(fCosThetaW1      , (bsWeight*sigma));
  hPhiW ->Fill(fPhiW1           , (bsWeight*sigma));
  hCosle->Fill(g_CosThetaLepton , (bsWeight*sigma));
  hPhile->Fill(g_PhiLepton      , (bsWeight*sigma));
  hMw   ->Fill(qw1              , (bsWeight*sigma));
  hCosU ->Fill(fCosThetaF1      , (bsWeight*sigma));
  hPhiU ->Fill(fPhiF1           , (bsWeight*sigma));
  hWDK  ->Fill((Double_t)fW1Mode, (bsWeight*sigma));
  hHelIn->Fill((Double_t)fJCombI, (bsWeight*sigma));
  hHelOt->Fill((Double_t)fJCombF, (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t NnBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------

  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());

  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  //--------------
  // CM -> N nu
  //--------------
  Double_t cosx  = cmbranch.GetCosTheta();
  Double_t phix  = cmbranch.GetPhi();
  Double_t mrnu2 = cmbranch.GetM12();
  Double_t mnu2  = cmbranch.GetM22();

  GENPhase2 phaseCM(qcm, mrnu2, mnu2, cmframe, cosx, phix, 0);
  ANL4DVector px1 = phaseCM.GetFourMomentum(0);
  fP[0] = phaseCM.GetFourMomentum(1); 
  fM[0] = TMath::Sqrt(mnu2); 

  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  //--------------
  // N -> W+ l
  //--------------
  GENBranch &x1branch = * cmbranch.GetBranchPtr(0);

  Double_t cosw  = x1branch.GetCosTheta();
  Double_t phiw  = x1branch.GetPhi();
  Double_t mw2   = x1branch.GetM12();
  Double_t mlep2 = x1branch.GetM22();

  GENFrame cmframe2 = phaseCM.GetFrame();
  GENPhase2 phaseX1(px1, mw2, mlep2, cmframe2, cosw, phiw, 1);
  fP[1] = phaseX1.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(mlep2);
  ANL4DVector pw = phaseX1.GetFourMomentum(0);
  Double_t betaw = phaseX1.GetBetaBar();
  if (betaw <= 0.) return 0;

  //----------------
  // W+ -> fu + fdb
  //----------------
  GENBranch &wbranch =* x1branch.GetBranchPtr(0);
  Double_t cosf = wbranch.GetCosTheta();
  Double_t phif = wbranch.GetPhi();
  Double_t mu2  = wbranch.GetM12();
  Double_t md2  = wbranch.GetM22();
  GENFrame x1frame = phaseX1.GetFrame();
  GENPhase2 phaseW(pw, mu2, md2, x1frame, cosf, phif, 1);
  fP[2] = phaseW.GetFourMomentum(0);
  fM[2] = TMath::Sqrt(mu2);
  fP[3] = phaseW.GetFourMomentum(1);
  fM[3] = TMath::Sqrt(md2);
  Double_t betaf = phaseW.GetBetaBar();
  if (betaf <= 0.) return 0;

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
#if 0
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
#endif
  ANL4DVector pcm = fP[0] + fP[1] + fP[2] + fP[3];
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
#ifndef __NOWDECAY__
  static const Int_t    kNbr  = 3;
#else
  static const Int_t    kNbr  = 2;
  betaf = 1.;
#endif
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,2*kNbr+3));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax * betaw * betaf; // phase space factor
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
Double_t NnBases::AmpSquared(GENBranch &cmbranch)
{
#ifndef __NOWDECAY__
  Double_t  color = f1Ptr->GetColor();
  Int_t     ig1   = f1Ptr->GetGenNo() - 1;
  Int_t     ig2   = f2Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig1][ig2],2);
#else
  Double_t  color = 1.;
  Double_t  mix   = 1.;
#endif
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t NnBases::FullAmplitude()
{
  //-----------------------------------------
  // e+e- -> nu_jnu N_jnu -> nu_jnu W+ l_il 
  //-----------------------------------------
  Double_t glw = -kGw*kSqh; //takubo
  Double_t grw = 0.; //takubo

  static const Bool_t kIsIncoming = kTRUE;  
  static const Bool_t kIsOutgoing = kFALSE;
  
  HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); 
  HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);
  
  HELFermion fn (fP[0], fM[0], fHelFinal [0], -1, kIsIncoming); // neutrino incoming line
  HELFermion fl (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing); // charged lepton outgoing line
#ifndef __NOWDECAY__
  HELFermion fu (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // up quark outgoing line
  HELFermion fdb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // down bar quark incoming line
  Double_t gamw  = fW1BosonPtr->GetWidth(); 
  HELVector  wp (fdb, fu, glw, grw, kM_w, gamw);                // (W+ -> u dbar) internal line
#else
  ANL4DVector pwp = fP[2] + fP[3];
  HELVector   wp (pwp, kM_w, fHelFinal[2], +1);               // W+
#endif

  Complex_t amp = 0.;
  Double_t mrnu     = fNRPtr->GetMass();   // R-neutrino mass
  Double_t gamrn    = fNRPtr->GetWidth();  // R-neutrino width
  Double_t glrnwl   = fNRPtr->GetGwl()[0];
  Double_t grrnwl   = fNRPtr->GetGwl()[1];

  HELFermion fnr(fl, wp, glrnwl, grrnwl, mrnu, gamrn);  // (N -> W+ l-) internal line
  Double_t glznrn   = fNRPtr->GetGzn(fGenNu)[0];
  Double_t grznrn   = fNRPtr->GetGzn(fGenNu)[1];
  Double_t glrnwe   = fNRPtr->GetGwe()[0];
  Double_t grrnwe   = fNRPtr->GetGwe()[1];
  Double_t glnwe    = fGenNu == 1 ? glw : 0.;
  Double_t grnwe    = fGenNu == 1 ? grw : 0.;

  amp += AmpEEtoNn(em, ep, fnr, fn, 
                   glznrn, grznrn,
                   glrnwe, grrnwe, 
                   glnwe , grnwe,
                   mrnu  , gamrn);
#if 0
    cerr << " glrnwl = " << glrnwl
         << " grrnwl = " << grrnwl
	 << " glznrn = " << glznrn
	 << " grznrn = " << grznrn
	 << " glnwe  = " << glnwe
	 << " grnwe  = " << grnwe
	 << " mrnu   = " << mrnu
	 << " gamrn  = " << gamrn
	 << endl;
    cerr << "j = " << j << " amp = " << amp << endl;

#endif
  return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNn()
// --------------------------
Complex_t NnBases::AmpEEtoNn(const HELFermion &em,
                             const HELFermion &ep,
                             const HELFermion &fnr,
			     const HELFermion &fn,
			     Double_t glznrn,
			     Double_t grznrn,
			     Double_t glrnwe,
			     Double_t grrnwe,
			     Double_t glnwe,
			     Double_t grnwe,
			     Double_t mrnu,
			     Double_t gamrn)

{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamw  = fW1BosonPtr->GetWidth();
   Double_t  gamz  = fZBosonPtr->GetWidth();
  
   //--------------------------------------------------
   // Calculate Amplitudes
   //--------------------------------------------------
 
   Complex_t amp = 0.;

   if (glznrn != 0. || grznrn != 0.) {
     // S-channel
     HELVector zs(em, ep, glze, grze, kM_z, gamz);
     HELVertex ampsch(fn, fnr, zs, glznrn, grznrn);
     amp += ampsch;
   } 

   if ((glnwe != 0. || grnwe != 0.) && 
       (glrnwe != 0. || grrnwe != 0.)) {
     // T-channel
     HELVector wt(em, fnr, glrnwe, grrnwe, kM_w, gamw);
     HELVertex amptch(fn, ep, wt, glnwe, grnwe);
     amp += -amptch;
   }

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void NnBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("NnBases.BeamstrahlungFilepath",
                                           "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("NnBases.BeamstrahlungFilename","trc500"));
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
  //  Initialize W/Z decay table
  // --------------------------------------------
  if (!fW1BosonPtr) fW1BosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fW1BosonPtr->GetEntries(); m++) {
    GENDecayMode *mp = fW1BosonPtr->GetMode(m); 
    if (mp && (m<fWpModesLo || m>fWpModesHi)) {
      mp->Lock();
    }
  }
  fW1BosonPtr->DebugPrint();

  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  fZBosonPtr->DebugPrint();

  if (!fNRPtr)  fNRPtr = new RNeutrino(fMass, fMass4, fGenLepton, fNkk, fMassHiggs);
  fNRPtr->DebugPrint();

  fFPtr = static_cast<GENPDTEntry *>(fNRPtr->GetMode(fGenLepton)->At(0));

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
#if 0
  H1Init( "hEcm"  , "Ecm"    , 50,     0., fEcmInit*1.1 );
  H1Init( "hCosNR", "Costh"  , 50, fXL[0],        fXU[0]);
  H1Init( "hPhiNR", "Phi"    , 50, fXL[1],        fXU[1]);
  H1Init( "hMNR"  , "MNR"    , 50,   299.,          301.);
  H1Init( "hCosW" , "CosthW" , 50, fXL[2],        fXU[2]);
  H1Init( "hPhiW" , "PhiW"   , 50, fXL[3],        fXU[3]);
  H1Init( "hCosle" , "Costhle" , 50, fXL[2],        fXU[2]);
  H1Init( "hPhiW" , "PhiW"   , 50, fXL[3],        fXU[3]);
  H1Init( "hMw"   , "Mw"     , 50,    60.,          100.);
  H1Init( "hCosU" , "CosthU" , 50, fXL[4],        fXU[4]);
  H1Init( "hPhiU" , "PhiU"   , 50, fXL[5],        fXU[5]);
  H1Init( "hWDK"  , "W Mode" , 12,     0.,           12.);
  H1Init( "hHelIn", "Helin" ,  2,      0.,            2.);
#endif
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void NnBases::Userout()
{
  cout << "End of NnBases----------------------------------- "     << endl
       << "Ecm                   = " << fEcmInit << " [GeV]   "    << endl
       << "Electron Polarization = " << fPole                      << endl
       << "Neutrino Generation   = " << fGenNu                     << endl
       << "Lepton Generation     = " << fGenLepton                 << endl
       << "NR mass               = " << fMass    << " [GeV]   "    << endl
       << "4-dim NR mass         = " << fMass4   << " [GeV]   "    << endl
       << "Beamstrahlung         = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremsstrahlung        = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section   = " << GetEstimate()  << " +/- "
                                     << GetError()     << " [fb]"  << endl
       << "Number of iterations  = " << GetNoOfIterate()           << endl;

  hEcm  ->SetMinimum(0.); hEcm  ->Write();
  hCosNR->SetMinimum(0.); hCosNR->Write();
  hPhiNR->SetMinimum(0.); hPhiNR->Write();
  hMNR  ->SetMinimum(0.); hMNR  ->Write();
  hCosW ->SetMinimum(0.); hCosW ->Write();
  hPhiW ->SetMinimum(0.); hPhiW ->Write();
  hCosle->SetMinimum(0.); hCosle->Write();
  hPhile->SetMinimum(0.); hPhile->Write();
  hMw   ->SetMinimum(0.); hMw   ->Write();
  hCosU ->SetMinimum(0.); hCosU ->Write();
  hPhiU ->SetMinimum(0.); hPhiU ->Write();
  hWDK  ->SetMinimum(0.); hWDK  ->Write();
  hHelIn->SetMinimum(0.); hHelIn->Write();
  hHelOt->SetMinimum(0.); hHelOt->Write();
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void NnBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					   {+1, -1}};
                                             
#ifndef __NOWDECAY__
   static const Int_t kNf = 1;
   //                                      nub,  l, fu, fdb 
   static const Int_t kFHelComb[kNf][4] = {{+1, -1, -1, +1}};
#else
   static const Int_t kNf = 3;
   //                                      nub,  l, W, unused 
   static const Int_t kFHelComb[kNf][4] = {{+1, -1, -1, +1},
                                           {+1, -1,  0, +1},
                                           {+1, -1, +1, +1}};
#endif
  
   Double_t helm = (1. - fPole)/2.;

   Int_t    cp   = fHelCombInitial < 0.5 ? +1 : -1; // cp flipping flag
   Double_t xhel = cp > 0 ? 2*fHelCombInitial : 2*fHelCombInitial-1;
   weight = 2;
   if (xhel < helm) {
      fJCombI = cp > 0 ? 0 : 1;

   } else {
      fJCombI = cp > 0 ? 1 : 0;
   }
   fHelInitial[0] = kIHelComb[fJCombI][0];
   fHelInitial[1] = kIHelComb[fJCombI][1];
   fJCombF = (Int_t)(fHelCombFinal*kNf);
   fJCombF = TMath::Min(fJCombF, kNf-1);
   fHelFinal  [0] = kFHelComb[fJCombF][0];
   fHelFinal  [1] = kFHelComb[fJCombF][1];
   fHelFinal  [2] = kFHelComb[fJCombF][2];
   fHelFinal  [3] = kFHelComb[fJCombF][3];
   weight *= kNf;
}

//_____________________________________________________________________________
// ==============================
//  class RNeutrino 
// ==============================
// Notice:
//   3 generations of NR are degenerate in the model assumed.
//   In order to simplify the calculation, the NR basis is so 
//   chosen that the final-state lepton only couples to NR_2.
//   Consequently the c-tor of this class requires the generation
//   number of the final-state lepton to be observed.
//_____________________________________________________________________________
// --------------------------
//  c-tor
// --------------------------
RNeutrino::RNeutrino(Double_t mnr,    // NR mass
                     Double_t m4,     // 4-dim NR mass
                     Int_t    gen,    // final-state lepton generation
                     Int_t    kkmode, // KK mode number
                     Double_t mh)     // Higgs mass
{
  fName    = TString("RN");
  fPID     = 20000000;
  fCharge  = 0.;
  fSpin    = 0.5;
  fMass    = mnr;
  fMass4   = m4;
  fGen     = gen;
  fN       = kkmode;
  fMassH   = mh;
  fIsoSpin = 0.;
  fColor   = 1.0;

  Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void RNeutrino::Initialize()
{
  Double_t a = -(kGw*kSqh)*(1/(kPi*(2.*fN-1.)/2.))*sqrt(2./fMass4);
  Double_t b = -(kGz/2)*(1/(kPi*(2.*fN-1.)/2.))*sqrt(2./fMass4);
  Double_t c = -(kGw/2)*(fMass/kM_w)*(1/(kPi*(2.*fN-1.)/2.))*sqrt(2./fMass4);

  Double_t x12, x22, x32;
  switch (fGen) { 
  // e
    case 1:
           x12 = TMath::Sqrt(kMnu[0]*kMNS[0][0]*kMNS[0][0] +
                             kMnu[1]*kMNS[0][1]*kMNS[0][1] +
                             kMnu[2]*kMNS[0][2]*kMNS[0][2]);
           x22 =            (kMnu[0]*kMNS[0][0]*kMNS[1][0] +
                             kMnu[1]*kMNS[0][1]*kMNS[1][1] +
                             kMnu[2]*kMNS[0][2]*kMNS[1][2])/x12;
           x32 =            (kMnu[0]*kMNS[0][0]*kMNS[2][0] +
                             kMnu[1]*kMNS[0][1]*kMNS[2][1] +
                             kMnu[2]*kMNS[0][2]*kMNS[2][2])/x12;
	   break;
  // mu
    case 2:  
           x22 = TMath::Sqrt(kMnu[0]*kMNS[1][0]*kMNS[1][0] +
                             kMnu[1]*kMNS[1][1]*kMNS[1][1] +
                             kMnu[2]*kMNS[1][2]*kMNS[1][2]);
           x32 =            (kMnu[0]*kMNS[1][0]*kMNS[2][0] +
                             kMnu[1]*kMNS[1][1]*kMNS[2][1] +
                             kMnu[2]*kMNS[1][2]*kMNS[2][2])/x22;
           x12 =            (kMnu[0]*kMNS[1][0]*kMNS[0][0] +
                             kMnu[1]*kMNS[1][1]*kMNS[0][1] +
                             kMnu[2]*kMNS[1][2]*kMNS[0][2])/x22;
	   break;
  // tau
    case 3: 
           x32 = TMath::Sqrt(kMnu[0]*kMNS[2][0]*kMNS[2][0] +
                             kMnu[1]*kMNS[2][1]*kMNS[2][1] +
                             kMnu[2]*kMNS[2][2]*kMNS[2][2]);
           x12 =            (kMnu[0]*kMNS[2][0]*kMNS[0][0] +
                             kMnu[1]*kMNS[2][1]*kMNS[0][1] +
                             kMnu[2]*kMNS[2][2]*kMNS[0][2])/x32;
           x22 =            (kMnu[0]*kMNS[2][0]*kMNS[1][0] +
                             kMnu[1]*kMNS[2][1]*kMNS[1][1] +
                             kMnu[2]*kMNS[2][2]*kMNS[1][2])/x32;
	   break;
    default:
           cerr << " Invalid fGen = " << fGen << endl
                << " abort! " << endl; 
	   ::abort();
           break;
  }
  // G_NR-W-l coupling at the decay vertex
  fGwl[0] = a*(fGen == 1 ? x12 : (fGen == 2 ? x22 : x32));
  fGwl[1] = 0.;

  // G_NR-W-e coupling at the (t-channel) production vertex
  fGwe[0] = a*x12;
  fGwe[1] = 0.;

  // G_NR-Z-n coupling at the (s-channel) production vertex
  fGzn[0][0] = b*x12;
  fGzn[1][0] = 0.; 
  fGzn[0][1] = b*x22;
  fGzn[1][1] = 0.;
  fGzn[0][2] = b*x32;
  fGzn[1][2] = 0.;
  //--
  // N --> l + W 
  //--
  Double_t gwl[3];
  gwl[0] = a*x12;
  gwl[1] = a*x22;
  gwl[2] = a*x32;
  for (Int_t ig=0; ig<3; ig++) {
    const Char_t   *namd = kName[0][1][ig];
          Int_t     pidd = kPID [0][1][ig];
	  Double_t  qfd  = -1;
	  Double_t  t3d  = -0.5;
	  Double_t  spin = 0.5;	  
	  Double_t  md   = kMass[0][1][ig];
	  Double_t  cf   = 1.;

    GENDecayMode *dmp;
    GENPDTEntry  *d1p;
    GENPDTEntry  *d2p;
    d1p = new GENPDTEntry(namd, pidd, qfd, spin, md, ig+1, t3d, cf);
    d2p = new GENPDTWBoson();
  
    Double_t m1    = d1p->GetMass();
    Double_t m2    = d2p->GetMass();
    Double_t ident = 1.;
    Double_t a     = gwl[ig];
    Double_t gam   = GamToFV(m1, m2, a)/ident;

    dmp = new GENDecayMode(gam);
    dmp->Add(d1p);
    dmp->Add(d2p);
    Add(dmp);
  }
  //--
  // N --> n + Z 
  //--
  Double_t gzn[3];
  gzn[0] = b*x12;
  gzn[1] = b*x22;
  gzn[2] = b*x32;
  for (Int_t ig=0; ig<3; ig++) {
    const Char_t   *namd = kName[0][0][ig];
          Int_t     pidd = kPID [0][0][ig];
	  Double_t  qfd  = 0;
	  Double_t  t3d  = +0.5;
	  Double_t  spin = 0.5;	  
	  Double_t  md   = kMass[0][0][ig];
	  Double_t  cf   = 1.;

    GENDecayMode *dmp;
    GENPDTEntry  *d1p;
    GENPDTEntry  *d2p;
    d1p = new GENPDTEntry(namd, pidd, qfd, spin, md, ig+1, t3d, cf);
    d2p = new GENPDTZBoson();
  
    Double_t m1    = d1p->GetMass();
    Double_t m2    = d2p->GetMass();
    Double_t ident = 1.;
    Double_t a     = gzn[ig];
    Double_t gam   = GamToFV(m1, m2, a)/ident;

    dmp = new GENDecayMode(gam);
    dmp->Add(d1p);
    dmp->Add(d2p);
    Add(dmp);
  }
  //--
  // N --> n + H 
  //--
  Double_t ghn[3];
  ghn[0] = c*x12;
  ghn[1] = c*x22;
  ghn[2] = c*x32;
  for (Int_t ig=0; ig<3; ig++) {
    const Char_t   *namd = kName[0][0][ig];
          Int_t     pidd = kPID [0][0][ig];
	  Double_t  qfd  = 0;
	  Double_t  t3d  = +0.5;
	  Double_t  spin = 0.5;	  
	  Double_t  md   = kMass[0][0][ig];
	  Double_t  cf   = 1.;
	  Double_t  mh   = fMassH;

    GENDecayMode *dmp;
    GENPDTEntry  *d1p;
    GENPDTEntry  *d2p;
    d1p = new GENPDTEntry(namd, pidd, qfd, spin, md, ig+1, t3d, cf);
    d2p = new GENPDTEntry( "h",   25,  0.,   0., mh,    1,-0.5,  1); // dummy H
  
    Double_t m1    = d1p->GetMass();
    Double_t m2    = d2p->GetMass();
    Double_t ident = 1.;
    Double_t a     = ghn[ig];
    Double_t gam   = GamToFS(m1, m2, a)/ident;

    dmp = new GENDecayMode(gam);
    dmp->Add(d1p);
    dmp->Add(d2p);
    Add(dmp);
  }
}

// --------------------------
//  GamToFV: N -> l W or n Z
// --------------------------
Double_t RNeutrino::GamToFV(Double_t m1, // 1st daughter mass (l/n)
			    Double_t m2, // 2nd daughter mass (W/Z)
			    Double_t a)  // coupling
{
  if (m1+m2 >= fMass) return 0.;
  Double_t x1   = TMath::Power(m1/fMass,2);
  Double_t x2   = TMath::Power(m2/fMass,2);
  Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

  if (beta <= 0.) return 0.;
  beta = TMath::Sqrt(beta);

#if 1
  Double_t gam = a*a*fMass/(32*kPi*x2)*((1-x1)*(1-x1)+(1+x1)*x2-2*x2*x2)*beta;
#else
  ANL4DVector px(fMass,0.,0.,0.);
  double p = fMass*beta/2;
  ANL4DVector pw(sqrt(p*p+m2*m2),0.,0.,p);
  ANL4DVector pl(sqrt(p*p+m1*m1),0.,0.,-p);

  static const Bool_t kIsIncoming = kTRUE;
  static const Bool_t kIsOutgoing = kFALSE;
  HELFermion fxl(px, fMass, -1, +1, kIsIncoming);
  HELFermion fxr(px, fMass, +1, +1, kIsIncoming);
  HELFermion fl (pl, m1   , -1, +1, kIsOutgoing);
  HELVector wm(pw, m2, -1, +1);
  HELVector w0(pw, m2,  0, +1);
  HELVector wp(pw, m2, +1, +1);

  Double_t ar = 0.;
  HELVertex amplm(fxl, fl, wm, a, ar);
  HELVertex ampl0(fxl, fl, w0, a, ar);
  HELVertex amplp(fxl, fl, wp, a, ar);
  HELVertex amprm(fxr, fl, wm, a, ar);
  HELVertex ampr0(fxr, fl, w0, a, ar);
  HELVertex amprp(fxr, fl, wp, a, ar);

  Double_t tta = TMath::Power(abs(amplm),2)
	       + TMath::Power(abs(ampl0),2)
	       + TMath::Power(abs(amplp),2)
               + TMath::Power(abs(amprm),2)
	       + TMath::Power(abs(ampr0),2)
	       + TMath::Power(abs(amprp),2);
  Double_t spin = 1./2;
  Double_t gam  = (1/(2*fMass))*spin*tta*(beta/(8*kPi));
#endif

  return gam;
}

// --------------------------
//  GamToFS: N -> n h
// --------------------------
Double_t RNeutrino::GamToFS(Double_t m1, // 1st daughter mass (n)
			    Double_t m2, // 2nd daughter mass (h)
			    Double_t a)  // coupling
{
  if (m1+m2 >= fMass) return 0.;
  Double_t x1   = TMath::Power(m1/fMass,2);
  Double_t x2   = TMath::Power(m2/fMass,2);
  Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

  if (beta <= 0.) return 0.;
  beta = TMath::Sqrt(beta);

  Double_t gam = a*a* fMass*beta*beta /(32*kPi);

  return gam;
}
