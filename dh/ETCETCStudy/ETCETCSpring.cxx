//*****************************************************************************
//* =====================
//*  ETCETCSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- --> XX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/17  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ETCETCSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PAHSESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ETCETCSpring)
ClassImp(ETCETCSpringBuf)
ClassImp(ETCETCBases)
ClassImp(ETCBoson)
ClassImp(ETIBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class ETCETCSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETCETCSpring::ETCETCSpring(const char      *name,
                           const char      *title,
                               ETCETCBases *bases)
            : JSFSpring(name, title, bases)
{
  fEventBuf = new ETCETCSpringBuf("ETCETCSpringBuf",
                                  "ETCETCSpring event buffer",
                                  this);
  if (!bases) { 
    SetBases(new ETCETCBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ETCETCSpring::~ETCETCSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class ETCETCSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ETCETCSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ETCETCBases  *bases   = (ETCETCBases*)((ETCETCSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 10;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];    // eta_I from eta-
  pv[6] = bases->fP[1];    // fub from W-
  pv[7] = bases->fP[2];    // fd  from W-
  pv[3] = pv[6] + pv[7];   // W-
  pv[4] = bases->fP[3];    // eta_I from eta+
  pv[8] = bases->fP[4];    // fu  from W+
  pv[9] = bases->fP[5];    // fdb from W+
  pv[5] = pv[8] + pv[9];   // W+
  pv[0] = pv[2] + pv[3];   // eta-
  pv[1] = pv[4] + pv[5];   // eta+

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2XX         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X1         = bases->GetQ2X1();
  fCosThetaW1   = bases->GetCosThetaW1();
  fPhiW1        = bases->GetPhiW1();
  fQ2W1         = bases->GetQ2W1();
  fCosThetaF1   = bases->GetCosThetaF1();
  fPhiF1        = bases->GetPhiF1();
  fQ2X2         = bases->GetQ2X2();
  fCosThetaW2   = bases->GetCosThetaW2();
  fPhiW2        = bases->GetPhiW2();
  fQ2W2         = bases->GetQ2W2();
  fCosThetaF2   = bases->GetCosThetaF2();
  fPhiF2        = bases->GetPhiF2();
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
  Int_t    idx     = bases->fXBosonPtr ->GetPID(); // PDG code for ETC? 
  Int_t    iddm    = bases->fDMBosonPtr->GetPID(); // PDG code for ETI? 

  Int_t    idw     = 24;                           // PDG code for W
  Int_t    id1     = bases->f1Ptr->GetPID();       // PDG code for 1st daughter
  Int_t    id2     = bases->f2Ptr->GetPID();       // PDG code for 2nd daughter
  Double_t m1      = bases->f1Ptr->GetMass();      // 1st daughter mass
  Double_t m2      = bases->f2Ptr->GetMass();      // 2nd daughter mass
  Double_t colorm  = bases->f1Ptr->GetColor();     // color factor
  Double_t chg1    = bases->f1Ptr->GetCharge();    // 1st daughter charge
  Double_t chg2    = bases->f2Ptr->GetCharge();    // 2nd daughter charge
  Int_t    islevm  = colorm > 1. ? 101 : 0;        // shower level
  Int_t    icfm    = 1;                            // color flux id

  Int_t    id3     = bases->f3Ptr->GetPID();       // PDG code for 1st daughter
  Int_t    id4     = bases->f4Ptr->GetPID();       // PDG code for 2nd daughter
  Double_t m3      = bases->f3Ptr->GetMass();      // 1st daughter mass
  Double_t m4      = bases->f4Ptr->GetMass();      // 2nd daughter mass
  Double_t colorp  = bases->f3Ptr->GetColor();     // color factor
  Double_t chg3    = bases->f3Ptr->GetCharge();    // 1st daughter charge
  Double_t chg4    = bases->f4Ptr->GetCharge();    // 2nd daughter charge
  Int_t    islevp  = colorp > 1. ? 201 : 0;        // shower level
  Int_t    icfp    = 2;                            // color flux id

  //Double_t rq2etm  = pv[0].Mag();
  //Double_t rq2etp  = pv[1].Mag();
  Double_t rq2wm   = pv[3].Mag();
  Double_t rq2wp   = pv[5].Mag();

  Double_t mass    = bases->GetMass();
  Double_t msdm    = bases->GetMassDM();

  //                               No. PID   Mass   Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton( 1, idx , mass ,   -1., *qp[0], 2, 3,  0, 0,   0,     0);
  new (partons[1]) JSFSpringParton( 2,-idx , mass ,   +1., *qp[1], 2, 5,  0, 0,   0,     0);
  new (partons[2]) JSFSpringParton( 3, iddm, msdm ,    0., *qp[2], 0, 0,  1, 0,   0,     0);
  new (partons[3]) JSFSpringParton( 4, idw , rq2wm,   -1., *qp[3], 2, 7,  1, 0,   0,     0);
  new (partons[4]) JSFSpringParton( 5, iddm, msdm ,    0., *qp[4], 0, 0,  1, 0,   0,     0);
  new (partons[5]) JSFSpringParton( 6,-idw , rq2wp,   +1., *qp[5], 2, 9,  2, 0,   0,     0);
  new (partons[6]) JSFSpringParton( 7, id1 ,    m1,  chg1, *qp[6], 0, 0,  4, 0,icfm,islevm);
  new (partons[7]) JSFSpringParton( 8, id2 ,    m2,  chg2, *qp[7], 0, 0,  4, 0,icfm,islevm);
  new (partons[8]) JSFSpringParton( 9, id3 ,    m3,  chg3, *qp[8], 0, 0,  6, 0,icfp,islevp);
  new (partons[9]) JSFSpringParton(10, id4 ,    m4,  chg4, *qp[9], 0, 0,  6, 0,icfp,islevp);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ETCETCBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETCETCBases::ETCETCBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 180.),
           fMassDM    (  60.),
           fEcmInit   ( 500.),
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
           fW1ModePtr (0),
           f1Ptr      (0),
           f2Ptr      (0),
           fQ2W1      (0.),
           fW2ModePtr (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fQ2W2      (0.),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaW1(0.),
           fPhiW1     (0.),
           fXQ2X1     (0.),
           fCosThetaF1(0.),
           fPhiF1     (0.),
           fXQ2W1     (0.),
           fCosThetaW2(0.),
           fPhiW2     (0.),
           fXQ2X2     (0.),
           fCosThetaF2(0.),
           fPhiF2     (0.),
           fXQ2W2     (0.),
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

  cout << "Init etcetcbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("ETCETCBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.CosthWmRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.PhiWmOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.CosthFmRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.PhiFmOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.CosthWpRange","-1.0 1.0"));
  ins >> fXL[6] >> fXU[6];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.PhiWpOverPiRange","0.0 2.0"));
  ins >> fXL[7] >> fXU[7];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.CosthFpRange","-1.0 1.0"));
  ins >> fXL[8] >> fXU[8];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.PhiFpOverPiRange","0.0 2.0"));
  ins >> fXL[9] >> fXU[9];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.MassX","180.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.MassDM","60.")); 	 // M_dm [GeV]
  ins >> fMassDM;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETCETCBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("ETCETCBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ETCETCBases.BeamstrahlungFilename","trc500"));
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
  DefineVariable(fW1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fW2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2X1         , 0., 1., 0, 1);
  DefineVariable(fXQ2X2         , 0., 1., 0, 1);
  DefineVariable(fXQ2W1         , 0., 1., 0, 1);
  DefineVariable(fXQ2W2         , 0., 1., 0, 1);
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

  DefineVariable(fCosTheta  , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi       , fXL[1], fXU[1], 0, 0);
  DefineVariable(fCosThetaW1, fXL[2], fXU[2], 0, 0);
  DefineVariable(fPhiW1     , fXL[3], fXU[3], 0, 0);
  DefineVariable(fCosThetaF1, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiF1     , fXL[5], fXU[5], 0, 1);
  DefineVariable(fCosThetaW2, fXL[6], fXU[6], 0, 0);
  DefineVariable(fPhiW2     , fXL[7], fXU[7], 0, 0);
  DefineVariable(fCosThetaF2, fXL[8], fXU[8], 0, 1);
  DefineVariable(fPhiF2     , fXL[9], fXU[9], 0, 1);

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
  SetIteration2(0.05, 50);

}
// --------------------------
//  D-tor
// --------------------------
ETCETCBases::~ETCETCBases()
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
Double_t ETCETCBases::Func()
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

  fW1ModePtr = fWBosonPtr->PickMode(fW1DecayMode, weight, fW1Mode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fW1ModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  fW2ModePtr = fWBosonPtr->PickMode(fW2DecayMode, weight, fW2Mode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fW2ModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fW2ModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < 2*fMassDM + m1 + m2 + m3 + m4) {
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
  Double_t qmin = fMassDM + m1 + m2;
  Double_t qmax = rs - (fMassDM + m3 + m4);
#ifndef __ZEROWIDTH__
  fQ2X1 = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2X1, weight);
#else
  fQ2X1 = TMath::Power(fXBosonPtr->GetMass(),2);
  weight = kPi*fXBosonPtr->GetMass()*fXBosonPtr->GetWidth();
#endif
  Double_t qx1 = TMath::Sqrt(fQ2X1);
  bsWeight *= weight;

#ifndef __ZEROWIDTH__
  qmin  = fMassDM + m3 + m4;
  qmax  = rs - qx1;
  fQ2X2 = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2X2, weight);
#else
  fQ2X2 = TMath::Power(fXBosonPtr->GetMass(),2);
  weight = kPi*fXBosonPtr->GetMass()*fXBosonPtr->GetWidth();
#endif
  Double_t qx2 = TMath::Sqrt(fQ2X2);
  bsWeight *= weight;

  qmin  = m1 + m2;
  qmax  = qx1 - fMassDM;
  fQ2W1 = fWBosonPtr->GetQ2BW(qmin, qmax, fXQ2W1, weight);
  bsWeight *= weight;

  qmin  = m3 + m4;
  qmax  = qx2 - fMassDM;
  fQ2W2 = fWBosonPtr->GetQ2BW(qmin, qmax, fXQ2W2, weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch w1branch(fQ2W1, fCosThetaF1, fPhiF1, m1*m1    , m2*m2          );
  GENBranch x1branch(fQ2X1, fCosThetaW1, fPhiW1, &w1branch, fMassDM*fMassDM);
  GENBranch w2branch(fQ2W2, fCosThetaF2, fPhiF2, m3*m3    , m4*m4          );
  GENBranch x2branch(fQ2X2, fCosThetaW2, fPhiW2, &w2branch, fMassDM*fMassDM);
  GENBranch cmbranch(fQ2XX, fCosTheta  , fPhi  , &x1branch, &x2branch      );

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X1), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaW1       , (bsWeight*sigma));
  Xh_fill( 6, fPhiW1            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2W1), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaF1       , (bsWeight*sigma));
  Xh_fill( 9, fPhiF1            , (bsWeight*sigma));
  Xh_fill(10, TMath::Sqrt(fQ2X2), (bsWeight*sigma));
  Xh_fill(11, fCosThetaW2       , (bsWeight*sigma));
  Xh_fill(12, fPhiW2            , (bsWeight*sigma));
  Xh_fill(13, TMath::Sqrt(fQ2W2), (bsWeight*sigma));
  Xh_fill(14, fCosThetaF2       , (bsWeight*sigma));
  Xh_fill(15, fPhiF2            , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fW1Mode , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fW2Mode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ETCETCBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  // CM -> X- X+
  Double_t q2x1  = cmbranch.GetM12();
  Double_t q2x2  = cmbranch.GetM22();
  Double_t cosx1 = cmbranch.GetCosTheta();
  Double_t phix1 = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x1, q2x2, cmframe, cosx1, phix1, 0);
  ANL4DVector px1 = phaseCM.GetFourMomentum(0);
  ANL4DVector px2 = phaseCM.GetFourMomentum(1);
  Double_t betax  = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  // X- -> W- X0
  GENBranch &x1branch = * cmbranch.GetBranchPtr(0);
  Double_t cosw1 = x1branch.GetCosTheta();
  Double_t phiw1 = x1branch.GetPhi     ();
  Double_t mw12  = x1branch.GetM12();
  Double_t md12  = x1branch.GetM22();
  GENPhase2 phaseX1(px1, mw12, md12, cmframe, cosw1, phiw1, 1);
  fP[0]           = phaseX1.GetFourMomentum(1);
  fM[0]           = TMath::Sqrt(md12);
  ANL4DVector pw1 = phaseX1.GetFourMomentum(0);
  Double_t betaw1 = phaseX1.GetBetaBar();
  if (betaw1 <= 0.) return 0.;

  // W- -> fub + fd
  GENBranch &w1branch = * x1branch.GetBranchPtr(0);
  Double_t cosf1   = w1branch.GetCosTheta();
  Double_t phif1   = w1branch.GetPhi     ();
  Double_t m12     = w1branch.GetM12();
  Double_t m22     = w1branch.GetM22();
  GENFrame x1frame = phaseX1.GetFrame();
  GENPhase2 phaseW1(pw1, m12, m22, x1frame, cosf1, phif1, 1);
  fP[1] = phaseW1.GetFourMomentum(0);
  fP[2] = phaseW1.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m12);
  fM[2] = TMath::Sqrt(m22);
  Double_t betaf1 = phaseW1.GetBetaBar();
  if (betaf1 <= 0.) return 0.;

  // X+ -> W+ X0
  GENBranch &x2branch = * cmbranch.GetBranchPtr(1);
  Double_t cosw2 = x2branch.GetCosTheta();
  Double_t phiw2 = x2branch.GetPhi     ();
  Double_t mw22  = x2branch.GetM12();
  Double_t md22  = x2branch.GetM22();
  GENPhase2 phaseX2(px2, mw22, md22, cmframe, cosw2, phiw2, 1);
  fP[3]           = phaseX2.GetFourMomentum(1);
  fM[3]           = TMath::Sqrt(md22);
  ANL4DVector pw2 = phaseX2.GetFourMomentum(0);
  Double_t betaw2 = phaseX2.GetBetaBar();
  if (betaw2 <= 0.) return 0.;

  // W- -> fub + fd
  GENBranch &w2branch = * x2branch.GetBranchPtr(0);
  Double_t cosf2   = w2branch.GetCosTheta();
  Double_t phif2   = w2branch.GetPhi     ();
  Double_t m32     = w2branch.GetM12();
  Double_t m42     = w2branch.GetM22();
  GENFrame x2frame = phaseX2.GetFrame();
  GENPhase2 phaseW2(pw2, m32, m42, x2frame, cosf2, phif2, 1);
  fP[4] = phaseW2.GetFourMomentum(0);
  fP[5] = phaseW2.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(m32);
  fM[5] = TMath::Sqrt(m42);
  Double_t betaf2 = phaseW2.GetBetaBar();
  if (betaf2 <= 0.) return 0.;

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
  ANL4DVector qw1 = fP[1] + fP[2];
  ANL4DVector qx1 = fP[0] + qw1;
  cerr << " qw1 = (" 
       << qw1.E () << ","
       << qw1.Px() << ","
       << qw1.Py() << ","
       << qw1.Pz() << ")" << endl;
  ANL4DVector qw2 = fP[4] + fP[5];
  ANL4DVector qx2 = fP[3] + qw2;
  cerr << " qw2 = (" 
       << qw2.E () << ","
       << qw2.Px() << ","
       << qw2.Py() << ","
       << qw2.Pz() << ")" << endl;
  cerr << " ---- " << endl;
  cerr << " rs = " << rs << endl;
  cerr << " md1 = " << fP[0].GetMass() << endl;
  cerr << " md2 = " << fP[3].GetMass() << endl;
  cerr << " mw1 = " << qw1.GetMass() << endl;
  cerr << " mw2 = " << qw2.GetMass() << endl;
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
  static const Int_t    kNbr  = 5;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                              // identical particle factor
  Double_t dPhase = kFact * betax * betaw1 * betaw2  // phase space factor
	                          * betaf1 * betaf2;
  Double_t flux   = 1./(2.* s * beta_e);             // beam flux factor
  Double_t spin   = 1./2.;                           // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ETCETCBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ETCETCBases::FullAmplitude()
{
   Double_t gamw   = fWBosonPtr->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0.;
   Double_t gxdmw  = kGw/2.; // imaginary coupling (i ommitted)

   Double_t mx     = fXBosonPtr->GetMass();
   Double_t gamx   = fXBosonPtr->GetWidth();

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  dm1(fP[0], +1);
   HELFermion fd (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
   HELFermion fub(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);
   HELVector  wm(fub, fd, glw, grw, kM_w, gamw);
   HELScalar  xm(wm, dm1, gxdmw, mx, gamx);

   HELScalar  dm2(fP[3], +1);
   HELFermion fu (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing);
   HELFermion fdb(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming);
   HELVector  wp(fdb, fu, glw, grw, kM_w, gamw);
   HELScalar  xp(wp, dm2, gxdmw, mx, gamx);
#ifdef __DEBUG__
   cerr << " dm1 = " << static_cast<Complex_t>(dm1) << " ";
   dm1.GetFourMomentum().DebugPrint();
   cerr << " wm : "; 
   wm.GetFourMomentum().DebugPrint();

   cerr << " dm2 = " << static_cast<Complex_t>(dm2) << " ";
   dm2.GetFourMomentum().DebugPrint();
   cerr << " wp : "; 
   wp.GetFourMomentum().DebugPrint();

   cerr << " xm = " << static_cast<Complex_t>(xm) << " ";
   xm.GetFourMomentum().DebugPrint();
   cerr << " xp = " << static_cast<Complex_t>(xp) << " ";
   xp.GetFourMomentum().DebugPrint();
#endif

   Complex_t amp = AmpEEtoXX(em, ep, xm, xp);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoXX()
// --------------------------
Complex_t ETCETCBases::AmpEEtoXX(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELScalar  &xm,
                                 const HELScalar  &xp)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  gamz  = fZBosonPtr->GetWidth();

   Double_t  gzxx  = glze;
   Double_t  gaxx  = ge;

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
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel   Z: backward-going
   HELVertex ampzxx(zs, xm, xp, gzxx);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma
   HELVertex ampaxx(as, xm, xp, gaxx);

   //--
   // Sum of the two
   //--
   Complex_t amp = ampzxx + ampaxx;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ETCETCBases::Userin()
{
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

   cerr << " Eta +/- Boson " << endl;
   if (!fXBosonPtr)  fXBosonPtr  = new ETCBoson(fMass,fMassDM);
   fXBosonPtr->DebugPrint();

   cerr << " Eta_I Boson " << endl;
   if (!fDMBosonPtr) fDMBosonPtr = new ETIBoson(fMassDM);
   fDMBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4,   100.,   300.,       50, "Mx-"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthW-");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiW-"  );
  Xh_init( 7,    60.,   100.,       50, "MW-"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthF1");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiF1"  );
  Xh_init(10,   100.,   300.,       50, "Mx+"    );
  Xh_init(11, fXL[2], fXU[2],       50, "CosthW+");
  Xh_init(12, fXL[3], fXU[3],       50, "PhiW+"  );
  Xh_init(13,    60.,   100.,       50, "MW+"    );
  Xh_init(14, fXL[4], fXU[4],       50, "CosthF2");
  Xh_init(15, fXL[5], fXU[5],       50, "PhiF2"  );
  Xh_init(16,     0.,     2.,        2, "Helin " );
  Xh_init(17,     0.,     2.,        2, "Helot " );
  Xh_init(18,     0.,    12.,       12, "W- mode");
  Xh_init(19,     0.,    12.,       12, "W+ mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ETCETCBases::Userout()
{
  cout << "End of ETCETCBases----------------------------------- "  << endl
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
void ETCETCBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][6] = {{0, +1, -1, 0, -1, +1}};
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

//-----------------------------------------------------------------------------
// ==============================
//  class ETCBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETCBoson::ETCBoson(Double_t m,
                   Double_t mdm)
        : fMassDM(mdm)
{
   fName    = TString("ETC");
   fPID     = 200000000;
   fCharge  = -1.0;
   fSpin    =  0.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin = -0.5;
   fColor   =  1.0;

   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void ETCBoson::Initialize()
{
   //--
   // Couplings
   //--
   Double_t  a  = kGw/2;
   //--
   // ETC --> ETI + W
   //--
   GENDecayMode *dmp;
   GENPDTEntry  *d1p, *d2p;
   d1p = new ETIBoson(fMassDM);
   d2p = new GENPDTWBoson();
   Double_t cf    = 1.;
   Double_t m1    = d1p->GetMass();
   Double_t m2    = d2p->GetMass();
   Double_t ident = 1.;
   Double_t gam = GamToSV(m1, m2, a, cf)/ident;
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
Double_t ETCBoson::GamToSV(Double_t m1, // 1st daughter mass
                           Double_t m2, // 2nd daughter mass
                           Double_t a,  // coupling
                           Double_t cf) // color factor
{
   Double_t x1   = TMath::Power(m1/fMass,2);
   Double_t x2   = TMath::Power(m2/fMass,2);
   Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

   if (beta <= 0.) return 0.;
   beta = TMath::Sqrt(beta);

#if 1
   Double_t p1p2 = (fMass*fMass - m1*m1 - m2*m2)/2.;
   Double_t tta = a*a*(-(4*m1*m1+m2*m2+4*p1p2)+TMath::Power(m2*m2+2*p1p2,2)/m2/m2);
#else
   ANL4DVector px(fMass,0.,0.,0.);
   double p = fMass*beta/2;
   ANL4DVector pw(sqrt(p*p+m2*m2),0.,0.,p);
   ANL4DVector pdm(sqrt(p*p+m1*m1),0.,0.,-p);
   HELScalar x(px, -1);
   HELScalar dm(pdm, +1);
   HELVector wm(pw, m2, -1, +1);
   HELVector w0(pw, m2,  0, +1);
   HELVector wp(pw, m2, +1, +1);
   HELVertex ampm(wm, x, dm, a);
   HELVertex amp0(w0, x, dm, a);
   HELVertex ampp(wp, x, dm, a);
   double tta = pow(abs(ampm),2) + pow(abs(amp0),2)+ pow(abs(ampp),2);
#endif
   Double_t fac = 1./(16.*kPi)/fMass;
   Double_t gam = fac*tta*beta*cf;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ETIBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETIBoson::ETIBoson(Double_t m)
{
   fName    = TString("ETI");
   fPID     = 200000001;
   fCharge  =  0.0;
   fSpin    =  0.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin = +0.5;
   fColor   =  1.0;
}