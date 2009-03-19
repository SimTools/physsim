//*****************************************************************************
//* =====================
//*  ETRETISpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- --> eta_R eta_I generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/17  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ETRETISpring.h"

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

ClassImp(ETRETISpring)
ClassImp(ETRETISpringBuf)
ClassImp(ETRETIBases)
ClassImp(ETRBoson)
ClassImp(ETIBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class ETRETISpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETRETISpring::ETRETISpring(const char      *name,
                           const char      *title,
                               ETRETIBases *bases)
            : JSFSpring(name, title, bases)
{
  fEventBuf = new ETRETISpringBuf("ETRETISpringBuf",
                                  "ETRETISpring event buffer",
                                  this);
  if (!bases) { 
    SetBases(new ETRETIBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ETRETISpring::~ETRETISpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ETRETISpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ETRETIBases *bs = static_cast<ETRETIBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ETRETIBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ETRETISpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ETRETISpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ETRETIBases  *bases   = (ETRETIBases*)((ETRETISpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 6;
  ANL4DVector pv[kNparton];
  pv[1] = bases->fP[0];    // eta_I
  pv[2] = bases->fP[3];    // eta_I from eta_R 
  pv[4] = bases->fP[1];    // f  from Z
  pv[5] = bases->fP[2];    // fb from Z
  pv[3] = pv[4] + pv[5];   // Z from eta_R
  pv[0] = pv[2] + pv[3];   // eta_R

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2XD         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X          = bases->GetQ2X();
  fCosThetaZ    = bases->GetCosThetaZ();
  fPhiZ         = bases->GetPhiZ();
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
  Int_t    idx     = bases->fXBosonPtr ->GetPID(); // PDG code for ETR? 
  Double_t mass    = bases->GetMass();

  Int_t    iddm    = bases->fDMBosonPtr->GetPID(); // PDG code for ETI? 
  Double_t msdm    = bases->GetMassDM();

  Int_t    idz     = 23;                           // PDG code for Z
  Int_t    idd     = bases->f1Ptr->GetPID();       // PDG code for Z daughter
  Double_t md      = bases->f1Ptr->GetMass();      // Z daughter mass
  Int_t    hel1    = bases->fHelFinal[1];          // 1st daughter helicity
  Int_t    hel2    = bases->fHelFinal[2];          // 2nd daughter helicity
  Double_t color   = bases->f1Ptr->GetColor();     // color factor
  Double_t chg     = bases->f1Ptr->GetCharge();    // Z daughter charge
  Int_t    islev   = color > 1. ? 101 : 0;         // shower level
  Int_t    icf     = 1;                            // color flux id
  Double_t rq2z    = pv[3].Mag();

  //                               No. PID   Mass Charge   pv   Nd 1st Mom   hel col shower
  new (partons[0]) JSFSpringParton( 1, idx , mass,    0., *qp[0], 2, 3,  0,   0,   0,     0);
  new (partons[1]) JSFSpringParton( 2, iddm, msdm,    0., *qp[1], 0, 0,  0,   0,   0,     0);
  new (partons[2]) JSFSpringParton( 3, iddm, msdm,    0., *qp[2], 0, 0,  1,   0,   0,     0);
  new (partons[3]) JSFSpringParton( 4, idz , rq2z,    0., *qp[3], 2, 5,  1,   0,   0,     0);
  new (partons[4]) JSFSpringParton( 5, idd , md  ,   chg, *qp[4], 0, 0,  4,hel1, icf, islev);
  new (partons[5]) JSFSpringParton( 6,-idd , md  ,  -chg, *qp[5], 0, 0,  4,hel2, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ETRETIBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETRETIBases::ETRETIBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 180.),
           fMassDM    (  60.),
           fEcmInit   ( 500.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fXBosonPtr ( 0),
           fDMBosonPtr( 0),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fPhotonPtr ( 0),
           fGluonPtr  ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2XD      (0.),
           fQ2X      (0.),
           fZModePtr (0),
           f1Ptr      (0),
           f2Ptr      (0),
           fQ2Z      (0.),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaZ(0.),
           fPhiZ     (0.),
           fXQ2X     (0.),
           fCosThetaF(0.),
           fPhiF     (0.),
           fXQ2Z     (0.),
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
  stringstream ins(gJSF->Env()->GetValue("ETRETIBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.CosthZRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.PhiZOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.MassX","180.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.MassDM","60.")); 	 // M_dm [GeV]
  ins >> fMassDM;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.BeamWidth","0.002")); // Beam width
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ETRETIBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fZDecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2X         , 0., 1., 0, 1);
  DefineVariable(fXQ2Z         , 0., 1., 0, 1);
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
  DefineVariable(fCosThetaZ , fXL[2], fXU[2], 0, 0);
  DefineVariable(fPhiZ      , fXL[3], fXU[3], 0, 0);
  DefineVariable(fCosThetaF , fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiF      , fXL[5], fXU[5], 0, 1);

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
ETRETIBases::~ETRETIBases()
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
Double_t ETRETIBases::Func()
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

  fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
#ifndef __NODECAY__
  bsWeight *= weight;
#endif
  f1Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < 2*fMassDM + m1 + m2) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average for e- is taken
  //  care of here
  SelectHelicities(weight);
#ifdef __NODECAY__
  weight = 1.;
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  fQ2XD = fEcmIP*fEcmIP;

  Double_t rs   = fEcmIP;
  Double_t qmin = fMassDM + m1 + m2;
  Double_t qmax = rs - fMassDM;
#ifndef __ZEROWIDTH__
  fQ2X = fXBosonPtr->GetQ2BW(qmin, qmax, fXQ2X, weight);
#else
  fQ2X = TMath::Power(fXBosonPtr->GetMass(),2);
  weight = kPi*fXBosonPtr->GetMass()*fXBosonPtr->GetWidth();
#ifdef __NODECAY__
  weight = 1.;
#endif
#endif
  Double_t qx = TMath::Sqrt(fQ2X);
  bsWeight *= weight;

#ifndef __ZEROWIDTH__
  qmin  = m1 + m2;
  qmax  = qx - fMassDM;
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
#else
  fQ2Z  = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#ifdef __NODECAY__
  weight = 1.;
#endif
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch zbranch (fQ2Z , fCosThetaF, fPhiF, m1*m1   , m2*m2          );
  GENBranch xbranch (fQ2X , fCosThetaZ, fPhiZ, &zbranch, fMassDM*fMassDM);
  GENBranch cmbranch(fQ2XD, fCosTheta , fPhi , &xbranch, fMassDM*fMassDM);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X) , (bsWeight*sigma));
  Xh_fill( 5, fCosThetaZ        , (bsWeight*sigma));
  Xh_fill( 6, fPhiZ             , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill( 8, fCosThetaF        , (bsWeight*sigma));
  Xh_fill( 9, fPhiF             , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(12, (Double_t)fZMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ETRETIBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  // CM -> XR XI
  Double_t q2x   = cmbranch.GetM12();
  Double_t q2dm  = cmbranch.GetM22();
  Double_t cosx  = cmbranch.GetCosTheta();
  Double_t phix  = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x, q2dm, cmframe, cosx, phix, 0);
  ANL4DVector px  = phaseCM.GetFourMomentum(0);
  fP[0]           = phaseCM.GetFourMomentum(1);
  fM[0]           = TMath::Sqrt(q2dm);
  Double_t betax  = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  // XR -> Z XI
  GENBranch &xbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosz  = xbranch.GetCosTheta();
  Double_t phiz  = xbranch.GetPhi     ();
  Double_t mz2   = xbranch.GetM12();
  Double_t mdm2  = xbranch.GetM22();
  GENPhase2 phaseX(px, mz2, mdm2, cmframe, cosz, phiz, 1);
  fP[3]          = phaseX.GetFourMomentum(1);
  fM[3]          = TMath::Sqrt(mdm2);
  ANL4DVector pz = phaseX.GetFourMomentum(0);
  Double_t betaz = phaseX.GetBetaBar();
  if (betaz <= 0.) return 0.;

  // Z -> f + fb
  GENBranch &zbranch = * xbranch.GetBranchPtr(0);
  Double_t cosf   = zbranch.GetCosTheta();
  Double_t phif   = zbranch.GetPhi     ();
  Double_t m12    = zbranch.GetM12();
  Double_t m22    = zbranch.GetM22();
  GENFrame xframe = phaseX.GetFrame();
  GENPhase2 phaseZ(pz, m12, m22, xframe, cosf, phif, 1);
  fP[1] = phaseZ.GetFourMomentum(0);
  fP[2] = phaseZ.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m12);
  fM[2] = TMath::Sqrt(m22);
  Double_t betaf = phaseZ.GetBetaBar();
  if (betaf <= 0.) return 0.;

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
  ANL4DVector qz = fP[1] + fP[2];
  ANL4DVector qx = fP[3] + qz;
  cerr << " qz = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  cerr << " ---- " << endl;
  cerr << " rs = " << rs << endl;
  cerr << " md1 = " << fP[0].GetMass() << endl;
  cerr << " md2 = " << fP[3].GetMass() << endl;
  cerr << " mz  = " << qz.GetMass() << endl;
  cerr << " mx  = " << qx.GetMass() << endl;

  ANL4DVector pcm = qx + fP[0];
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
  betaz = 1.;
  betaf = 1.;
#endif

  Double_t identp = 1.;                              // identical particle factor
  Double_t dPhase = kFact * betax * betaz * betaf;   // phase space factor
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
Double_t ETRETIBases::AmpSquared(GENBranch &cmbranch)
{
#ifndef __NODECAY__
  Double_t  color = f1Ptr->GetColor();
#else
  Double_t  color = 1.;
#endif
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ETRETIBases::FullAmplitude()
{
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf     = f1Ptr->GetCharge();
   Double_t t3f    = f1Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   Double_t gzxd   = kGz/2.; // imaginary coupling (i ommitted)

   Double_t mx     = fXBosonPtr->GetMass();
   Double_t gamx   = fXBosonPtr->GetWidth();

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

#ifndef __NODECAY__
   HELScalar  dm(fP[0], +1); // eta_I 

   HELScalar  dmx(fP[3], +1); // eta_I from eta_R
   HELFermion f (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing);
   HELFermion fb(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming);
   HELVector  z(fb, f, glz, grz, kM_z, gamz);
   HELScalar  x(z, dmx, gzxd, mx, gamx);             // eta_R
#ifdef __DEBUG__
   cerr << " z : "; 
   z.GetFourMomentum().DebugPrint();
   cerr << " dmx = " << static_cast<Complex_t>(dmx) << " ";
   dmx.GetFourMomentum().DebugPrint();
#endif
#else
   ANL4DVector px  = fP[1] + fP[2] + fP[3]; 
   HELScalar  x (px   , +1); // eta_R
   HELScalar  dm(fP[0], +1); // eta_I 
#endif
#ifdef __DEBUG__
   cerr << " dm = " << static_cast<Complex_t>(dm) << " ";
   dm.GetFourMomentum().DebugPrint();
   cerr << " x = " << static_cast<Complex_t>(x) << " ";
   x.GetFourMomentum().DebugPrint();
#endif

   Complex_t amp = AmpEEtoXD(em, ep, x, dm);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoXD()
// --------------------------
Complex_t ETRETIBases::AmpEEtoXD(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELScalar  &x,
                                 const HELScalar  &dm)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);
   Double_t  gamz  = fZBosonPtr->GetWidth();
   Double_t  gzxd  = kGz/2.;

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
   // XD Production Amplitude
   //---------------------------
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel   Z: backward-going
   HELVertex ampzxd(zs, x, dm, gzxd);

   Complex_t amp = ampzxd;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ETRETIBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ETRETIBases.BeamstrahlungFilepath",
                                           "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ETRETIBases.BeamstrahlungFilename","trc500"));
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

  if (!fXBosonPtr)  fXBosonPtr  = new ETRBoson(fMass,fMassDM);
  fXBosonPtr->DebugPrint();
  if (!fDMBosonPtr) fDMBosonPtr = new ETIBoson(fMassDM);
  fDMBosonPtr->DebugPrint();

  Double_t mx = fMass;

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4, mx-10., mx+10.,       50, "Mxr"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthZ" );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiZ"   );
  Xh_init( 7,    70.,   110.,       50, "MZ"     );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthF");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiF"  );
  Xh_init(10,     0.,     2.,        2, "Helin " );
  Xh_init(11,     0.,     2.,        2, "Helot " );
  Xh_init(12,     0.,    12.,       12, "Z mode" );
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ETRETIBases::Userout()
{
  cout << "End of ETRETIBases----------------------------------- "  << endl
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
void ETRETIBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   //                                      dm   f  fb  dm
   static const Int_t kFHelComb[kNf][4] = {{0, -1, +1, 0},
                                           {0, +1, -1, 0}};
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
//  class ETRBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ETRBoson::ETRBoson(Double_t m,
                   Double_t mdm)
        : fMassDM(mdm)
{
   fName    = TString("ETR");
   fPID     = 200000002;
   fCharge  =  0.0;
   fSpin    =  0.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin =  0.5;
   fColor   =  1.0;

   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void ETRBoson::Initialize()
{
   //--
   // Couplings
   //--
   Double_t  a  = kGz/2;
   //--
   // ETR --> ETI + Z
   //--
   GENDecayMode *dmp;
   GENPDTEntry  *d1p, *d2p;
   d1p = new ETIBoson(fMassDM);
   d2p = new GENPDTZBoson();
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
Double_t ETRBoson::GamToSV(Double_t m1, // 1st daughter mass
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
   ANL4DVector pz(sqrt(p*p+m2*m2),0.,0.,p);
   ANL4DVector pdm(sqrt(p*p+m1*m1),0.,0.,-p);
   HELScalar x(px, -1);
   HELScalar dm(pdm, +1);
   HELVector zm(pz, m2, -1, +1);
   HELVector z0(pz, m2,  0, +1);
   HELVector zp(pz, m2, +1, +1);
   HELVertex ampm(zm, x, dm, a);
   HELVertex amp0(z0, x, dm, a);
   HELVertex ampp(zp, x, dm, a);
   double tta = TMath::Power(abs(ampm),2) + TMath::Power(abs(amp0),2)+ TMath::Power(abs(ampp),2);
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
#if 0
   fPID     = 200000001;
#else
   fPID     = 220000; // LSP code for JSFHadronizer.
#endif
   fCharge  =  0.0;
   fSpin    =  0.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin = +0.5;
   fColor   =  1.0;
}
