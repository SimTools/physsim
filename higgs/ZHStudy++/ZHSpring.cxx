//*****************************************************************************
//* =====================
//*  ZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZH generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZHSpring.h"

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
#define TEMP_H
#ifdef TEMP_H
#if 0
static TH1F *hMh       = 0;
static TH1F *hRSH      = 0;
static TH1F *hEsum     = 0;
static TH1F *hBSweight = 0;
#endif
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZHSpring)
ClassImp(ZHSpringBuf)
ClassImp(ZHBases)

Bool_t ZHBases::fgEnableHtoAA = kFALSE;

//-----------------------------------------------------------------------------
// ==============================
//  class ZHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHSpring::ZHSpring(const char      *name,
                       const char      *title,
                             ZHBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new ZHSpringBuf("ZHSpringBuf",
                              "ZHSpring event buffer",
                              this);
  if (!bases) { 
    SetBases(new ZHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZHSpring::~ZHSpring()
{
  //delete fEventBuf;   // JSFSpring takes care of deleting these
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZHBases *bs = static_cast<ZHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class ZHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZHBases      *bases   = (ZHBases*)((ZHSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  if (!ZHBases::fgEnableHtoAA) {
    fNparton = 4;
  } else {
    fNparton = 10;
  }

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[2] = bases->fP[1];
  pv[3] = bases->fP[2];
  pv[1] = pv[2] + pv[3];

  if (ZHBases::fgEnableHtoAA) {
    GENFrame    cmframe;
    Double_t cosp = gRandom->Uniform(-1.,+1.);
    Double_t phip = gRandom->Uniform(0.,2.*TMath::Pi());
    Double_t mp   = bases->GetMassA();
    Double_t m52  = mp*mp;
    Double_t m62  = m52;
    GENPhase2 phaseH(pv[0], m52, m62, cmframe, cosp, phip, 1);
    pv[4] = phaseH.GetFourMomentum(0);
    pv[5] = phaseH.GetFourMomentum(1);

    Double_t cosa = gRandom->Uniform(-1.,+1.);
    Double_t phia = gRandom->Uniform(0.,2.*TMath::Pi());
    GENPhase2 phaseP1(pv[4], 0., 0., phaseH.GetFrame(), cosa, phia, 1);
    pv[6] = phaseP1.GetFourMomentum(0);
    pv[7] = phaseP1.GetFourMomentum(1);

    cosa = gRandom->Uniform(-1.,+1.);
    phia = gRandom->Uniform(0.,2.*TMath::Pi());
    GENPhase2 phaseP2(pv[5], 0., 0., phaseH.GetFrame(), cosa, phia, 1);
    pv[8] = phaseP2.GetFourMomentum(0);
    pv[9] = phaseP2.GetFourMomentum(1);
  }

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2ZH         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2Z          = bases->GetQ2Z();
  fCosThetaF    = bases->GetCosThetaF();
  fPhiF         = bases->GetPhiF();
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
  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    idp     = 36;                          // PDG code for A
  Int_t    ida     = 22;                          // PDG code for photon
  Int_t    idf     = bases->f3Ptr->GetPID   ();   // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge();   // F charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // F mass
  Double_t m4      = m3;                          // F mass
  Int_t    hel1    = bases->fHelFinal[0];         // f1 helicity
  Int_t    hel2    = bases->fHelFinal[1];         // f2 helicity
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	  // shower level
  Int_t    icf     = 2;                           // color flux id
  Double_t rq2z    = pv[1].Mag();

  Double_t mass    = bases->GetMass();
  Double_t pmass   = bases->GetMassA();
#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idh << " m=" << mass  << " Q=" << 0.
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 2 pid=" << idz << " m=" << rq2z << " Q=" << 0.
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
  TVector qcm(4), qh(4), qz(4);
  qh  = *qp[0];
  qz  = *qp[2] + *qp[3];
  qcm = qh + qz;
  cerr << " ph=(" << qh[0] << ", "
                  << qh[1] << ", "
                  << qh[2] << ", "
                  << qh[3] << ") " << endl;
  cerr << " pz=(" << qz[0] << ", "
                  << qz[1] << ", "
                  << qz[2] << ", "
                  << qz[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  if (!ZHBases::fgEnableHtoAA) {
  //                                No. PID  Mass  Charge   pv    Nd 1st Mom hel  col shower
    new (partons[0]) JSFSpringParton(1, idh, mass,    0., *qp[0], 0, 0,  0,    0,   0,     0);
    new (partons[1]) JSFSpringParton(2, idz, rq2z,    0., *qp[1], 2, 3,  0,    0,   0,     0);
    new (partons[2]) JSFSpringParton(3, idf,   m3,  chrg, *qp[2], 0, 0,  2, hel1, icf, islev);
    new (partons[3]) JSFSpringParton(4,-idf,   m4, -chrg, *qp[3], 0, 0,  2, hel2, icf, islev);
  } else {
  //                                 No. PID  Mass  Charge   pv    Nd 1st Mom hel  col shower
    new (partons[0]) JSFSpringParton( 1, idh, mass,    0., *qp[0], 2, 5,  0,    0,   0,     0);
    new (partons[1]) JSFSpringParton( 2, idz, rq2z,    0., *qp[1], 2, 3,  0,    0,   0,     0);
    new (partons[2]) JSFSpringParton( 3, idf,   m3,  chrg, *qp[2], 0, 0,  2, hel1, icf, islev);
    new (partons[3]) JSFSpringParton( 4,-idf,   m4, -chrg, *qp[3], 0, 0,  2, hel2, icf, islev);
    new (partons[4]) JSFSpringParton( 5, idp,pmass,    0., *qp[4], 2, 7,  1,    0,   0,     0);
    new (partons[5]) JSFSpringParton( 6, idp,pmass,    0., *qp[5], 2, 9,  1,    0,   0,     0);
    new (partons[6]) JSFSpringParton( 7, ida,   0.,    0., *qp[6], 0, 0,  5,    0,   0,     0);
    new (partons[7]) JSFSpringParton( 8, ida,   0.,    0., *qp[7], 0, 0,  5,    0,   0,     0);
    new (partons[8]) JSFSpringParton( 9, ida,   0.,    0., *qp[8], 0, 0,  6,    0,   0,     0);
    new (partons[9]) JSFSpringParton(10, ida,   0.,    0., *qp[9], 0, 0,  6,    0,   0,     0);
  }
  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHBases::ZHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 120.),
           fMassA     (  0.2),
           fLambda    (1000.),
           fA         (   0.),
           fB         (   0.),
           fBtilde    (   0.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2ZH      (0.),
           fQ2Z       (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
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

  cout << "Init zhbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("ZHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.MassH","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.MassA","0.2")); 	 // M_A [GeV]
  ins >> fMassA;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Lambda","1000.")); 	 // Lambda [GeV]
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.A","0.")); 	 // a
  ins >> fA;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.B","0.")); 	 // b
  ins >> fB;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Btilde","0.")); 	 // btilde
  ins >> fBtilde;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.EnableHtoAA","0")); 	 // H --> AA
  Int_t HtoAA;
  ins >> HtoAA;
  if (HtoAA) fgEnableHtoAA = kTRUE;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.BeamWidth","0.002")); // Beam energy spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2Z          , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaF, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiF     , fXL[3], fXU[3], 0, 1);

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
ZHBases::~ZHBases()
{
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZHBases::Func()
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
  fQ2ZH = fEcmIP*fEcmIP;

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

  GENBranch zbranch (fQ2Z , fCosThetaF, fPhiF, m3*m3, m4*m4);
  GENBranch cmbranch(fQ2ZH, fCosTheta , fPhi , fMass*fMass, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

#ifdef TEMP_H
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);
  ANL4DVector qz  = fP[1] + fP[2];
  qz.Boost(boostv);
  ANL4DVector qcm(fEcmInit,0.,0.,0.);
  ANL4DVector qmm = qcm - qz;

#if 1
  H1Fill("hMh"      , qmm.Mag()        , (bsWeight*sigma));
  H1Fill("hCosth"   , fCosTheta        , (bsWeight*sigma));
  H1Fill("hPhi"     , fPhi             , (bsWeight*sigma));
  H1Fill("hCosthF"  , fCosThetaF       , (bsWeight*sigma));
  H1Fill("hPhiF"    , fPhiF            , (bsWeight*sigma));
  H1Fill("hMz"      , TMath::Sqrt(fQ2Z), (bsWeight*sigma));
  H1Fill("hRSH"     , fEcmIP           , (bsWeight*sigma));
  H1Fill("hEsum"    , eplus+eminus     , (bsWeight*sigma));
  H1Fill("hBSweight", eplus+eminus     , (bsWeight));
#endif
#endif
#if 0
  Xh_fill( 1, fEcmIP           , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta        , (bsWeight*sigma));
  Xh_fill( 3, fPhi             , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2Z), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaF       , (bsWeight*sigma));
  Xh_fill( 6, fPhiF            , (bsWeight*sigma));
  Xh_fill( 7, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill( 8, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill( 9, (Double_t)fZMode , (bsWeight*sigma));
  Xh_fill(10, qmm.Mag()        , (bsWeight*sigma));
#endif

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZHBases::DSigmaDX(GENBranch &cmbranch)
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
  fP[0] = px;
  fM[0] = fMass;

  GENBranch &zbranch = * cmbranch.GetBranchPtr(1);
  Double_t cosf = zbranch.GetCosTheta();
  Double_t phif = zbranch.GetPhi     ();
  Double_t m32  = zbranch.GetM12();
  Double_t m42  = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m32, m42, cmframe, cosf, phif, 1);
  fP[1] = phaseZ.GetFourMomentum(0);
  fP[2] = phaseZ.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m32);
  fM[2] = TMath::Sqrt(m42);
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
  ANL4DVector qz = fP[1] + fP[2];
  cerr << " qz = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector pcm = fP[0] + qz;
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
  static const Int_t    kNbr  = 2;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax * betaf;         // phase space factor
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
Double_t ZHBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

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
Complex_t ZHBases::FullAmplitude()
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

   HELScalar  xf(fP[0]);

   HELFermion f (fP[1], fM[1], fHelFinal [0], +1, kIsOutgoing);
   HELFermion fb(fP[2], fM[2], fHelFinal [1], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoZH(em, ep, xf, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZH()
// --------------------------
Complex_t ZHBases::AmpEEtoZH(const HELFermion &em,
                             const HELFermion &ep,
                             const HELScalar  &xf,
                             const HELVector  &zf)
{
   Double_t  qe    = -1.;
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
   // Higgs Production Amplitude
   //---------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Double_t gzzh   = kGz*kM_z;
   Double_t g1     = gzzh + 2 * kM_z * kM_z * (fA/fLambda);
   Double_t g2     = -2 * (fB/fLambda);
   Double_t g3     = -4 * (fBtilde/fLambda);

#ifndef ANOM_ZZH
   Complex_t amp = HELVertex(zs, zf, xf, gzzh);
#else
   Complex_t amp = HELVertex(zs, zf, xf, g1, g2, g3);
#endif
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZHBases::Userin()
{
  TDirectory *last = gDirectory;
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZHBases.BeamstrahlungFilename","trc500"));
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
  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZBosonPtr->GetMode(m); 
     if (mp && (m<fZModesLo || m>fZModesHi)) {
        mp->Lock();
     }
  }
  fZBosonPtr->DebugPrint();

  last->cd();
  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
#if 0
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh" );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init( 4,    70.,   110.,       50, "Mz"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthF");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF"  );
  Xh_init( 7,     0.,     2.,        2, "Helin ");
  Xh_init( 8,     0.,     2.,        2, "Helot ");
  Xh_init( 9,     0.,    12.,       12, "Z mode");
#ifdef TEMP_H
  Xh_init(10, fMass-10., fMass+30.,100, "M_rec ");
#endif
#endif
#ifdef TEMP_H
#if 0
  if (!hMh)   hMh   = new TH1F("hMh"  ,"", 200,115.,135.);
  if (!hRSH)  hRSH  = new TH1F("hRSH" ,"",1100,  0.,fEcmInit*1.1);
  if (!hEsum) hEsum = new TH1F("hEsum","",1100,  0.,fEcmInit*1.1);
  if (!hBSweight) hBSweight = new TH1F("hBSweight","",1100,0,fEcmInit*1.1);
#else
  //H1Init("hMh"      ,"", 200,fMass-10.,    fMass+30.);
  H1Init("hMh"      ,"", 500,     100.,         200.);
  H1Init("hCosth"   ,"", 100,      -1.,          +1.);
  H1Init("hPhi"     ,"", 100,       0.,         k2Pi);
  H1Init("hCosthF"  ,"", 100,      -1.,          +1.);
  H1Init("hPhiF"    ,"", 100,       0.,         k2Pi);
  H1Init("hMz"      ,"", 200,      70.,         110.);
  H1Init("hRSH"     ,"",1100,       0., fEcmInit*1.1);
  H1Init("hEsum"    ,"",1100,       0., fEcmInit*1.1);
  H1Init("hBSweight","",1100,       0., fEcmInit*1.1);
#endif
#endif
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZHBases::Userout()
{
  cout << "End of ZHBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremsstrahlung       = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
#ifdef TEMP_H
#if 0
  hMh  ->Write();
  hRSH ->Write();
  hEsum->Write();
  hBSweight->Write();
#else
#endif
#endif
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void ZHBases::SelectHelicities(Double_t &weight)
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
   fHelFinal  [0] = kFHelComb[fJCombF][0];
   fHelFinal  [1] = kFHelComb[fJCombF][1];
   weight = kNf;
}
