//*****************************************************************************
//* =====================
//*  WHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> WH generator
//*
//* (Update Record)
//*    2014/01/28  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "WHSpring.h"

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
  fNparton = 4;
  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[2] = bases->fP[1];
  pv[3] = bases->fP[2];
  pv[1] = pv[2] + pv[3];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2WH         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2W          = bases->GetQ2W();
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
  Int_t    cp      = bases->GetCP();               // (W+H-,W-H+)=(+1,-1)
  Int_t    idh     = 37*cp;                        // PDG code for H+
  Int_t    idw     = 24*cp;                        // PDG code for W+
  Int_t    idf3    = bases->f3Ptr->GetPID   ()*cp; // PDG code for f
  Double_t chrg3   = bases->f3Ptr->GetCharge()*cp; // F charge
  Double_t m3      = bases->f3Ptr->GetMass  ();    // F mass
  Int_t    idf4    = bases->f4Ptr->GetPID   ()*cp; // PDG code for f
  Double_t chrg4   = bases->f4Ptr->GetCharge()*cp; // F charge
  Double_t m4      = bases->f4Ptr->GetMass  ();    // F mass
  Int_t    hel1    = bases->fHelFinal[0];          // f1 helicity
  Int_t    hel2    = bases->fHelFinal[1];          // f2 helicity
  Double_t color   = bases->f3Ptr->GetColor();     // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	   // shower level
  Int_t    icf     = 2;                            // color flux id
  Double_t rq2w    = pv[1].Mag();

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

  //                               No. PID  Mass  Charge   pv    Nd 1st  Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, idh, mass,     -cp, *qp[0], 0, 0,  0,    0,   0,     0);
  new (partons[1]) JSFSpringParton(2, idw, rq2w,      cp, *qp[1], 2, 3,  0,    0,   0,     0);
  new (partons[2]) JSFSpringParton(3, idf3,   m3,  chrg3, *qp[2], 0, 0,  2, hel1, icf, islev);
  new (partons[3]) JSFSpringParton(4, idf4,   m4,  chrg4, *qp[3], 0, 0,  2, hel2, icf, islev);

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
         : JSFBases   (name, title), 
           fMass      ( 120.),
           fFhwz      ( 1.),
           fFhwa      ( 0.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fPolp      (0.),
	   fFixCP     ( 1),
           fWModesLo  ( 1),
           fWModesHi  (12),
           fWBosonPtr ( 0),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2WH      (0.),
           fQ2W       (0.),
	   fCP        (1),
           fWModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fXQ2W      (0.),
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

  cout << "Init whbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("WHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.MassH","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Fhwz","1.")); 	       // HWZ form factor
  ins >> fFhwz;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.Fhwa","1.")); 	       // HWA form factor
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
  ins.str(gJSF->Env()->GetValue("WHBases.WModesLo","1"));      // Z decay mode lo
  ins >> fWModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("WHBases.WModesHi","12"));     // Z decay mode hi
  ins >> fWModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fWDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2W          , 0., 1., 0, 1);
  DefineVariable(fCPFinal       , 0., 1., 0, 1);
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
  delete fWBosonPtr;
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
  GENDecayMode *fWModePtr = fWBosonPtr->PickMode(fWDecayMode, weight, fWMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fWModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fWModePtr->At(1));
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
  //  Notice that spin average is taken
  //  care of here
  SelectHelicities(weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  fQ2WH = fEcmIP*fEcmIP;

  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - fMass;
#ifndef __ZEROWIDTH__
  fQ2W = fWBosonPtr->GetQ2BW(qmin, qmax, fXQ2W, weight);
#else
  fQ2W = TMath::Power(fWBosonPtr->GetMass(),2);
  weight = kPi*fWBosonPtr->GetMass()*fWBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch wbranch (fQ2W , fCosThetaF, fPhiF, m3*m3, m4*m4);
  GENBranch cmbranch(fQ2WH, fCosTheta , fPhi , fMass*fMass, &wbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

#ifdef TEMP_H
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  TVector3 boostv(0.,0.,fZBoost/elab);
  ANL4DVector qw  = fP[1] + fP[2];
  qw.Boost(boostv);
  ANL4DVector qcm(fEcmInit,0.,0.,0.);
  ANL4DVector qmm = qcm - qw;

  hMh  ->Fill(qmm.Mag()   , (bsWeight*sigma));
  hRSH ->Fill(fEcmIP      , (bsWeight*sigma));
  hEsum->Fill(eplus+eminus, (bsWeight*sigma));
  hBSweight->Fill(eplus+eminus, (bsWeight));
#endif
  Xh_fill( 1, fEcmIP           , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta        , (bsWeight*sigma));
  Xh_fill( 3, fPhi             , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2W), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaF       , (bsWeight*sigma));
  Xh_fill( 6, fPhiF            , (bsWeight*sigma));
  Xh_fill( 7, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill( 8, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill( 9, (Double_t)fWMode , (bsWeight*sigma));
  Xh_fill(10, (Double_t)(fCP+1)/2, (bsWeight*sigma));

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
  ANL4DVector px = phaseCM.GetFourMomentum(0);
  ANL4DVector pw = phaseCM.GetFourMomentum(1);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;
  fP[0] = px;
  fM[0] = fMass;

  GENBranch &wbranch = *cmbranch.GetBranchPtr(1);
  Double_t cosf = wbranch.GetCosTheta();
  Double_t phif = wbranch.GetPhi     ();
  Double_t m32  = wbranch.GetM12();
  Double_t m42  = wbranch.GetM22();
  GENPhase2 phaseW(pw, m32, m42, cmframe, cosf, phif, 1);
  fP[1] = phaseW.GetFourMomentum(0);
  fP[2] = phaseW.GetFourMomentum(1);
  fM[1] = TMath::Sqrt(m32);
  fM[2] = TMath::Sqrt(m42);
  Double_t betaf = phaseW.GetBetaBar();
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
  ANL4DVector qw = fP[1] + fP[2];
  cerr << " qw = (" 
       << qw.E () << ","
       << qw.Px() << ","
       << qw.Py() << ","
       << qw.Pz() << ")" << endl;
  ANL4DVector pcm = fP[0] + qw;
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
  Double_t  color = f3Ptr->GetColor();
  Int_t     ig1   = f3Ptr->GetGenNo() - 1;
  Int_t     ig2   = f4Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig1][ig2],2);
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
   Double_t gamw   = fWBosonPtr->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  xf(fP[0]);

   Int_t i1 = 1;
   Int_t i2 = 2;
   if (fCP < 0) {
      i1 = 2;
      i2 = 1;
   }
   HELFermion f (fP[i1], fM[i1], fHelFinal [i1-1], +1, kIsOutgoing);
   HELFermion fb(fP[i2], fM[i2], fHelFinal [i2-1], -1, kIsIncoming);
   HELVector  wf(fb, f, glw, grw, kM_w, gamw);

   Complex_t amp = AmpEEtoWH(em, ep, xf, wf);

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

   Double_t  gamw  = fWBosonPtr->GetWidth();
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
  //  Initialize W decay table
  // --------------------------------------------
  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  for (Int_t m=1; m<=fWBosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fWBosonPtr->GetMode(m); 
     if (mp && (m<fWModesLo || m>fWModesHi)) {
        mp->Lock();
     }
  }
  fWBosonPtr->DebugPrint();
  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
  if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
  fZBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh" );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init( 4,    60.,   100.,       50, "Mw"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthF");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF"  );
  Xh_init( 7,     0.,     2.,        2, "Helin ");
  Xh_init( 8,     0.,     2.,        2, "Helot ");
  Xh_init( 9,     0.,    12.,       12, "W mode");
  Xh_init(10,     0.,     2.,        2, "CP    ");
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
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][2] = {{-1, +1}};
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
   weight *= kNf;
}
