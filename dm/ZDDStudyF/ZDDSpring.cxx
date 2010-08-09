//*****************************************************************************
//* =====================
//*  ZDDSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZDD generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZDDSpring.h"

#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(ZDDSpring)
ClassImp(ZDDSpringBuf)
ClassImp(ZDDBases)

//-----------------------------------------------------------------------------
// ==============================
//  class ZDDSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZDDSpring::ZDDSpring(const char      *name,
                     const char      *title,
                            ZDDBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new ZDDSpringBuf("ZDDSpringBuf",
                               "ZDDSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new ZDDBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZDDSpring::~ZDDSpring()
{
  delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZDDSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZDDBases *bs = static_cast<ZDDBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZDDBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZDDSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZDDSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZDDBases     *bases   = (ZDDBases*)((ZDDSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 5;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // d1
  pv[1] = bases->fP[1]; // d2
  pv[3] = bases->fP[2]; // f1 from Z
  pv[4] = bases->fP[3]; // f2 from Z
  pv[2] = pv[3] + pv[4];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2ZDD        = fEcmIP*fEcmIP;
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
  Int_t    iddm    = bases->fDMFermionPtr->GetPID(); // PDG code for AH? 
  Int_t    idz     = 23;                          // PDG code for Z
  Int_t    idf     = bases->f3Ptr->GetPID   ();   // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge();   // F charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // F mass
  Double_t m4      = m3;                          // F mass
  Double_t color   = bases->f3Ptr->GetColor();    // color factor
  Int_t    islev   = color > 1. ? 201 : 0;  	  // shower level
  Int_t    icf     = 2;                           // color flux id
  Double_t rq2z    = pv[2].Mag();
  Int_t    hel3    = bases->fHelFinal[2];
  Int_t    hel4    = bases->fHelFinal[3];

  Double_t mass    = bases->GetMass();
  Int_t    hel1    = bases->fHelFinal[0];
  Int_t    hel2    = bases->fHelFinal[1];

  //                                No. PID  Mass  Charge   pv  Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, iddm, mass,    0., *qp[0], 0, 0,  0, hel1,   0,     0);
  new (partons[1]) JSFSpringParton(2, iddm, mass,    0., *qp[1], 0, 0,  0, hel2,   0,     0);
  new (partons[2]) JSFSpringParton(3, idz , rq2z,    0., *qp[2], 2, 4,  0,    0,   0,     0);
  new (partons[3]) JSFSpringParton(4, idf ,   m3,  chrg, *qp[3], 0, 0,  3, hel3, icf, islev);
  new (partons[4]) JSFSpringParton(5,-idf ,   m4, -chrg, *qp[4], 0, 0,  3, hel4, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZDDBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZDDBases::ZDDBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 62.),
           fMassH     ( 120.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fPolp      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
           fCf        (6.86),
           fLambda    (1000.),
           fDMFermionPtr( 0),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2ZDD     (0.),
           fQ2Z       (0.),
           fQ2HH      (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fXQ2Z      (0.),
           fCosThetaF (0.),
           fPhiF      (0.),
           fXQ2HH     (0.),
           fCosThetaH (0.),
           fPhiH      (0.),
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
  stringstream ins(gJSF->Env()->GetValue("ZDDBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.CosthHRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.PhiHOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.MassH","120.")); 	 // M_h [GeV]
  ins >> fMassH;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.MassDM","62.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Cf","6.86")); 	 // Cf
  ins >> fCf;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Lambda","1000.")); 	 // Lambda
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Ecm","500."));        // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.BeamWidth","0.002")); // Beam spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Pole","0."));         // electron polarization
  ins >> fPole;
  cout << "Pole: " << fPole << endl;


  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.Polp","0."));         // positron polarization
  ins >> fPolp;
  cout << "Polp: " << fPolp << endl;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZDDBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fXQ2HH         , 0., 1., 1, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2Z          , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaF, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiF     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaH, fXL[4], fXU[4], 1, 1);
  DefineVariable(fPhiH     , fXL[5], fXU[5], 1, 1);

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
ZDDBases::~ZDDBases()
{
  delete fDMFermionPtr;
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t ZDDBases::Func()
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
  if (fEcmIP < 2*fMass + m3 + m4) {
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
  fQ2ZDD = fEcmIP*fEcmIP;

  // Z
  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - 2*fMass;
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // DD system
  Double_t qhd2min = TMath::Power(2*fMass,2);
  Double_t qhd2max = TMath::Power(rs - TMath::Sqrt(fQ2Z),2);
           fQ2HH   = qhd2min + (qhd2max - qhd2min)*fXQ2HH;
  bsWeight *= qhd2max - qhd2min;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch zbranch (fQ2Z  , fCosThetaF, fPhiF, m3*m3, m4*m4);
  GENBranch hhbranch(fQ2HH , fCosThetaH, fPhiH, fMass*fMass, fMass*fMass);
  GENBranch cmbranch(fQ2ZDD, fCosTheta , fPhi , &hhbranch, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill( 5, fCosThetaF        , (bsWeight*sigma));
  Xh_fill( 6, fPhiF             , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2HH), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaH        , (bsWeight*sigma));
  Xh_fill( 9, fPhiH             , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(12, (Double_t)fZMode  , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZDDBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  Double_t q2hh  = cmbranch.GetM12();
  Double_t q2z   = cmbranch.GetM22();
  Double_t coshh = cmbranch.GetCosTheta();
  Double_t phihh = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2hh, q2z, cmframe, coshh, phihh, 0);
  ANL4DVector phh = phaseCM.GetFourMomentum(0);
  ANL4DVector pz  = phaseCM.GetFourMomentum(1);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  GENBranch &hhbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosh = hhbranch.GetCosTheta();
  Double_t phih = hhbranch.GetPhi     ();
  Double_t mh12 = hhbranch.GetM12();
  Double_t md22 = hhbranch.GetM22();
  GENFrame hhframe = phaseCM.GetFrame();
  GENPhase2 phaseHH(phh, mh12, md22, hhframe, cosh, phih, 1);
  ANL4DVector ph1 = phaseHH.GetFourMomentum(0);
  ANL4DVector pd2 = phaseHH.GetFourMomentum(1);
  fP[0] = ph1;
  fM[0] = TMath::Sqrt(mh12);
  fP[1] = pd2;
  fM[1] = TMath::Sqrt(md22);
  Double_t betah = phaseHH.GetBetaBar();
  if (betah <= 0.) return 0.;

  GENBranch &zbranch = * cmbranch.GetBranchPtr(1);
  Double_t cosf = zbranch.GetCosTheta();
  Double_t phif = zbranch.GetPhi     ();
  Double_t m32  = zbranch.GetM12();
  Double_t m42  = zbranch.GetM22();
  GENPhase2 phaseZ(pz, m32, m42, cmframe, cosf, phif, 1);
  fP[2] = phaseZ.GetFourMomentum(0);
  fP[3] = phaseZ.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
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
  cerr << " fP[3] = (" 
       << fP[3].E () << ","
       << fP[3].Px() << ","
       << fP[3].Py() << ","
       << fP[3].Pz() << ")" << endl;
  ANL4DVector qz = fP[2] + fP[3];
  cerr << " qz = (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector pcm = fP[0] + fP[1] + qz;
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

  Double_t identp = 1./2.;                         // identical particle factor
  Double_t dPhase = kFact * betax * betaf * betah; // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor
  //  Double_t spin   = 1./2.;                         // spin average for e+
  Double_t spin   = 1.;                         // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ZDDBases::AmpSquared(GENBranch &cmbranch)
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
Complex_t ZDDBases::FullAmplitude()
{
   Double_t gamz   = fZBosonPtr->GetWidth();
   Double_t mdm    = fDMFermionPtr->GetMass();

   Double_t qf     = f3Ptr->GetCharge();
   Double_t t3f    = f3Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELFermion d1(fP[0], mdm  , fHelFinal [0], -1, kIsIncoming);  // DM1
   HELFermion d2(fP[1], mdm  , fHelFinal [1], +1, kIsOutgoing);  // DM2

   HELFermion f (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
   HELFermion fb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoZDD(em, ep, d1, d2, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZDD()
// --------------------------
Complex_t ZDDBases::AmpEEtoZDD(const HELFermion &em,
                               const HELFermion &ep,
                               const HELFermion &d1,
                               const HELFermion &d2,
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

   Double_t v     = 2.*kM_w/kGw;
   Double_t gddhl = -fCf * v / fLambda;
   Double_t gddhr = gddhl;

   Double_t gzzh  = kGz*kM_z;

   HELScalar hh(d1, d2, gddhl, gddhr, fMassH, 0.);
   HELVertex amp1(zs, zf, hh, gzzh);         // HHH self-coupling

   Complex_t amp  = amp1;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZDDBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZDDBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZDDBases.BeamstrahlungFilename","trc500"));
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

  if (!fDMFermionPtr) fDMFermionPtr = new DMFermion(fMass);
  fDMFermionPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh" );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init( 4,    70.,   110.,       50, "Mz"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthF");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiF"  );
  Xh_init( 7,     0., fEcmInit,     50, "Mhh"   );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthH");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiH"  );
  Xh_init(10,     0.,     2.,        2, "Helin ");
  Xh_init(11,     0.,     8.,        8, "Helot ");
  Xh_init(12,     0.,    12.,       12, "Z mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZDDBases::Userout()
{
  cout << "End of ZDDBases----------------------------------- "   << endl
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
void ZDDBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 8;
   static const Int_t kFHelComb[kNf][4] = {{-1, -1, -1, +1},
                                           {-1, -1, +1, -1},
                                           {-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1},
                                           {+1, +1, -1, +1},
                                           {+1, +1, +1, -1}};
   Double_t helm = (1. - fPole)/2.;
   if (fHelCombInitial < helm) {
      fJCombI = 0;
      weight = (1. + fPolp)/2.;
   } else {
      fJCombI = 1;
      weight = (1. - fPolp)/2.;
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

//-----------------------------------------------------------------------------
// ==============================
//  class DMFermion
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
DMFermion::DMFermion(Double_t m)
{
   fName    = TString("DMV");
#if 0
   fPID     = 200000001;
#else
   fPID     = 220000; // LSP code for JSFHadronizer.
#endif
   fCharge  =  0.0;
   fSpin    =  1.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin =   .0;
   fColor   =   .0;
}
