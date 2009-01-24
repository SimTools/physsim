//*****************************************************************************
//* =====================
//*  ZHAHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> ZH AH generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "ZHAHSpring.h"

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

ClassImp(ZHAHSpring)
ClassImp(ZHAHSpringBuf)
ClassImp(ZHAHBases)
ClassImp(ZHBoson)
ClassImp(AHBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class ZHAHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHAHSpring::ZHAHSpring(const char      *name,
                       const char      *title,
                               ZHAHBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new ZHAHSpringBuf("ZHAHSpringBuf",
                                "ZHAHSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new ZHAHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
ZHAHSpring::~ZHAHSpring()
{
  delete fEventBuf;
  delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t ZHAHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    ZHAHBases *bs = static_cast<ZHAHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> ZHAHBases written to file" << endl;
  }

  return kTRUE;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHAHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t ZHAHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  ZHAHBases  *bases     = static_cast<ZHAHBases *>
                         (static_cast<ZHAHSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 4;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];    // h  from ZH
  pv[3] = bases->fP[1];    // AH from ZH
  pv[1] = bases->fP[2];    // AH
  pv[0] = pv[2] + pv[3];   // ZH

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2XD         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X          = bases->GetQ2X();
  fCosThetaH    = bases->GetCosThetaH();
  fPhiH         = bases->GetPhiH();
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
  Int_t    idx     = bases->fXBosonPtr ->GetPID(); // PDG code for ZH? 
  Double_t mass    = bases->GetMass();

  Int_t    iddm    = bases->fDMBosonPtr->GetPID(); // PDG code for AH? 
  Double_t msdm    = bases->GetMassDM();
  Double_t mh      = bases->GetMassHiggs();

  Int_t    idh     = 25;                           // PDG code for H

  //                               No. PID   Mass Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton( 1, idx , mass,    0., *qp[0], 2, 3,  0, 0,   0,     0);
  new (partons[1]) JSFSpringParton( 2, iddm, msdm,    0., *qp[1], 0, 0,  0, 0,   0,     0);
  new (partons[2]) JSFSpringParton( 3, idh , mh  ,    0., *qp[2], 0, 0,  1, 0,   0,     0);
  new (partons[3]) JSFSpringParton( 4, iddm, msdm,    0., *qp[3], 0, 0,  1, 0,   0,     0);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHAHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHAHBases::ZHAHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
	   fF         ( 580.),
	   fKappaL    (  0.5),
           fMass      ( 368.),
           fMassDM    ( 81.9),
           fMassT     ( 410.),
           fMassHiggs ( 134.),
           fEcmInit   (1000.),
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
           fQ2XD      (0.),
           fQ2X      (0.),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaH(0.),
           fPhiH     (0.),
           fXQ2X     (0.),
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

  cout << "Init ZHAHBases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("ZHAHBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.CosthHRange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.PhiHOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.F","580.")); 	 // F [GeV]
  ins >> fF;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.KappaL","0.5")); 	 // kappa_l
  ins >> fKappaL;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.MassHiggs","134."));  // M_h [GeV]
  ins >> fMassHiggs;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("ZHAHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fXQ2X          , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();

  DefineVariable(fCosTheta  , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi       , fXL[1], fXU[1], 0, 0);
  DefineVariable(fCosThetaH , fXL[2], fXU[2], 0, 0);
  DefineVariable(fPhiH      , fXL[3], fXU[3], 0, 0);

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
ZHAHBases::~ZHAHBases()
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
Double_t ZHAHBases::Func()
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

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < fMassHiggs + 2*fMassDM) {
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
  fQ2XD = fEcmIP*fEcmIP;

  Double_t rs   = fEcmIP;
  Double_t qmin = fMassHiggs + fMassDM;
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
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch xbranch (fQ2X , fCosThetaH, fPhiH, fMassHiggs*fMassHiggs, fMassDM*fMassDM);
  GENBranch cmbranch(fQ2XD, fCosTheta , fPhi , &xbranch             , fMassDM*fMassDM);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta         , (bsWeight*sigma));
  Xh_fill( 3, fPhi              , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X) , (bsWeight*sigma));
  Xh_fill( 5, fCosThetaH        , (bsWeight*sigma));
  Xh_fill( 6, fPhiH             , (bsWeight*sigma));
  Xh_fill( 7, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill( 8, (Double_t)fJCombF , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t ZHAHBases::DSigmaDX(GENBranch &cmbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t rs = TMath::Sqrt(cmbranch.GetQ2());
  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  // CM -> X DM
  Double_t q2x   = cmbranch.GetM12();
  Double_t q2dm  = cmbranch.GetM22();
  Double_t cosx  = cmbranch.GetCosTheta();
  Double_t phix  = cmbranch.GetPhi     ();
  GENPhase2 phaseCM(qcm, q2x, q2dm, cmframe, cosx, phix, 0);
  ANL4DVector px  = phaseCM.GetFourMomentum(0);
  fP[2]           = phaseCM.GetFourMomentum(1);
  fM[2]           = TMath::Sqrt(q2dm);
  Double_t betax  = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;

  // X -> H DM
  GENBranch &xbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosz  = xbranch.GetCosTheta();
  Double_t phiz  = xbranch.GetPhi     ();
  Double_t mh2   = xbranch.GetM12();
  Double_t mdm2  = xbranch.GetM22();
  GENPhase2 phaseX(px, mh2, mdm2, cmframe, cosz, phiz, 1);
  fP[0]          = phaseX.GetFourMomentum(0);
  fM[0]          = TMath::Sqrt(mh2);
  fP[1]          = phaseX.GetFourMomentum(1);
  fM[1]          = TMath::Sqrt(mdm2);
  Double_t betah = phaseX.GetBetaBar();
  if (betah <= 0.) return 0.;

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
  ANL4DVector qx = fP[0] + fP[1];
  ANL4DVector qh = fP[0];
  cerr << " qh = (" 
       << qh.E () << ","
       << qh.Px() << ","
       << qh.Py() << ","
       << qh.Pz() << ")" << endl;
  cerr << " ---- " << endl;
  cerr << " rs = " << rs << endl;
  cerr << " md1 = " << fP[2].GetMass() << endl;
  cerr << " md2 = " << fP[1].GetMass() << endl;
  cerr << " mh  = " << qh.GetMass() << endl;
  cerr << " mx  = " << qx.GetMass() << endl;

  ANL4DVector pcm = qx + fP[2];
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
  static const Int_t    kNbr  = 2;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
#else
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3+1));
  betah = 1.;
#endif

  Double_t identp = 1.;                      // identical particle factor
  Double_t dPhase = kFact * betax * betah;   // phase space factor
  Double_t flux   = 1./(2.* s * beta_e);     // beam flux factor
  Double_t spin   = 1./2.;                   // spin average for e+

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  return sigma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t ZHAHBases::AmpSquared(GENBranch &cmbranch)
{
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2);

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t ZHAHBases::FullAmplitude()
{
   Double_t mx     = fXBosonPtr ->GetMass();
   Double_t gamx   = fXBosonPtr ->GetWidth();
   Double_t mdm    = fDMBosonPtr->GetMass();

   Double_t Ch     = fXBosonPtr->GetCh();
   Double_t Sh     = fXBosonPtr->GetSh();
   Double_t Cf     = fXBosonPtr->GetCf();
   Double_t Sf     = fXBosonPtr->GetSf();
   Double_t gxhdm  = -kGw*kGw*Sf*Cf*fF/(kCos2W*2.*TMath::Sqrt(2.))
	                  *(Ch*kCosW+Sh*kSinW)*(Sh*kCosW-Ch*kSinW);
#if 0
   cerr << " Ch = " << Ch << " Sh = " << Sh << " Cf = " << Cf << " Sf = " << Sf << endl;
   cerr << " gxhdm = " << gxhdm << " kGw = " << kGw << endl;
#endif

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

#ifndef __NODECAY__
   HELVector  dm (fP[2], mdm  , fHelFinal[2], +1);  // AH
   HELScalar  h  (fP[0], +1);                       // H  from ZH
   HELVector  dmd(fP[1], mdm  , fHelFinal[1], +1);  // AH from ZH
   HELVector  x(dmd, h , gxhdm, mx, gamx);           // ZH
#else
   ANL4DVector px  = fP[0] + fP[1];
   ANL4DVector pdm = fP[2];
   HELVector   x(px  , mx , fHelFinal[1], +1);
   HELVector   dm(pdm, mdm, fHelFinal[2], +1);
#endif
   Complex_t amp = AmpEEtoXD(em, ep, x, dm);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoXD()
// --------------------------
Complex_t ZHAHBases::AmpEEtoXD(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &x,
                               const HELVector  &dm)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t Ch     = fXBosonPtr->GetCh();
   Double_t Sh     = fXBosonPtr->GetSh();
   Double_t  glehzh = kGw*(Ch-Sh*kSinW/(5.*kCosW))/2.;
   Double_t  grehzh = 0.;
   Double_t  glehah = kGw*(Sh+Ch*kSinW/(5.*kCosW))/2.;
   Double_t  grehah = 0.;
#if 0
   static Int_t ncall = 0;
   if (!ncall) {
      ncall = 1;
      cerr << " gehzh = " << glehzh << endl;
      cerr << " gehah = " << glehah << endl;
   }
#endif

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
   // T-channel eH
   //--
   Double_t dummy = 5.; // dummy width for eH
   HELFermion eh1(em, x, glehzh, grehzh, fMassT, dummy);
   HELVertex  amptee1(eh1,ep,dm,glehah,grehah);

   HELFermion eh2(em, dm, glehah, grehah, fMassT, dummy);
   HELVertex  amptee2(eh2,ep,x,glehzh,grehzh);

   Complex_t amp = amptee1 + amptee2;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void ZHAHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("ZHAHBases.BeamstrahlungFilepath",
                                           "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("ZHAHBases.BeamstrahlungFilename","trc500"));
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

  cerr << " ZH Boson "      << endl;
  if (!fXBosonPtr)  fXBosonPtr  = new ZHBoson(fF, fMassHiggs);
  fXBosonPtr->DebugPrint();
  fMass   = fXBosonPtr->GetMass();
  fMassDM = fXBosonPtr->GetMassDM();

  fMassT  = TMath::Sqrt(2.)*fKappaL*fF;
  cerr << " M_ZH = " << fMass   << " [GeV]" << endl
       << " M_AH = " << fMassDM << " [GeV]" << endl
       << " M_eH = " << fMassT  << " [GeV]" << endl;

  cerr << " AH Boson " << endl;
  if (!fDMBosonPtr) fDMBosonPtr = new AHBoson(fMassDM);
  fDMBosonPtr->DebugPrint();

  Double_t mx = fMass;

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"    );
  Xh_init( 2, fXL[0], fXU[0],       50, "Costh"  );
  Xh_init( 3, fXL[1], fXU[1],       50, "Phi"    );
  Xh_init( 4, mx-10., mx+10.,       50, "Mxr"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthH" );
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiH"   );
  Xh_init( 7,     0.,     2.,        2, "Helin " );
  Xh_init( 8,     0.,     9.,        9, "Helot " );
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void ZHAHBases::Userout()
{
  cout << "End of ZHAHBases----------------------------------- "  << endl
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
void ZHAHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 9;
   //                                       h  dm  dm  
   static const Int_t kFHelComb[kNf][3] = {{0, -1, -1},
                                           {0, -1,  0},
                                           {0, -1, +1},
                                           {0,  0, -1},
                                           {0,  0,  0},
                                           {0,  0, +1},
                                           {0, +1, -1},
                                           {0, +1,  0},
                                           {0, +1, +1}};
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
   weight = kNf;
}

//-----------------------------------------------------------------------------
// ==============================
//  class ZHBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
ZHBoson::ZHBoson(Double_t f,
                 Double_t mh)
	: fF(f), fMassHiggs(mh)
{
   fCf          = 1.-4.*kCos2W*kSin2W*kM_z*kM_z/TMath::Power(kGe*f,2);
   fSf          = TMath::Sqrt((1.-fCf)*(1.+fCf));
   Double_t g2  = kGw;
   Double_t gp  = kGe/kCosW;
   Double_t maa = g2*g2*f*f*(fCf*fCf+7.)/8.;
   Double_t mab = g2*gp*f*f*(1.-fCf)*(1.+fCf)/8.;
   Double_t mbb = gp*gp*f*f*(5.*fCf*fCf+3.)/40.;

   fName    = TString("ZH");
   fPID     = 200000003;
   fCharge  =  0.0;
   fSpin    =  1.0;
   fMass    = TMath::Sqrt(0.5*(maa+mbb+TMath::Sqrt(TMath::Power(maa-mbb,2)+4*mab*mab))); 
   fGen     =    0;
   fIsoSpin =  0.0;
   fColor   =  1.0;

   fSh      = -2.*mab/TMath::Sqrt(
                  TMath::Power(maa-mbb+TMath::Sqrt(TMath::Power(maa-mbb,2)+4*mab*mab),2)
                  +4.*mab*mab);
   fCh      = TMath::Sqrt((1.-fSh)*(1.+fSh));
   fMassDM  = TMath::Sqrt(0.5*(maa+mbb-TMath::Sqrt(TMath::Power(maa-mbb,2)+4.*mab*mab)));

   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void ZHBoson::Initialize()
{
   //--
   // Couplings
   //--
   Double_t  a  = -kGw*kGw*fSf*fCf*fF/(kCos2W*2.*TMath::Sqrt(2.))
	                  *(fCh*kCosW+fSh*kSinW)*(fSh*kCosW-fCh*kSinW);
   //--
   // ZH --> H + AH
   //--
   GENDecayMode *dmp;
   GENPDTEntry  *d1p, *d2p;
   Int_t    idh   = 25;
   Double_t chgh  = 0.0;
   Double_t spinh = 0.0;
   d1p = new GENPDTEntry("H",idh,chgh,spinh,fMassHiggs);
   d2p = new AHBoson(fMassDM);
   Double_t m1    = d1p->GetMass();
   Double_t m2    = d2p->GetMass();
   Double_t ident = 1.;
   Double_t gam = GamToSV(m1, m2, a)/ident;
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
Double_t ZHBoson::GamToSV(Double_t m1, // 1st daughter mass
                          Double_t m2, // 2nd daughter mass
                          Double_t a)  // coupling
{
   Double_t x1   = TMath::Power(m1/fMass,2);
   Double_t x2   = TMath::Power(m2/fMass,2);
   Double_t beta = 1. - 2.*(x1+x2) + TMath::Power((x1-x2),2);

   if (beta <= 0.) return 0.;
   beta = TMath::Sqrt(beta);

   ANL4DVector px(fMass,0.,0.,0.);
   double p = fMass*beta/2;
   ANL4DVector ph (sqrt(p*p+m1*m1),0.,0., p);
   ANL4DVector pdm(sqrt(p*p+m2*m2),0.,0.,-p);
   HELScalar h(ph, +1);
   HELVector *xPtr[3], *dPtr[3];
   for (Int_t i=0; i<3; i++) {
     xPtr[i] = new HELVector(px , fMass, i-1, +1);
     dPtr[i] = new HELVector(pdm, m2   , i-1, +1);
   } 
   Double_t amp2 = 0.;
   for (Int_t jx=0; jx<3; jx++) {
      for (Int_t jd=0; jd<3; jd++) {
         HELVertex amp(*xPtr[jx], *dPtr[jd], h, a);
         amp2 += pow(abs(amp),2);
      }
   }
   Double_t spin = 2*fSpin+1;
   Double_t fac  = 1./(16.*kPi)/fMass/spin;
   Double_t gam  = fac*amp2*beta;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class AHBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
AHBoson::AHBoson(Double_t m)
{
   fName    = TString("AH");
#if 0
   fPID     = 200000001;
#else
   fPID     = 220000; // LSP code for JSFHadronizer.
#endif
   fCharge  =  0.0;
   fSpin    =  1.0;
   fMass    =    m;
   fGen     =    0;
   fIsoSpin =  0.0;
   fColor   =  1.0;
}
