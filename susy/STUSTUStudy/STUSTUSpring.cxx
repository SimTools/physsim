//*****************************************************************************
//* =====================
//*  STUSTUSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> STUSTU generator
//*
//* (Update Record)
//*    2011/02/13  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "STUSTUSpring.h"

#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(STUSTUSpring)
ClassImp(STUSTUSpringBuf)
ClassImp(STUSTUBases)

//Bool_t STUSTUBases::fgEnableStauDecay = kTRUE;
Bool_t STUSTUBases::fgEnableStauDecay = kFALSE;

//-----------------------------------------------------------------------------
// ==============================
//  class STUSTUSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
STUSTUSpring::STUSTUSpring(const char      *name,
                           const char      *title,
                               STUSTUBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new STUSTUSpringBuf("STUSTUSpringBuf",
                                  "STUSTUSpring event bustustuer",
                                  this);
  if (!bases) { 
    SetBases(new STUSTUBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
STUSTUSpring::~STUSTUSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class STUSTUSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t STUSTUSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  STUSTUBases  *bases   = (STUSTUBases*)((STUSTUSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  if (!STUSTUBases::fgEnableStauDecay) {
    fNparton = 2;
  } else {
    fNparton = 6;
  }

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[1] = bases->fP[1];

  Double_t mstau = bases->GetMassStau();        // stau mass
  Double_t mtau  = bases->fFPtr->GetMass();     // tau mass
  Double_t mdm   = bases->GetMassDM();          // dark matter mass

  if (STUSTUBases::fgEnableStauDecay) {
    GENFrame    cmframe;
    Double_t cosm = gRandom->Uniform(-1.,+1.);
    Double_t phim = gRandom->Uniform(0.,2.*TMath::Pi());
    Double_t m32  = mdm*mdm;
    Double_t m42  = mtau*mtau;
    GENPhase2 phaseStauM(pv[0], m32, m42, cmframe, cosm, phim, 1);
    pv[2] = phaseStauM.GetFourMomentum(0);
    pv[3] = phaseStauM.GetFourMomentum(1);
    Double_t cosp = gRandom->Uniform(-1.,+1.);
    Double_t phip = gRandom->Uniform(0.,2.*TMath::Pi());
    Double_t m52  = mdm*mdm;
    Double_t m62  = mtau*mtau;
    GENPhase2 phaseStauP(pv[1], m52, m62, cmframe, cosp, phip, 1);
    pv[4] = phaseStauP.GetFourMomentum(0);
    pv[5] = phaseStauP.GetFourMomentum(1);
  }

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
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

  Int_t    idstu   = 1000015;                     // PDG code for stau-
  //Int_t    iddm    = 1000039;                   // PDG code for gravitino
  Int_t    iddm    = 1000022;                     // PDG code for neutralino
  Int_t    idtu    = bases->fFPtr->GetPID();      // PDG code for tau-
  Double_t chrg    = bases->fFPtr->GetCharge();   // stau charge
  Double_t helm    = (1. - bases->fHelTauM)/2.;   // tau- helicity
  Double_t help    = (1. - bases->fHelTauP)/2.;   // tau- helicity
  Double_t htm     = 0;
  Double_t htp     = 0;
  if (gRandom->Uniform() < helm) {
    htm = -1.;
  } else {
    htm = +1.;
  }
  if (gRandom->Uniform() < help) {
    htp = -1.;
  } else {
    htp = +1.;
  }
  Double_t ctau   = bases->fCTau;

  //                                 No. PID   Mass  Charge   pv    Nd 1st Mom hel col shower
  if (!STUSTUBases::fgEnableStauDecay) {
    new (partons[0]) JSFSpringParton(1, idstu, mstau,  chrg, *qp[0], 0, 0,  0,  0., 0, 0);
    new (partons[1]) JSFSpringParton(2,-idstu, mstau, -chrg, *qp[1], 0, 0,  0,  0., 0, 0);
  } else {
    new (partons[0]) JSFSpringParton(1, idstu, mstau,  chrg, *qp[0], 2, 3,  0,  0., 0, 0, ctau);
    new (partons[1]) JSFSpringParton(2,-idstu, mstau, -chrg, *qp[1], 2, 5,  0,  0., 0, 0, ctau);
    new (partons[2]) JSFSpringParton(3,  iddm,   mdm,    0., *qp[3], 0, 0,  1,  0., 0, 0);
    new (partons[3]) JSFSpringParton(4,  idtu,  mtau,  chrg, *qp[2], 0, 0,  1, htm, 0, 0);
    new (partons[4]) JSFSpringParton(5,  iddm,   mdm,    0., *qp[5], 0, 0,  2,  0., 0, 0);
    new (partons[5]) JSFSpringParton(6, -idtu,  mtau, -chrg, *qp[4], 0, 0,  2, htp, 0, 0);
  }

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class STUSTUBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
STUSTUBases::STUSTUBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMassStau  (150.),
           fMassDM    (  1.),
           fThetaMix  (  0.),
           fHelTauM   (+1.),
           fHelTauP   (-1.),
           fCTau      (1.e-2), // [cm]
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fFPtr      (0),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fCosTheta  (0.),
           fPhi       (0.),
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

  cout << "Init stustubases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("STUSTUBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.StauMass","150."));
  ins >> fMassStau;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.DarkMatterMass","1."));
  ins >> fMassDM;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.ThetaMix","0."));
  ins >> fThetaMix;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.HelTauMinus","1."));
  ins >> fHelTauM;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.HelTauPlus","-1."));
  ins >> fHelTauP;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.CTau","1.e-2")); // ctau [cm]
  ins >> fCTau;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.Ecm","500."));       // E_cm (0.5TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.BeamWidth","0.002")); // Beam spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("STUSTUBases.PolElectron","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("STUSTUBases.BeamstrahlungFilepath",
                                  "/proj/soft/data5/samples/gen/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("STUSTUBases.BeamstrahlungFilename","500_nominal"));
    ins >> bsfilename;

    TDirectory *lastDir = gDirectory;
    
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
           << " which is distustuerent from fEcm/2: "    << (fEcmInit/2)
           << endl;        
    } // check the energy homogeneity
    lastDir->cd();
  } // if beamstrahlung is on

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);

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
  SetIteration1(0.05, 5);
  SetIteration2(0.05, 50);

}
// --------------------------
//  D-tor
// --------------------------
STUSTUBases::~STUSTUBases()
{
  delete fZBosonPtr;
  delete fFPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t STUSTUBases::Func()
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
  if (fEcmIP < 2*fMassStau) {
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
  //  Calcuate distustuerential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX();

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  H1Fill( "h01", fEcmIP           , (bsWeight*sigma));
  H1Fill( "h02", fCosTheta        , (bsWeight*sigma));
  H1Fill( "h03", fPhi             , (bsWeight*sigma));
  H1Fill( "h04", (Double_t)fJCombI, (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t STUSTUBases::DSigmaDX()
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  // fEcmIP, fCosTheta, fPhi, fFPtr);
  Double_t rs   = fEcmIP;
  Double_t cosx = fCosTheta;
  Double_t phix = fPhi;
  Double_t mst  = fMassStau;
  Double_t mst2 = mst*mst;

  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  GENPhase2 phaseCM(qcm, mst2, mst2, cmframe, cosx, phix, 0);
  Double_t betax = phaseCM.GetBetaBar();
  if (betax <= 0.) return 0.;
  fP[0] = phaseCM.GetFourMomentum(0);
  fP[1] = phaseCM.GetFourMomentum(1);

  Double_t eb     = rs/2.;
  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

  // --------------------------------------------
  //  Calcuate distustuerential cross section
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
  ANL4DVector pcm = fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
#endif

  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();

  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 1;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax;                 // phase space factor
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
Double_t STUSTUBases::AmpSquared()
{
  Double_t  color = fFPtr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t STUSTUBases::FullAmplitude()
{
   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  stm(fP[0], +1);
   HELScalar  stp(fP[1], +1);

   Complex_t amp = AmpEEtoSTUSTU(em, ep, stm, stp);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoSTUSTU()
// --------------------------
Complex_t STUSTUBases::AmpEEtoSTUSTU(const HELFermion &em,
                                     const HELFermion &ep,
			             const HELScalar  &stm,
			             const HELScalar  &stp)
{
   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);
   Double_t  ge    = qe*kGe;

   Double_t  gamz  = fZBosonPtr->GetWidth();

   Double_t  qf    = fFPtr->GetCharge();
   Double_t  t3f   = fFPtr->GetISpin();
   Double_t  glzf  = -kGz*(t3f - qf*kSin2W);
   Double_t  grzf  = -kGz*(    - qf*kSin2W);
   Double_t  gf    = qf*kGe;
   Double_t  gast  = gf;
   Double_t  gzst  = glzf*TMath::Sin(fThetaMix) + grzf*TMath::Cos(fThetaMix);

   //--------------------------------------------------
   // Calculate Amplitudes
   //--------------------------------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Complex_t ampz = HELVertex(zs, stm, stp, gzst);

   HELVector as(em, ep, ge, ge, 0., 0.);
   Complex_t ampa = HELVertex(as, stm, stp, gast);

   Complex_t amp = ampz + ampa;

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void STUSTUBases::Userin()
{
  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
   if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
   fZBosonPtr->DebugPrint();
   Int_t ig = 2;
   Int_t i3 = 1;
   Int_t lq = 0;
   const Char_t   *nam  = kName[lq][i3][ig];
         Int_t     pid  = kPID [lq][i3][ig];
	 Double_t  qf   = lq ? (i3 ? -1./3. : 2./3.) : (i3 ? -1. : 0.);
	 Double_t  spin = 0.5;
	 Double_t  mf   = kMass[lq][i3][ig];
	 Double_t  t3   = 0.5 - i3;
	 Double_t  cf   = lq ? 3. : 1.; 
   if (!fFPtr) fFPtr = new GENPDTEntry(nam, pid, qf, spin, mf, ig+1,  t3, cf);

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  H1Init( "h01", "Ecm"   , 50,     0., fEcmInit*1.1 );
  H1Init( "h02", "Costh" , 50, fXL[0],        fXU[0]);
  H1Init( "h03", "Phi"   , 50, fXL[1],        fXU[1]);
  H1Init( "h04", "Helin" ,  2,     0.,            2.);
  H1Init( "h05", "Helot" ,  2,     0.,            2.);
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void STUSTUBases::Userout()
{
  cout << "End of STUSTUBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "ostustu")  << endl
       << "Bremsstrahlung       = " << (fISR     ? "on" : "ostustu")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void STUSTUBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   Double_t helm = (1. - fPole)/2.;
   if (fHelCombInitial < helm) {
      fJCombI = 0;
   } else {
      fJCombI = 1;
   }
   fHelInitial[0] = kIHelComb[fJCombI][0];
   fHelInitial[1] = kIHelComb[fJCombI][1];
   weight = 1.;
}
