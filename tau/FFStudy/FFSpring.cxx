//*****************************************************************************
//* =====================
//*  FFSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> FF generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "FFSpring.h"

#include "TRandom.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(FFSpring)
ClassImp(FFSpringBuf)
ClassImp(FFBases)

//-----------------------------------------------------------------------------
// ==============================
//  class FFSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
FFSpring::FFSpring(const char      *name,
                       const char      *title,
                             FFBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new FFSpringBuf("FFSpringBuf",
                              "FFSpring event buffer",
                              this);
  if (!bases) { 
    SetBases(new FFBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
FFSpring::~FFSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class FFSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t FFSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  FFBases      *bases   = (FFBases*)((FFSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 2;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0];
  pv[1] = bases->fP[1];

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
  Int_t    idf     = bases->fFPtr->GetPID   ();   // PDG code for f
  Double_t mf      = bases->fFPtr->GetMass  ();   // F mass
  Double_t chrg    = bases->fFPtr->GetCharge();   // F charge
  Double_t color   = bases->fFPtr->GetColor ();   // color factor
  Double_t hf      = bases->fHelFinal[0];         // F helicity
  Double_t hfb     = bases->fHelFinal[1];         // F helicity
  Int_t    islev   = color > 1. ? 201 : 0;     	  // shower level
  Int_t    icf     = 2;                           // color flux id

  //                                No. PID  Mass  Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton(1, idf,   mf,  chrg, *qp[0], 0, 0,  0, hf , icf, islev);
  new (partons[1]) JSFSpringParton(2,-idf,   mf, -chrg, *qp[1], 0, 0,  0, hfb, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class FFBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
FFBases::FFBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fPID       (15),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
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

  cout << "Init ffbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("FFBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("FFBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("FFBases.Ecm","500."));       // E_cm (0.5TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("FFBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("FFBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("FFBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("FFBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("FFBases.BeamstrahlungFilename","trc500"));
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
  SetNoOfSample(2000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 5);
  SetIteration2(0.05, 10);

}
// --------------------------
//  D-tor
// --------------------------
FFBases::~FFBases()
{
  delete fZBosonPtr;
  delete fFPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t FFBases::Func()
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
  Double_t mf   = fFPtr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < 2*mf) {
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

  Double_t sigma = DSigmaDX();

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  H1Fill( "h01", fEcmIP           , (bsWeight*sigma));
  H1Fill( "h02", fCosTheta        , (bsWeight*sigma));
  H1Fill( "h03", fPhi             , (bsWeight*sigma));
  H1Fill( "h04", (Double_t)fJCombI, (bsWeight*sigma));
  H1Fill( "h05", (Double_t)fJCombF, (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t FFBases::DSigmaDX()
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  // fEcmIP, fCosTheta, fPhi, fFPtr);
  Double_t rs   = fEcmIP;
  Double_t cosx = fCosTheta;
  Double_t phix = fPhi;
  Double_t mf   = fFPtr->GetMass();
  Double_t mf2  = mf*mf;

  ANL4DVector qcm(rs, 0., 0., 0.);
  GENFrame    cmframe;

  GENPhase2 phaseCM(qcm, mf2, mf2, cmframe, cosx, phix, 0);
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
Double_t FFBases::AmpSquared()
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
Complex_t FFBases::FullAmplitude()
{
   Double_t mf     = fFPtr->GetMass();

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELFermion f (fP[0], mf, fHelFinal [0], +1, kIsOutgoing);
   HELFermion fb(fP[1], mf, fHelFinal [1], -1, kIsIncoming);

   Complex_t amp = AmpEEtoFF(em, ep, f, fb);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoFF()
// --------------------------
Complex_t FFBases::AmpEEtoFF(const HELFermion &em,
                             const HELFermion &ep,
			     const HELFermion &f,
			     const HELFermion &fb)
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

   //--------------------------------------------------
   // Calculate Amplitudes
   //--------------------------------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Complex_t ampz = HELVertex(f, fb, zs, glzf, grzf);

   HELVector as(em, ep, ge, ge, 0., 0.);
   Complex_t ampa = HELVertex(f, fb, as, gf, gf);

   Complex_t amp = ampz + ampa;

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void FFBases::Userin()
{
  // --------------------------------------------
  //  Initialize Z decay table
  // --------------------------------------------
   if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
   fZBosonPtr->DebugPrint();
   Int_t ig = 2;
   Int_t i3 = 1;
   Int_t lq = 0;
   switch (fPID) {
     case 11:
        ig = 0;
	break;
     case 13:
        ig = 1;
	break;
     case 15:
        ig = 2;
	break;
     default:
        cerr << " PID = " << fPID << " not supported. Quit." << endl;
	abort();
	break;
   }
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
void FFBases::Userout()
{
  cout << "End of FFBases----------------------------------- "  << endl
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
void FFBases::SelectHelicities(Double_t &weight)
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
