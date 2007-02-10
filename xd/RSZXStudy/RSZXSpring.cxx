//*****************************************************************************
//* =====================
//*  RSZXSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- --> ZX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "RSZXSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __HIGGS__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
const Double_t   kSqh     = TMath::Sqrt(0.5); // sqrt(1/2)
const Double_t   kPi      = TMath::Pi();      // Pi
const Double_t   k2Pi     = 2*TMath::Pi();    // 2*Pi
const Double_t   k4Pi     = 4*TMath::Pi();    // 4*Pi
const Double_t   k8Pi     = 8*TMath::Pi();    // 8*Pi
const Double_t   kGeV2fb  = 0.389379292e12;   // GeV to fb
const Double_t   kAlpha   = 1./128.;          // alpha(mz)  = 1/128.
const Double_t   kAlpha0  = 1./137.0359895;   // alpha(q=0) = 1/137.
const Double_t   kAlphaS  = 0.12;             // alphaS(mz) = 0.12
const Double_t   kM_e     = 0.510998902e-3;   // electron mass [GeV]
const Double_t   kM_z     = 91.188;           // Z mass [GeV]
const Double_t   kSin2W   = 0.23;                     // sin^2(theta_W)
const Double_t   kSinW    = TMath::Sqrt(kSin2W);      // sin(theta_W)
const Double_t   kCos2W   = (1. - kSinW)*(1. + kSinW);// cos^2(theta_W)
const Double_t   kCosW    = TMath::Sqrt(kCos2W);      // cos(theta_W)
const Double_t   kSinCosW = kSinW*kCosW;              // sin(2theta_W)/2
const Double_t   kGe      = TMath::Sqrt(k4Pi*kAlpha);  // e
const Double_t   kGw      = kGe/kSinW;                 // gw
const Double_t   kGz      = kGw/kCosW;                 // gz
  
const Char_t   *kName[2][2][3] = {{{"nu_e", "nu_mu"  , "nu_tau"},
                                   {"e"   , "mu"     , "tau"   }},
                                  {{"up"  , "charm"  , "top"   },
                                   {"down", "strange", "bottom"}}};
const Int_t     kPID [2][2][3] = {{{    12,        14,       16},
                                   {    11,        13,       15}},
                                  {{     2,         4,        6},
                                   {     1,         3,        5}}};
const Double_t  kChrg[2][2][3] = {{{    0.,        0.,       0.},
                                   {   -1.,       -1.,      -1.}},
                                  {{  2/3.,      2/3.,     2/3.},
                                   { -1/3.,     -1/3.,    -1/3.}}};

const Double_t  kMass[2][2][3] = {{{0.000000, 0.00000,   0.0000},
                                   {0.511e-3, 0.10566,   1.7770}},
                                  {{0.04    , 1.5    , 175.    },
                                   {0.04    , 0.1    ,   4.7   }}};

ClassImp(TVectorC)

ClassImp(RSZXSpring)
ClassImp(RSZXSpringBuf)
ClassImp(RSZXBases)

ClassImp(HELFermion)
ClassImp(HELVector)
ClassImp(HELScalar)

ClassImp(GENBranch)
ClassImp(GENFrame)
ClassImp(GENPhase2)
ClassImp(GENDecayMode)
ClassImp(GENModePicker)
ClassImp(GENPDTEntry)
ClassImp(GENPDTZBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class RSZXSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSZXSpring::RSZXSpring(const char      *name,
                       const char      *title,
                             RSZXBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new RSZXSpringBuf("RSZXSpringBuf",
                                "RSZXSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new RSZXBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
RSZXSpring::~RSZXSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class RSZXSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t RSZXSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  RSZXBases    *bases   = (RSZXBases*)((RSZXSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  const Int_t kNparton = 6;
  ANL4DVector pv[kNparton];
  pv[2] = bases->fP[0];
  pv[3] = bases->fP[1];
  pv[4] = bases->fP[2];
  pv[5] = bases->fP[3];
  pv[0] = pv[2] + pv[3];
  pv[1] = pv[4] + pv[5];

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fQ2ZX         = fEcmIP*fEcmIP;
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fQ2X          = bases->GetQ2X();
  fCosThetaA    = bases->GetCosThetaA();
  fPhiA         = bases->GetPhiA();
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
#if 0
  Int_t    idx     = 25; 	// PDG code for higgs
#else
  Int_t    idx     = 20000000; 	// PDG code for X??
#endif
  Int_t    ida     = 22; 	// PDG code for gamma
  Int_t    idz     = 23; 	// PDG code for Z
  Int_t    idf     = bases->f3Ptr->GetPID   (); // PDG code for f
  Double_t chrg    = bases->f3Ptr->GetCharge(); // F charge
  Double_t m3      = bases->f3Ptr->GetMass  (); // F mass
  Double_t m4      = m3;                        // F mass
  Double_t color   = bases->f3Ptr->GetColor();  // color factor
  Int_t    islev   = color > 1. ? 101 : 0;  	// shower level
  Int_t    icf     = 1;  	// color flux id
  Double_t rq2z    = pv[1].Mag();

  Double_t mass    = bases->GetMass();
#if 0
//#ifdef __DEBUG__
  cerr << " -------------------------- " << endl;
  cerr << " 1 pid=" << idx << " m=" << mass  << " Q=" << 0.
       << " pv=(" << (*qp[0])(0) << ", "
                  << (*qp[0])(1) << ", "
                  << (*qp[0])(2) << ", "
                  << (*qp[0])(3) << ") " << endl;
  cerr << " 2 pid=" << idz << " m=" << rq2z << " Q=" << 0.
       << " pv=(" << (*qp[1])(0) << ", "
                  << (*qp[1])(1) << ", "
                  << (*qp[1])(2) << ", "
                  << (*qp[1])(3) << ") " << endl;
  cerr << " 3 pid=" << ida << " m=" << 0. << " Q=" << 0.
       << " pv=(" << (*qp[2])(0) << ", "
                  << (*qp[2])(1) << ", "
                  << (*qp[2])(2) << ", "
                  << (*qp[2])(3) << ") " << endl;
  cerr << " 4 pid=" << ida << " m=" << 0. << " Q=" << 0.
       << " pv=(" << (*qp[3])(0) << ", "
                  << (*qp[3])(1) << ", "
                  << (*qp[3])(2) << ", "
                  << (*qp[3])(3) << ") " << endl;
  cerr << " 5 pid=" << idf << " m=" << m3 << " Q=" << chrg
       << " pv=(" << (*qp[4])(0) << ", "
                  << (*qp[4])(1) << ", "
                  << (*qp[4])(2) << ", "
                  << (*qp[4])(3) << ") " << endl;
  cerr << " 6 pid=" <<-idf << " m=" << m4 << " Q=" << -chrg
       << " pv=(" << (*qp[5])(0) << ", "
                  << (*qp[5])(1) << ", "
                  << (*qp[5])(2) << ", "
                  << (*qp[5])(3) << ") " << endl;
  TVector qcm(4), qx(4), qz(4);
  qx  = *qp[2] + *qp[3];
  qz  = *qp[4] + *qp[5];
  qcm = qx +qz;
  cerr << " px=(" << qx[0] << ", "
                  << qx[1] << ", "
                  << qx[2] << ", "
                  << qx[3] << ") " << endl;
  cerr << " pz=(" << qz[0] << ", "
                  << qz[1] << ", "
                  << qz[2] << ", "
                  << qz[3] << ") " << endl;
  cerr << " pcm=(" << qcm[0] << ", "
                   << qcm[1] << ", "
                   << qcm[2] << ", "
                   << qcm[3] << ") " << endl;
#endif

  //                              No. PID  Mass  Charge   pv   Nd 1st Mom hel col shower
  new (partons[0]) JSFSpringParton(1, idx, mass,    0., *qp[0], 2, 3,  0, 0,   0,     0);
  new (partons[1]) JSFSpringParton(2, idz, rq2z,    0., *qp[1], 2, 5,  0, 0,   0,     0);
  new (partons[2]) JSFSpringParton(3, ida,   0.,    0., *qp[2], 0, 0,  1, 0,   0,     0);
  new (partons[3]) JSFSpringParton(4, ida,   0.,    0., *qp[3], 0, 0,  1, 0,   0,     0);
  new (partons[4]) JSFSpringParton(5, idf,   m3,  chrg, *qp[4], 0, 0,  2, 0, icf, islev);
  new (partons[5]) JSFSpringParton(6,-idf,   m4, -chrg, *qp[5], 0, 0,  2, 0, icf, islev);
  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class RSZXBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSZXBases::RSZXBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fLambda    (1000.),
           fC0        (0.),
           fC1        (1.),
           fC2        (1.),
           fC3        (1.),
           fMass      ( 120.),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fPole      (0.),
           fZBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
           fQ2ZX      (0.),
           fQ2X       (0.),
           fQ2Z       (0.),
           f3Ptr      (0),
           f4Ptr      (0),
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaA (0.),
           fPhiA      (0.),
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

  cout << "Init rszxbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("RSZXBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.CosthARange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiAOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.CosthFRange","-1.0 1.0"));
  ins >> fXL[4] >> fXU[4];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.PhiFOverPiRange","0.0 2.0"));
  ins >> fXL[5] >> fXU[5];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Lambda","1000."));    // Lambda [GeV]
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C0","0.")); 	  	 // C_0
  ins >> fC0;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C1","1.")); 		 // C_1
  ins >> fC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C2","1.")); 		 // C_2
  ins >> fC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.C3","1.")); 		 // C_3
  ins >> fC3;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.MassX","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Bremstrahlung","1")); // ISR (on)
  ins >> fISR;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSZXBases.Pole","0."));         // electron polarization
  ins >> fPole;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSZXBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSZXBases.BeamstrahlungFilename","trc500"));
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
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2Z          , 0., 1., 0, 1);
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[2]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();
  fXL[5] = fXL[5]*TMath::Pi();
  fXU[5] = fXU[5]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaA, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiA     , fXL[3], fXU[3], 0, 1);
  DefineVariable(fCosThetaF, fXL[4], fXU[4], 0, 1);
  DefineVariable(fPhiF     , fXL[5], fXU[5], 0, 1);

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
RSZXBases::~RSZXBases()
{
  delete fZBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t RSZXBases::Func()
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
    bsWeight = 1.;
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
  } // if Bremstrahlung


  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung

  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  Double_t weight = 1;
  GENDecayMode *zdkmode = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;

  f3Ptr = static_cast<GENPDTEntry *>(zdkmode->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(zdkmode->At(1));
  Double_t m1   = 0.;
  Double_t m2   = 0.;
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  if (fEcmIP < fMass + m3 + m4) {
    bsWeight = 0.;
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
  fQ2ZX = fEcmIP*fEcmIP;
  fQ2X  = fMass*fMass; // narrow width approx. for X

  Double_t rs   = fEcmIP;
  Double_t qmin = m3 + m4;
  Double_t qmax = rs - fMass;
  fQ2Z = fZBosonPtr->GetQ2BW(qmin, qmax, fXQ2Z, weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  GENBranch xbranch (fQ2X , fCosThetaA, fPhiA, m1, m2);
  GENBranch zbranch (fQ2Z , fCosThetaF, fPhiF, m3, m4);
  GENBranch cmbranch(fQ2ZX, fCosTheta , fPhi , &xbranch, &zbranch);

  Double_t sigma = DSigmaDX(cmbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP           , (bsWeight*sigma));
  Xh_fill( 2, fCosTheta        , (bsWeight*sigma));
  Xh_fill( 3, fPhi             , (bsWeight*sigma));
  Xh_fill( 4, TMath::Sqrt(fQ2X), (bsWeight*sigma));
  Xh_fill( 5, fCosThetaA       , (bsWeight*sigma));
  Xh_fill( 6, fPhiA            , (bsWeight*sigma));
  Xh_fill( 7, TMath::Sqrt(fQ2Z), (bsWeight*sigma));
  Xh_fill( 8, fCosThetaF       , (bsWeight*sigma));
  Xh_fill( 9, fPhiF            , (bsWeight*sigma));
  Xh_fill(10, (Double_t)fJCombI, (bsWeight*sigma));
  Xh_fill(11, (Double_t)fJCombF, (bsWeight*sigma));
  Xh_fill(12, (Double_t)fZMode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t RSZXBases::DSigmaDX(GENBranch &cmbranch)
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

  GENBranch &xbranch = * cmbranch.GetBranchPtr(0);
  Double_t cosa = xbranch.GetCosTheta();
  Double_t phia = xbranch.GetPhi     ();
  Double_t m12  = xbranch.GetM12();
  Double_t m22  = xbranch.GetM22();
  GENPhase2 phaseX(px, m12, m22, cmframe, cosa, phia, 1);
  fP[0] = phaseX.GetFourMomentum(0);
  fP[1] = phaseX.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betaa = phaseX.GetBetaBar();

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

  Double_t eb     = rs/2.;
  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

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
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k8Pi*k4Pi*k2Pi,kNbr));

  Double_t identp = 1.;                            // identical particle factor
  Double_t dPhase = kFact * betax * betaa * betaf; // phase space factor
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
Double_t RSZXBases::AmpSquared(GENBranch &cmbranch)
{
  Double_t  color = f3Ptr->GetColor();
  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;
#if 1
  amp2 *= k2Pi * k8Pi; // branch factor correction: no phase space for X
#endif

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
Complex_t RSZXBases::FullAmplitude()
{
   Double_t   gamz    = fZBosonPtr->GetWidth();

   Double_t qf     = f3Ptr->GetCharge();
   Double_t t3f    = f3Ptr->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELScalar  xf(fP[0]+fP[1], kFALSE);

   HELFermion f (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
   HELFermion fb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming);
   HELVector  zf(fb, f, glz, grz, kM_z, gamz);

   Complex_t amp = AmpEEtoZX(em, ep, xf, zf);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoZX()
// --------------------------
Complex_t RSZXBases::AmpEEtoZX(const HELFermion &em,
                               const HELFermion &ep,
                               const HELScalar  &xf,
                               const HELVector  &zf)
{
   //-------------------
   // Coupling consts.
   //-------------------
   Double_t  a1    = fC1*kCos2W + fC2*kSin2W;
   Double_t  a2    = (fC1 - fC2) * kSinCosW;
   Double_t  a3    = fC1*kSin2W + fC2*kCos2W;
   Double_t  gaax  = a1*kSqh/fLambda;
   Double_t  gazx  = a2*kSqh/fLambda;
   Double_t  gzzx  = a3*kSqh/fLambda;

   Double_t  qe    = -1.;
   Double_t  ge    = -qe*kGe;
   Double_t  glae  = ge;
   Double_t  grae  = ge;

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
#ifdef __HIGGS__
   //---------------------------
   // Higgs Production Amplitude
   //---------------------------
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   Double_t gzzh  = kGz*kM_z;
   Complex_t amp  = zs[0]*zf[0] - zs[1]*zf[1] - zs[2]*zf[2] - zs[3]*zf[3];
             amp *= gzzh;
#else
   //---------------------------
   // ZX Production Amplitude
   //---------------------------
#if 0
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz);
   ANL4DVector k1b  = zs.GetFourMomentum();      // s-channel   Z: backward-going
   ANL4DVector k1   = ANL4DVector(-k1b.E(),      // s-channel   Z: forward-going
                                  -k1b.Px(), 
                                  -k1b.Py(),
                                  -k1b.Pz());
   ANL4DVector k2   = zf.GetFourMomentum();      // final state Z: forward-going 

   Double_t    k1k2 = k1*k2;
   Complex_t   k1zf = k1(0)*zf[0] - k1(1)*zf[1] - k1(2)*zf[2] - k1(3)*zf[3];
   Complex_t   zszf = conj(zs[0])*zf[0] - conj(zs[1])*zf[1]
                                        - conj(zs[2])*zf[2]
                                        - conj(zs[3])*zf[3];
   Complex_t   k2zs = conj(zs[0])*k2(0) - conj(zs[1])*k2(1)
                                        - conj(zs[2])*k2(2)
                                        - conj(zs[3])*k2(3);
   Complex_t   ampzzx = gzzx * (k1k2*zszf - k1zf*k2zs);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma

   Complex_t   aszf = conj(as[0])*zf[0] - conj(as[1])*zf[1]
                                        - conj(as[2])*zf[2]
                                        - conj(as[3])*zf[3];
   Complex_t   k2as = conj(as[0])*k2(0) - conj(as[1])*k2(1)
                                        - conj(as[2])*k2(2)
                                        - conj(as[3])*k2(3);
   Complex_t   ampazx = gazx * (k1k2*aszf - k1zf*k2as);
#else
   //--
   // S-channel Z
   //--
   HELVector zs(em, ep, glze, grze, kM_z, gamz); // s-channel   Z: backward-going
   ANL4DVector k1   = zs.GetFourMomentum();      // s-channel   Z: backward-going
   ANL4DVector k2   = zf.GetFourMomentum();      // final state Z: forward-going 

   Double_t    k1k2 = k1*k2;
   Complex_t   k1zf = k1(0)*zf[0] - k1(1)*zf[1] - k1(2)*zf[2] - k1(3)*zf[3];
   Complex_t   zszf = zs[0]*zf[0] - zs[1]*zf[1] - zs[2]*zf[2] - zs[3]*zf[3];
   Complex_t   k2zs = zs[0]*k2(0) - zs[1]*k2(1) - zs[2]*k2(2) - zs[3]*k2(3);
   Complex_t   ampzzx = gzzx * (k1k2*zszf - k1zf*k2zs);

   //--
   // S-channel A
   //--
   HELVector as(em, ep, glae, grae,   0.,   0.); // s-channel gamma

   Complex_t   aszf = as[0]*zf[0] - as[1]*zf[1] - as[2]*zf[2] - as[3]*zf[3];
   Complex_t   k2as = as[0]*k2(0) - as[1]*k2(1) - as[2]*k2(2) - as[3]*k2(3);
   Complex_t   ampazx = gazx * (k1k2*aszf - k1zf*k2as);
#endif

   //--
   // Sum of the two
   //--
   Complex_t amp = ampzzx + ampazx;
#endif /* end __HIGGS__ */
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void RSZXBases::Userin()
{
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
  Xh_init( 4,   100.,   140.,       50, "Mx"    );
  Xh_init( 5, fXL[2], fXU[2],       50, "CosthA");
  Xh_init( 6, fXL[3], fXU[3],       50, "PhiA"  );
  Xh_init( 7,    70.,   110.,       50, "Mz"    );
  Xh_init( 8, fXL[4], fXU[4],       50, "CosthF");
  Xh_init( 9, fXL[5], fXU[5],       50, "PhiF"  );
  Xh_init(10,     0.,     2.,        2, "Helin ");
  Xh_init(11,     0.,     2.,        2, "Helot ");
  Xh_init(12,     0.,    12.,       12, "Zdecay");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void RSZXBases::Userout()
{
  cout << "End of RSZXBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Bremstrahlung        = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void RSZXBases::SelectHelicities(Double_t &weight)
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
   fJCombF = TMath::Min(fJCombF, 1);
   fHelFinal  [0] = 0;
   fHelFinal  [1] = 0;
   fHelFinal  [2] = kFHelComb[fJCombF][0];
   fHelFinal  [3] = kFHelComb[fJCombF][1];
   weight = kNf;
}


//-----------------------------------------------------------------------------
// ==============================
//  class HELFermion
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------------------
// IXXXXX() and OXXXXX()
//----------------------
HELFermion::HELFermion(const ANL4DVector &p,
                             Double_t     m,
                             Int_t        hel,
                             Int_t        nsf,
                             Bool_t       isincom)
          : TVectorC(4),
            fP(nsf*p.E(), nsf*p.Px(), nsf*p.Py(), nsf*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSF(nsf),
            fIsIncoming(isincom)
{
   Double_t  sf[2], sfomeg[2], omega[2], pp, pp3, sqpop3, sqm;
   Complex_t chi[2];
   Int_t nh = hel*nsf;
   Int_t in = isincom ? 1 : -1;
   if (m == 0.) {
      sqpop3 = nsf*TMath::Sqrt(TMath::Max(p(0)+p(3), 0.));
      chi[0] = sqpop3;
      if (sqpop3 == 0.) chi[1] = -hel * TMath::Sqrt(2.*p(0));
      else              chi[1] = Complex_t(nh*p(1), in*p(2))/sqpop3;
      Int_t iu = (1-in)/2;
      Int_t id = (1+in)/2;
      if (nh == 1*in) {
         (*this)[0] = 0.;
         (*this)[1] = 0.;
         (*this)[2] = chi[iu];
         (*this)[3] = chi[id];
      } else {
         (*this)[0] = chi[id];
         (*this)[1] = chi[iu];
         (*this)[2] = 0.;
         (*this)[3] = 0.;
      }
   } else {
      pp = TMath::Min(p(0), p.Vect().Mag());
      if (pp == 0.) {
         sqm = TMath::Sqrt(m);
         Int_t ip =    (1+in*nh)/2;
         Int_t im = in*(1-in*nh)/2;
         (*this)[0] = ip       * sqm;
         (*this)[1] = im * nsf * sqm;
         (*this)[2] = ip * nsf * sqm;
         (*this)[3] = im       * sqm;
      } else {
         sf[0] = (1+nsf+(1-nsf)*nh)*0.5;
         sf[1] = (1+nsf-(1-nsf)*nh)*0.5;
         omega[0] = TMath::Sqrt(p(0)+pp);
         omega[1] = m/omega[0];
         Int_t ip = (1+nh)/2;
         Int_t im = (1-nh)/2;
         sfomeg[0] = sf[0] * omega[ip];
         sfomeg[1] = sf[1] * omega[im];
         pp3    = TMath::Max(pp+p(3), 0.);
         chi[0] = TMath::Sqrt(pp3*0.5/pp);
         if (pp3 == 0.) chi[1] = -nh;
         else           chi[1] = Complex_t(nh*p(1), in*p(2))/TMath::Sqrt(2.*pp*pp3);
         Int_t iu = (1-in)/2;
         Int_t id = (1+in)/2;
         (*this)[0] = sfomeg[iu] * chi[im];
         (*this)[1] = sfomeg[iu] * chi[ip];
         (*this)[2] = sfomeg[id] * chi[im];
         (*this)[3] = sfomeg[id] * chi[ip];
      }
   }
#ifdef __DEBUG__ 
	   cerr << (isincom ? "fin" : "fot")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
#endif
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELVector
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// VXXXXX()
//----------
HELVector::HELVector(const ANL4DVector &p,
                           Double_t     m,
                           Int_t        hel,
                           Int_t        nsv)
          : TVectorC(4),
            fP(nsv*p.E(), nsv*p.Px(), nsv*p.Py(), nsv*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSV(nsv)
{
   Double_t  hel0, pt, pt2, pp, pzpt, emp;
   Int_t     nsvahl;
   nsvahl = nsv*TMath::Abs(hel);
   pt2 = p.GetPt2();
   pp  = TMath::Min(p(0), p.Vect().Mag());
   pt  = TMath::Min(pp, TMath::Sqrt(pt2));
   if (m == 0.) {
      pp = p(0);
      pt = p.GetPt();
      (*this)[0] = 0.;
      (*this)[3] = hel*pt/pp*kSqh;
      if (pt != 0.) {
         pzpt = p(3)/(pp*pt)*kSqh*hel;
         (*this)[1] = Complex_t(-p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
         (*this)[2] = Complex_t(-p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
      } else {
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
      }
   } else {
      hel0 = 1. - TMath::Abs(hel);
      if (pp == 0.) {
         (*this)[0] = 0.;
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*kSqh);
         (*this)[3] = hel0;
      } else {
         emp = p(0)/(m*pp);
	 (*this)[0] = hel0*pp/m;
         (*this)[3] = hel0*p(3)*emp + hel*pt/pp*kSqh;
	 if (pt != 0.) {
            pzpt = p(3)/(pp*pt)*kSqh*hel;
	    (*this)[1] = Complex_t(hel0*p(1)*emp - p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
	    (*this)[2] = Complex_t(hel0*p(2)*emp - p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
	 } else {
	    (*this)[1] = -hel*kSqh;
	    (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
	 }
      }
   }
}

//----------
// JIOXXX()
//----------
HELVector::HELVector(const HELFermion &fin,
                     const HELFermion &fout,
                           Double_t    glv,
                           Double_t    grv,
                           Double_t    mv,
                           Double_t    gamv)
          : TVectorC(4),
            fP(fout.fP - fin.fP),
	    fM(mv)
{
   Complex_t c0, c1, c2, c3, cs, d;
   Double_t  q2, vm2, dd;
   q2  = fP.Mag2();
   vm2 = mv*mv;
   if (mv == 0.) {
      dd = 1./q2;
      if (grv == 0.) { // purely left-handed
         dd *= glv;
	 (*this)[0] = ( fout[2]*fin[0] + fout[3]*fin[1]) * dd;
	 (*this)[1] = (-fout[2]*fin[1] - fout[3]*fin[0]) * dd;
	 (*this)[2] = ( fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., dd);
	 (*this)[3] = (-fout[2]*fin[0] + fout[3]*fin[1]) * dd;
      } else {
         (*this)[0] = (  glv * ( fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] + fout[1]*fin[3])) * dd;
         (*this)[1] = (- glv * ( fout[2]*fin[1] + fout[3]*fin[0])
                       + grv * ( fout[0]*fin[3] + fout[1]*fin[2])) * dd;
         (*this)[2] = (  glv * ( fout[2]*fin[1] - fout[3]*fin[0])
                       + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., dd);
         (*this)[3] = (  glv * (-fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] - fout[1]*fin[3])) * dd;
      }
   } else {
      d = 1./Complex_t(q2 - vm2, TMath::Max(TMath::Sign(mv*gamv, q2), 0.));
      if (grv == 0.) { // purely left-handed
         d *= glv;
	 c0 =  fout[2]*fin[0] + fout[3]*fin[1];
	 c1 = -fout[2]*fin[1] - fout[3]*fin[0];
	 c2 = (fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., 1.);
	 c3 = -fout[2]*fin[0] + fout[3]*fin[1];
      } else {
         c0 =   glv * ( fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] + fout[1]*fin[3]);
         c1 = - glv * ( fout[2]*fin[1] + fout[3]*fin[0])
              + grv * ( fout[0]*fin[3] + fout[1]*fin[2]);
         c2 =  (glv * ( fout[2]*fin[1] - fout[3]*fin[0])
              + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., 1.);
         c3 =   glv * (-fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] - fout[1]*fin[3]);
      }
      cs = (fP(0)*c0 - fP(1)*c1 - fP(2)*c2 - fP(3)*c3)/vm2;
      (*this)[0] = (c0 - cs*fP(0))*d;
      (*this)[1] = (c1 - cs*fP(1))*d;
      (*this)[2] = (c2 - cs*fP(2))*d;
      (*this)[3] = (c3 - cs*fP(3))*d;
   }
#ifdef __DEBUG__ 
	   cerr << " jio =(" << (*this)[0] << ", "
                             << (*this)[1] << ", "
                             << (*this)[2] << ", "
                             << (*this)[3] << ") " << endl;
#endif
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELScalar
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// SXXXXX()
//----------
HELScalar::HELScalar(const ANL4DVector &p,
                           Int_t        nss)
          : Complex_t(1.,0.),
            fP(nss*p.E(), nss*p.Px(), nss*p.Py(), nss*p.Pz()),
            fNSS(nss)
{
}


//-----------------------------------------------------------------------------
// ==============================
//  class GENDecayMode
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  DebugPrint
// --------------------------
void GENDecayMode::DebugPrint(const Option_t *)
{
   using namespace std;
   cerr << " Gamma = " << setw(6) << setprecision(4) 
                       << setiosflags(ios::fixed)
                       << setiosflags(ios::showpoint)
        << GetGamma() 
        << " BR = "    << setw(6) << setprecision(4)
                       << setiosflags(ios::fixed)
                       << setiosflags(ios::showpoint)
        << GetBR()     <<  " : --> ";
   TIter next(this);
   GENPDTEntry *ep;
   while ((ep = static_cast<GENPDTEntry *>(next()))) {
      cerr <<  ep->GetName() << (ep->GetPID() > 0 ? " " : "bar ");
   }
   cerr << endl;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENModePicker
// ==============================
//_____________________________________________________________________________
// --------------------------
//  Add
// --------------------------
void GENModePicker::Add(GENDecayMode *mp)
{
   TObjArray::Add(mp);
   fDone = kFALSE;
}

//_____________________________________________________________________________
// --------------------------
//  PickMode
// --------------------------
GENDecayMode *GENModePicker::PickMode(Double_t x,
                                      Double_t &weight,
                                      Int_t    &mode)
{
   if (!fDone) Update();
   mode = 0;
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      if (x <= mp->fCumBR) {
         weight = 1./mp->fBR; 
         break;
      }
      mode++;
   }
   return static_cast<GENDecayMode *>(mp);
}

//_____________________________________________________________________________
// --------------------------
//  Update
// --------------------------
void GENModePicker::Update()
{
   if (fDone) return;
   fDone = kTRUE;

   fGamma = 0.;
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      fGamma += mp->fGamma;
   }

   Double_t cum = 0.;
   next.Reset();
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      mp->fBR = mp->fGamma/fGamma;
      cum    += mp->fBR;
      mp->fCumBR = cum;
   }
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTEntry
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTEntry::GENPDTEntry(const Char_t     *name,
                               Int_t       pid,
                               Double_t    charge,
                               Double_t    spin,
                               Double_t    mass,
                               Int_t       gen,
                               Double_t    ispin,
                               Double_t    color)
           : fName(name),
             fPID(pid),
             fCharge(charge),
             fSpin(spin),
             fMass(mass),
             fGen(gen),
             fIsoSpin(ispin),
             fColor(color)
{
}

//_____________________________________________________________________________
// --------------------------
//  GetQ2BW
// --------------------------
Double_t GENPDTEntry::GetQ2BW(Double_t    qmin,   // Q_min
                              Double_t    qmax,   // Q_max
                              Double_t       x,   // integration variable
                              Double_t &weight)   // Jacobian weight
{
   Double_t m    = GetMass ();
   Double_t gm   = GetWidth();
   Double_t mgm  = m*gm;
   Double_t m2   = m*fMass;
   Double_t mgm2 = mgm*mgm;

   Double_t thmin = TMath::ATan((qmin-m)*(qmin+m)/mgm);
   Double_t thmax = TMath::ATan((qmax-m)*(qmax+m)/mgm);
   Double_t theta = thmin + (thmax - thmin)*x;
   Double_t q2    = mgm * TMath::Tan(theta) + m2;

   weight = (thmax - thmin)*(TMath::Power(q2-m2,2) + mgm2)/mgm;

   return q2;
}

//_____________________________________________________________________________
// --------------------------
//  DebugPrint
// --------------------------
void GENPDTEntry::DebugPrint(const Option_t *opt)
{
   using namespace std;
   cerr << " ---------------------------------------------------- " << endl;
   Update();
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      mp->DebugPrint(opt);
   }
   cerr << " ---------------------------------------------------- " << endl
        << " Gamma_tot = " << GetWidth() << " [GeV]"                << endl;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTZBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTZBoson::GENPDTZBoson()
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTZBoson::~GENPDTZBoson()
{
   TIter next(this);
   GENDecayMode *dmp;
   while ((dmp = static_cast<GENDecayMode *>(next()))) {
      dmp->Delete();
   }
   Delete();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENPDTZBoson::Initialize()
{
   fName   = TString("Z");
   fPID    = 23;
   fCharge = 0.;
   fSpin   = 1.;
   fMass   = kM_z;
   //--
   //  Set decay modes
   //--
   for (Int_t ic=0; ic<2; ic++) {
      for (Int_t ig=0; ig<3; ig++) {
         for (Int_t it=0; it<2; it++) { 
            const Char_t   *name = kName[ic][it][ig];
                  Int_t     pid  = kPID [ic][it][ig];  
                  Double_t  t3   = (1 - 2*it)/2.;
                  Double_t  qf   = kChrg[ic][it][ig];
                  Double_t  spin = 0.5;
                  Double_t  cf   = 2*ic + 1;
                  Double_t  mass = kMass[ic][it][ig];
            GENDecayMode *dmp;
            GENPDTEntry  *d1p, *d2p;
            if (ic) cf  *= 1 + kAlphaS/kPi;
            d1p  = new GENPDTEntry(name, pid, qf, spin, mass, ig+1, t3, cf);
            d2p  = new GENPDTEntry(name,-pid,-qf, spin, mass, ig+1, t3, cf);
            Double_t gam  = GamToFF(t3, qf, cf, mass);
            if (gam == 0.) continue;
            dmp = new GENDecayMode(gam);
            dmp->Add(d1p);
            dmp->Add(d2p);
            Add(dmp); 
         }
      }
   }
}

//_____________________________________________________________________________
// --------------------------
//  GamToFF
// --------------------------
Double_t GENPDTZBoson::GamToFF(Double_t t3, // weak isospin
                               Double_t qf, // charge
                               Double_t cf, // color factor
                               Double_t m)  // mass
{
   Double_t mz2 = fMass*fMass;
   Double_t p1  = mz2/4 - m*m;

   if (p1 <= 0.) return 0.;

   p1 = TMath::Sqrt(p1);
   Double_t gv  =  t3/2 - kSin2W*qf;
   Double_t ga  = -t3/2;
   Double_t tta = 2*((gv*gv + ga*ga)*(mz2 - 4*p1*p1/3) - 4*ga*ga*m*m); 
   Double_t fac = (kAlpha/kSin2W/kCos2W)/2;
   Double_t gam = fac*tta*cf*p1/mz2;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENBranch
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENBranch::GENBranch(Double_t q2,
                     Double_t costh,
                     Double_t phi,
                     Double_t m12,
                     Double_t m22)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(m12),
           fM22(m22),
           fBR1Ptr(0),
           fBR2Ptr(0)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

GENBranch::GENBranch(Double_t   q2,
                     Double_t   costh,
                     Double_t   phi,
                     GENBranch *br1p,
                     GENBranch *br2p)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(br1p->GetQ2()),
           fM22(br2p->GetQ2()),
           fBR1Ptr(br1p),
           fBR2Ptr(br2p)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

// ==============================
//  class GENFrame
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENFrame::GENFrame()
{
   fEV[0].SetXYZ(1., 0., 0.);
   fEV[1].SetXYZ(0., 1., 0.);
   fEV[2].SetXYZ(0., 0., 1.);
}

GENFrame::GENFrame(const ANL4DVector &q, const GENFrame &eb)
{
   fEV[2] = q.Vect().Unit();
   fEV[1] = eb.fEV[2].Cross(fEV[2]);
   Double_t ae2  = fEV[1].Mag();
   static const Double_t kXmin = 1.e-12;
   if (ae2 < kXmin) {
      fEV[0] = eb.fEV[0]; fEV[1] = eb.fEV[1]; fEV[2] = eb.fEV[2];
      Double_t csth = fEV[2] * eb.fEV[2];
      if (csth <= 0.) {
         fEV[2] = -eb.fEV[2];
      }
      return;
   } else {
      fEV[1] = fEV[1].Unit();
   }
   fEV[0] = fEV[1].Cross(fEV[2]).Unit();
}

//_____________________________________________________________________________
// --------------------------
//  Transform
// --------------------------
ANL4DVector GENFrame::Transform(const ANL4DVector &pb)
{
   TVector3 pb3v = pb.X()*fEV[0] + pb.Y()*fEV[1] + pb.Z()*fEV[2];
   return ANL4DVector(pb.E(),pb3v.Px(), pb3v.Py(), pb3v.Pz());
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPhase2
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPhase2::GENPhase2(const ANL4DVector &q,
                           Double_t       m12,
                           Double_t       m22,
                     const GENFrame      &eb,
                           Double_t       costh,
                           Double_t       phi,
                           Int_t          mode)
         : fQ(q), fM12(m12), fM22(m22), 
           fEb(eb), fEa(eb), 
           fCosTheta(costh), fPhi(phi),
           fBetaBar(0.), 
           fMode(mode),
           fDone(kFALSE)
{
}

//_____________________________________________________________________________
// --------------------------
//  GetFourMomentum
// --------------------------
ANL4DVector GENPhase2::GetFourMomentum(Int_t i)
{
   if (!fDone) Update(); return i ? fP2 : fP1;
}

//_____________________________________________________________________________
// --------------------------
//  GetFrame
// --------------------------
GENFrame GENPhase2::GetFrame(Int_t i)
{
   if (!fDone) Update(); return i ? fEa : fEb;
}

//_____________________________________________________________________________
// --------------------------
//  GetBetaBar
// --------------------------
Double_t GENPhase2::GetBetaBar()
{
   if (!fDone) Update(); return fBetaBar;
}

//_____________________________________________________________________________
// --------------------------
//  Beta2
// --------------------------
Double_t GENPhase2::Beta2(Double_t x1, Double_t x2)
{
   return 1. - 2*(x1+x2) + (x1-x2)*(x1-x2);
}

//_____________________________________________________________________________
// --------------------------
//  Update
// --------------------------
void GENPhase2::Update()
{
   // -------------------------
   //  Check if update needed
   // -------------------------
   if (fDone) return;
   fDone = kTRUE;

   // -------------------------
   //  Calculate Beta_bar
   // -------------------------
   Double_t amq2 = fQ.Mag2();
   if (amq2 <= 0.) {
      cerr << " >>>>> Error in GENPhase2::GENPhase2() >>>>>>>> " << endl
           << " q = ("  << fQ.E() << ", " 
                        << fQ.X() << ", " 
                        << fQ.Y() << ", "
                        << fQ.Z() << ")" << endl
           << " q2  = " << amq2   
           << " m12 = " << fM12
           << " m22 = " << fM22 << endl;
      fBetaBar = 0.;
      return;
   }

   Double_t amq = TMath::Sqrt(amq2);
   fBetaBar = Beta2(fM12/amq2, fM22/amq2);   
   if (fBetaBar < 0.) {
      fBetaBar = 0.;
      return;
   }
   fBetaBar = TMath::Sqrt(fBetaBar);

   // -------------------------
   //  Daughter momenta
   // -------------------------
   Double_t ap1  = (amq/2) * fBetaBar;
   Double_t snth = TMath::Sqrt((1.-fCosTheta)*(1.+fCosTheta));
   fP1.SetXYZT(ap1*snth*TMath::Cos(fPhi),
               ap1*snth*TMath::Sin(fPhi),
               ap1*fCosTheta,
               TMath::Sqrt(ap1*ap1+fM12));
   fP2.SetXYZT(-fP1.X(), -fP1.Y(), -fP1.Z(), amq - fP1.E());

   // -------------------------
   //  Boost them to lab. frame
   // -------------------------
   if (!fMode) {
      fEa = fEb;
   } else {
      fEa = GENFrame(fQ, fEb);
      fP1 = fEa.Transform(fP1);
      fP2 = fEa.Transform(fP2);
      TVector3 boostv = fQ.BoostVector();
      fP1.Boost(boostv);
      fP2.Boost(boostv);
   }

   // -------------------------
   //  Fix round-off errors
   // -------------------------
   if (fP1.E() <= 0. || fP2.E() <= 0.) {
      fBetaBar = 0;
   } else {
      Double_t ap = fP1.Vect().Mag();
      if (fP1.E() < ap) fP1.SetE(TMath::Sqrt(ap*ap+fM12));
               ap = fP2.Vect().Mag();
      if (fP2.E() < ap) fP2.SetE(TMath::Sqrt(ap*ap+fM22));
   }
}
