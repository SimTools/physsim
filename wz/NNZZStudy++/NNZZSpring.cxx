//*****************************************************************************
//* =====================
//*  NNZZSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> NNZZ generator
//*
//* (Update Record)
//*    2014/09/20  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "NNZZSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(NNZZSpring)
ClassImp(NNZZSpringBuf)
ClassImp(NNZZBases)

//-----------------------------------------------------------------------------
// ==============================
//  class NNZZSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNZZSpring::NNZZSpring(const char      *name,
                       const char      *title,
                             NNZZBases *bases)
           : JSFSpring(name, title, bases)
{
  fEventBuf = new NNZZSpringBuf("NNZZSpringBuf",
                                "NNZZSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new NNZZBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
NNZZSpring::~NNZZSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t NNZZSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    NNZZBases *bs = static_cast<NNZZBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> NNZZBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class NNZZSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t NNZZSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  NNZZBases    *bases   = static_cast<NNZZBases *>(
                          static_cast<NNZZSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 8;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // ne
  pv[1] = bases->fP[1]; // neb
  pv[2] = bases->fP[2]; // f3
  pv[3] = bases->fP[3]; // f4b
  pv[4] = bases->fP[4]; // f5
  pv[5] = bases->fP[5]; // f6b
  pv[6] = pv[2] + pv[3];// Z1
  pv[7] = pv[4] + pv[5];// Z2

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------
  fZBoost       = bases->GetZBoost();
  fEcmIP        = bases->GetEcmIP();
  fXi           = bases->GetXi();
  fEta1         = bases->GetEta1();
  fEta2         = bases->GetEta2();
  fPhi1         = bases->GetPhi1();
  fPhi2         = bases->GetPhi2();
  fQ2ZZ         = bases->GetQ2ZZ();
  fCosZ1        = bases->GetCosZ1 ();
  fPhiZ1        = bases->GetPhiZ1 ();
  fQ2Z1         = bases->GetQ2Z1  ();
  fCosZ1F       = bases->GetCosZ1F();
  fPhiZ1F       = bases->GetPhiZ1F();
  fQ2Z2         = bases->GetQ2Z2  ();
  fCosZ2F       = bases->GetCosZ2F();
  fPhiZ2F       = bases->GetPhiZ2F();
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
  Int_t    idz     = 23;                          // PDG code for Z2
  Int_t    idne    = 12;                          // PDG code for nu_e

  // Z1
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f3 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f4 helicity

  Int_t    islevz1 = color3 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfz1   = 2;                           // color flux id
  Double_t rq2z1   = pv[6].Mag();

  // Z2
  Int_t    idf5    = bases->f5Ptr->GetPID   ();   // PDG code for f5
  Double_t chrg5   = bases->f5Ptr->GetCharge();   // f5 charge
  Double_t m5      = bases->f5Ptr->GetMass  ();   // f5 mass
  Int_t    hel5    = bases->fHelFinal[4];         // f5 helicity
  Double_t color5  = bases->f5Ptr->GetColor();    // color factor for f5

  Int_t    idf6    = bases->f6Ptr->GetPID   ();   // PDG code for f6
  Double_t chrg6   = bases->f6Ptr->GetCharge();   // f6 charge
  Double_t m6      = bases->f6Ptr->GetMass  ();   // f6 mass
  Int_t    hel6    = bases->fHelFinal[5];         // f6 helicity

  Int_t    islevz2 = color5 > 1. ? 201 : 0;  	  // shower level
  Int_t    icfz2   = 3;                           // color flux id
  Double_t rq2z2   = pv[7].Mag();
#ifdef __DEBUG__
  cerr << endl;
  ANL4DVector qcm;
  for (Int_t ip=0; ip<fNparton; ip++) {
    cerr << "pv[" << ip << "] = (" 
         <<  pv[ip].E() << ", "
         <<  pv[ip].Px() << ", "
         <<  pv[ip].Py() << ", "
         <<  pv[ip].Pz() << ") "
	 << "m = " << pv[ip].GetMass() << endl;
    if (ip<fNparton-2) qcm += pv[ip];
  }
  cerr << "qcm = ("
       <<  qcm.E() << ", "
       <<  qcm.Px() << ", "
       <<  qcm.Py() << ", "
       <<  qcm.Pz() << ") " << endl;
  cerr << "----" << endl;
#endif

  //                              No. PID   Mass  Charge   pv    Nd 1st Mom hel  col shower
  new (partons[0]) JSFSpringParton(1, idne,   0.,    0., *qp[0], 0, 0,  0,    0,    0,      0);
  new (partons[1]) JSFSpringParton(2,-idne,   0.,    0., *qp[1], 0, 0,  0,    0,    0,      0);
  new (partons[2]) JSFSpringParton(3,  idz,rq2z1,    0., *qp[6], 2, 5,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton(4,  idz,rq2z2,    0., *qp[7], 2, 7,  0,    0,    0,      0);
  new (partons[4]) JSFSpringParton(5, idf3,   m3, chrg3, *qp[2], 0, 0,  3, hel3,icfz1,islevz1);
  new (partons[5]) JSFSpringParton(6, idf4,   m4, chrg4, *qp[3], 0, 0,  3, hel4,icfz1,islevz1);
  new (partons[6]) JSFSpringParton(7, idf5,   m5, chrg5, *qp[4], 0, 0,  4, hel5,icfz2,islevz2);
  new (partons[7]) JSFSpringParton(8, idf6,   m6, chrg6, *qp[5], 0, 0,  4, hel6,icfz2,islevz2);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class NNZZBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNZZBases::NNZZBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 125.),
           fWidth     ( 0.004),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fZ1ModesLo ( 1),
           fZ1ModesHi (12),
           fZ2ModesLo ( 1),
           fZ2ModesHi (12),
	   fNCALL     (80000),
	   fACC1      (0.05),
	   fACC2      (0.05),
	   fITMX1     (20),
	   fITMX2     (40),
           fZ1BosonPtr( 0),
           fZ2BosonPtr( 0),
           fWBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2ZZ      (0.),
	   fCosZ1     (0.),
	   fPhiZ1     (0.),
	   fQ2Z1      (0.),
	   fCosZ1F    (0.),
	   fPhiZ1F    (0.),
	   fQ2Z2      (0.),
	   fCosZ2F    (0.),
	   fPhiZ2F    (0.),
           fZ1ModePtr (0),
           f3Ptr      (0),
           f4Ptr      (0),
           fZ2ModePtr (0),
           f5Ptr      (0),
           f6Ptr      (0),
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

  cout << "Init nnhbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("NNZZBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.WidthH","0.006")); // M_x [GeV]
  ins >> fWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Ecm","500."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Z1ModesLo","1"));      // Z1 decay mode lo
  ins >> fZ1ModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Z1ModesHi","12"));     // Z1 decay mode hi
  ins >> fZ1ModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Z2ModesLo","1"));      // Z2 decay mode lo
  ins >> fZ2ModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.Z2ModesHi","12"));     // Z2 decay mode hi
  ins >> fZ2ModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZZBases.NCALL","80000"));
  ins >> fNCALL;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fZ1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fZ2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2ZZ         , 0., 1., 1, 1);
  DefineVariable(fXQ2Z1         , 0., 1., 0, 1);
  DefineVariable(fXQ2Z2         , 0., 1., 0, 1);
  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1
  //--
  DefineVariable(fXXi    , 0., 1., 1, 1);
  DefineVariable(fXEta1  , 0., 1., 1, 1);
  DefineVariable(fXEta2  , 0., 1., 1, 1);
  DefineVariable(fXPhi1  , 0., 1., 1, 1);
  DefineVariable(fXPhi21 , 0., 1., 1, 1);
  DefineVariable(fXQ2ZZ  , 0., 1., 1, 1);
  DefineVariable(fXCosZ1 , 0., 1., 0, 1);
  DefineVariable(fXPhiZ1 , 0., 1., 0, 1);
  DefineVariable(fXQ2Z1  , 0., 1., 1, 1);
  DefineVariable(fXCosZ1F, 0., 1., 0, 1);
  DefineVariable(fXPhiZ1F, 0., 1., 0, 1);
  DefineVariable(fXQ2Z2  , 0., 1., 1, 1);
  DefineVariable(fXCosZ2F, 0., 1., 0, 1);
  DefineVariable(fXPhiZ2F, 0., 1., 0, 1);

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
  SetNoOfSample(fNCALL);

  SetTuneValue (1.5);
  SetIteration1(fACC1, fITMX1);
  SetIteration2(fACC2, fITMX2);
}
// --------------------------
//  D-tor
// --------------------------
NNZZBases::~NNZZBases()
{
  delete fZ1BosonPtr;
  delete fZ2BosonPtr;
  delete fWBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t NNZZBases::Func()
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
    eplus   *= fEcmInit/2.;
    eminus  *= fEcmInit/2.;
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

  GENDecayMode *fZ1ModePtr = fZ1BosonPtr->PickMode(fZ1DecayMode, weight, fZ1Mode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  GENDecayMode *fZ2ModePtr = fZ2BosonPtr->PickMode(fZ2DecayMode, weight, fZ2Mode);
  bsWeight *= weight;
  f5Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(0));
  f6Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(1));
  Double_t m5   = f5Ptr->GetMass();
  Double_t m6   = f6Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = 0.; // Mne
  Double_t m2   = 0.; // Mneb
  if (fEcmIP < m1 + m2 + m3 + m4 + m5 + m6) {
    return 0.;
  }

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average for e- is taken
  //  care of here
  SelectHelicities(weight);
  if (weight == 0.) return 0.;
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q^2 of internal lines
  // --------------------------------------------
  Double_t s      = fEcmIP*fEcmIP;
  Double_t rs     = fEcmIP;
  // Z1
  Double_t qz1min = m3 + m4;
  Double_t qz1max = rs - (m1 + m2 + m5 + m6);
#ifndef __ZEROWIDTH__
  fQ2Z1 = fZ1BosonPtr->GetQ2BW(qz1min, qz1max, fXQ2Z1, weight);
#else
  fQ2Z1 = TMath::Power(fZ1BosonPtr->GetMass(),2);
  weight = kPi*fZ1BosonPtr->GetMass()*fZ1BosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // Z2
  Double_t rq2z1  = TMath::Sqrt(fQ2Z1);
  Double_t qz2min = m5 + m6;
  Double_t qz2max = rs - (rq2z1 + m1 + m2);
#ifndef __ZEROWIDTH__
  fQ2Z2 = fZ2BosonPtr->GetQ2BW(qz2min, qz2max, fXQ2Z2, weight);
#else
  fQ2Z2 = TMath::Power(fZ2BosonPtr->GetMass(),2);
  weight = kPi*fZ2BosonPtr->GetMass()*fZ2BosonPtr->GetWidth();
#endif
  bsWeight *= weight;
  Double_t rq2z2  = TMath::Sqrt(fQ2Z2);

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t qzz2mn = TMath::Power(rq2z1 + rq2z2,2);
  Double_t qzz2mx = TMath::Power(rs - (m1+m2),2);
  fQ2ZZ     = qzz2mn + (qzz2mx-qzz2mn)*fXQ2ZZ;
  bsWeight *= qzz2mx - qzz2mn;

  Double_t qzz   = TMath::Sqrt(fQ2ZZ);
  Double_t xilo  = TMath::Log(qzz*(qzz+m1+m2)/s);
  Double_t xihi  = TMath::Log(1.-2.*TMath::Min(m1,m2)/rs);
  fXi       = xilo + (xihi-xilo)*fXXi;
  bsWeight *= xihi - xilo;

  Double_t mw   = fWBosonPtr->GetMass();
  Double_t dm1  = mw*mw/s;
  Double_t dp1  = 1.;
  Double_t dm2  = dp1;
  Double_t dp2  = dm1;

  Double_t etlo = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi =  TMath::Log((1.+dp1)/dm1)/2.;
  fEta1     = etlo + (ethi-etlo)*fXEta1;
  bsWeight *= ethi - etlo;

           etlo = -TMath::Log((1.+dm2)/dp2)/2.;
           ethi =  TMath::Log((1.+dp2)/dm2)/2.;
  fEta2     = etlo + (ethi-etlo)*fXEta2;
  bsWeight *= ethi - etlo;

  fPhi1  = k2Pi * fXPhi1;
  fPhi2  = k2Pi * fXPhi21 + fPhi1;
  fPhi2  = fPhi2 > k2Pi ? fPhi2 - k2Pi : fPhi2;
  bsWeight *= k2Pi*k2Pi;

  fCosZ1 = -1. + 2.*fXCosZ1;
  fPhiZ1 = k2Pi * fXPhiZ1;
  bsWeight *= 2.*k2Pi;

  fCosZ1F = -1. + 2.*fXCosZ1F;
  fPhiZ1F = k2Pi * fXPhiZ1F;
  bsWeight *= 2.*k2Pi;

  fCosZ2F = -1. + 2.*fXCosZ2F;
  fPhiZ2F = k2Pi * fXPhiZ2F;
  bsWeight *= 2.*k2Pi;

  fM[0] = m1;
  fM[1] = m2;

  GENBranch z1branch(fQ2Z1, fCosZ1F, fPhiZ1F, m3*m3    , m4*m4);
  GENBranch z2branch(fQ2Z2, fCosZ2F, fPhiZ2F, m5*m5    , m6*m6);
  GENBranch zzbranch(fQ2ZZ, fCosZ1 , fPhiZ1 , &z1branch, &z2branch);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX(zzbranch);

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill( 1, fEcmIP            , (bsWeight*sigma));
  Xh_fill( 2, fXi               , (bsWeight*sigma));
  Xh_fill( 3, fEta1             , (bsWeight*sigma));
  Xh_fill( 4, fEta2             , (bsWeight*sigma));
  Xh_fill( 5, fPhi1             , (bsWeight*sigma));
  Xh_fill( 6, fPhi2             , (bsWeight*sigma));
  Xh_fill( 7, k2Pi*fXPhi21      , (bsWeight*sigma));
  Xh_fill( 8, TMath::Sqrt(fQ2ZZ), (bsWeight*sigma));
  Xh_fill( 9, fCosZ1            , (bsWeight*sigma));
  Xh_fill(10, fPhiZ1            , (bsWeight*sigma));
  Xh_fill(11, TMath::Sqrt(fQ2Z1), (bsWeight*sigma));
  Xh_fill(12, fCosZ1F           , (bsWeight*sigma));
  Xh_fill(13, fPhiZ1F           , (bsWeight*sigma));
  Xh_fill(14, TMath::Sqrt(fQ2Z2), (bsWeight*sigma));
  Xh_fill(15, fCosZ2F           , (bsWeight*sigma));
  Xh_fill(16, fPhiZ2F           , (bsWeight*sigma));
  Xh_fill(17, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(18, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(19, (Double_t)fZ1Mode , (bsWeight*sigma));
  Xh_fill(20, (Double_t)fZ2Mode , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t NNZZBases::DSigmaDX(GENBranch &zzbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs   = fEcmIP;
  Double_t s    = rs*rs; 
  Double_t eb   = rs/2.;

  Double_t qzz2 = zzbranch.GetQ2();  // q^2_{ZZ}
  Double_t qzz  = TMath::Sqrt(qzz2);
  Double_t q2z1 = zzbranch.GetM12(); // Z1
  Double_t q2z2 = zzbranch.GetM22(); // Z2

  Double_t rxi  = TMath::Exp(fXi);
  weight *= rxi;

  Double_t m1  = fM[0]; // Mne
  Double_t m2  = fM[1]; // Mneb
  Double_t mw  = fWBosonPtr->GetMass();
  Double_t dm1 = mw*mw/s;
  Double_t dp1 = 1.;
  Double_t dm2 = dp1;
  Double_t dp2 = dm1;

  Double_t exp2et1 = TMath::Exp(2.*fEta1);
  Double_t sh1     = TMath::Min(TMath::Sqrt((1.+dp1+dm1)/(1.+exp2et1)   - dm1), 1.);
  Double_t ch1     = TMath::Min(TMath::Sqrt((1.+dp1+dm1)/(1.+1./exp2et1)- dp1), 1.);
  Double_t sn1     = 2.*sh1*ch1;
  Double_t cs1     = (1.+dp1+dm1)*TMath::TanH(fEta1) - dp1 + dm1;
  Double_t fi1     = fPhi1;
  weight *= 4.* (1.+dp1+dm1)/((1.+exp2et1)*(1+1/exp2et1));

  Double_t exp2et2 = TMath::Exp(2.*fEta2);
  Double_t sh2     = TMath::Min(TMath::Sqrt((1.+dp2+dm2)/(1.+exp2et2)   - dm2), 1.);
  Double_t ch2     = TMath::Min(TMath::Sqrt((1.+dp2+dm2)/(1.+1./exp2et2)- dp2), 1.);
  Double_t sn2     = 2.*sh2*ch2;
  Double_t cs2     = (1.+dp2+dm2)*TMath::TanH(fEta2) - dp2 + dm2;
  Double_t fi2     = fPhi2;
  weight *= 4.* (1.+dp2+dm2)/((1.+exp2et2)*(1+1/exp2et2));

  Double_t cs12    = cs1*cs2 + sn1*sn2*TMath::Cos(fi2-fi1);

  Double_t x1, x2;
  if (fEta1 > -fEta2) {
    Double_t rximn = (qzz-m1+m2)*(qzz+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - qzz2/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - qzz2/s);
  } else {
    Double_t rximn = (qzz+m1-m2)*(qzz+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - qzz2/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - qzz2/s);
  }

  Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
  Double_t beta_e = pb/eb;
  fK[0].SetXYZT(0., 0., pb, eb);
  fK[1].SetXYZT(0., 0.,-pb, eb);

  Double_t e1 = eb*x1;
  if (e1 < m1) return 0.;
  Double_t ap1 = TMath::Sqrt((e1-m1)*(e1+m1));
  fP[0].SetXYZT(ap1*sn1*TMath::Cos(fi1), ap1*sn1*TMath::Sin(fi1), ap1*cs1, e1);

  Double_t e2 = eb*x2;
  if (e2 < m2) return 0.;
  Double_t ap2 = TMath::Sqrt((e2-m2)*(e2+m2));
  fP[1].SetXYZT(ap2*sn2*TMath::Cos(fi2), ap2*sn2*TMath::Sin(fi2), ap2*cs2, e2);

  ANL4DVector qcm(rs, 0.,0.,0.);
  ANL4DVector qvzz = qcm - fP[0] - fP[1];
  GENFrame    cmframe;
  GENPhase2   phaseZZ(qvzz, q2z1, q2z2, cmframe, fCosZ1, fPhiZ1, 1);
  ANL4DVector pz1   = phaseZZ.GetFourMomentum(0);
  ANL4DVector pz2   = phaseZZ.GetFourMomentum(1);
  Double_t    betaz = phaseZZ.GetBetaBar();
  if (betaz <= 0.) return 0.;
  weight *= betaz;

  GENBranch &z1branch = *zzbranch.GetBranchPtr(0);
  Double_t cosz1f  = z1branch.GetCosTheta();
  Double_t phiz1f  = z1branch.GetPhi     ();
  Double_t m32     = z1branch.GetM12();
  Double_t m42     = z1branch.GetM22();
  GENFrame zzframe = phaseZZ.GetFrame(1);
  GENPhase2 phaseZ1(pz1, m32, m42, zzframe, cosz1f, phiz1f, 1);
  fP[2] = phaseZ1.GetFourMomentum(0);
  fP[3] = phaseZ1.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betaz1f = phaseZ1.GetBetaBar();
  if (betaz1f <= 0.) return 0.;
  weight *= betaz1f;

  GENBranch &z2branch = *zzbranch.GetBranchPtr(1);
  Double_t cosz2f = z2branch.GetCosTheta();
  Double_t phiz2f = z2branch.GetPhi     ();
  Double_t m52    = z2branch.GetM12();
  Double_t m62    = z2branch.GetM22();
  GENPhase2 phaseZ2(pz2, m52, m62, zzframe, cosz2f, phiz2f, 1);
  fP[4] = phaseZ2.GetFourMomentum(0);
  fP[5] = phaseZ2.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(m52);
  fM[5] = TMath::Sqrt(m62);
  Double_t betaz2f = phaseZ2.GetBetaBar();
  if (betaz2f <= 0.) return 0.;
  weight *= betaz2f;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
#ifdef __DEBUG__
  for (Int_t i=0; i<6; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qz1 = fP[2] + fP[3];
  ANL4DVector qz2 = fP[4] + fP[5];
  cerr << " qz1 = (" 
       << qz1.E () << ","
       << qz1.Px() << ","
       << qz1.Py() << ","
       << qz1.Pz() << ")" << endl;
  cerr << " qz2 = (" 
       << qz2.E () << ","
       << qz2.Px() << ","
       << qz2.Py() << ","
       << qz2.Pz() << ")" << endl;
  ANL4DVector pcm = qz1 + qz2 + fP[0] + fP[1];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
  cerr << " z1mass = " << qz1.GetMass() << endl;
  cerr << " z2mass = " << qz2.GetMass() << endl;
#endif
  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();
  
  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 5;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1./2.;                         // identical particle factor
  Double_t dPhase = kFact * weight;                // phase space factor
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
Double_t NNZZBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor() * f5Ptr->GetColor();

  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t NNZZBases::FullAmplitude()
{
   Double_t mz     = fZ1BosonPtr->GetMass();
   Double_t gamz   = fZ1BosonPtr->GetWidth();

   Double_t qf3    = f3Ptr->GetCharge();
   Double_t t3f3   = f3Ptr->GetISpin();
   Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
   Double_t grzf3  = -kGz*(     - qf3*kSin2W);

   Double_t qf5    = f5Ptr->GetCharge();
   Double_t t3f5   = f5Ptr->GetISpin();
   Double_t glzf5  = -kGz*(t3f5 - qf5*kSin2W);
   Double_t grzf5  = -kGz*(     - qf5*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELFermion ne (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing);
   HELFermion neb(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);

   HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // f3
   HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // f4b
   HELVector  z1(f4b, f3, glzf3, grzf3, mz, gamz);               // Z1

   HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // f5
   HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // f6b
   HELVector  z2(f6b, f5, glzf5, grzf5, mz, gamz);               // Z2

   Complex_t amp = AmpEEtoNNZZ(em, ep, ne, neb, z1, z2);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNZZ()
// --------------------------
Complex_t NNZZBases::AmpEEtoNNZZ(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &ne,
                                 const HELFermion &neb,
                                 const HELVector  &z1,
                                 const HELVector  &z2)
{
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
   Double_t  mz     = fZ1BosonPtr->GetMass();
   Double_t  mw     = fWBosonPtr->GetMass();
   Double_t  gamw   = fWBosonPtr->GetWidth();

   Double_t  glwf  = -kGw*kSqh;
   Double_t  grwf  = 0.;

   Double_t  qe    = -1.;
   Double_t  t3e   = -1./2.;
   Double_t  glze  = -kGz*(t3e - qe*kSin2W);
   Double_t  grze  = -kGz*(    - qe*kSin2W);

   Double_t  t3n   = +1./2.;
   Double_t  glzn  = -kGz*t3n;
   Double_t  grzn  = 0.;

   Double_t  mn    = kMass[0][0][0];
   Double_t  me    = kMass[0][1][0];
   Double_t  gamn  = 0.;
   Double_t  game  = 0.;

   Double_t  mh    = fMass;
   Double_t  gamh  = fWidth;

   Double_t  gwwz  = kGw*kCosW;
   Double_t  gwwh  = kGw*mw;
   Double_t  gzzh  = kGz*mz;

   //---------------
   // Internal lines
   //---------------
   // -- virtual e- with 1 Z
   const HELFermion em01 (em  , z1, glze, grze, me , game );
   const HELFermion em02 (em  , z2, glze, grze, me , game );
   // -- virtual e- with 2 Zs
   const HELFermion em03 (em01, z2, glze, grze, me , game );
   const HELFermion em04 (em02, z1, glze, grze, me , game );
   // -- virtual e+ with 1 Z
   const HELFermion ep05 (ep  , z1, glze, grze, me , game );
   const HELFermion ep06 (ep  , z2, glze, grze, me , game );
   // -- virtual e+ with 2 Zs
   const HELFermion ep07 (ep05, z2, glze, grze, me , game );
   const HELFermion ep08 (ep06, z1, glze, grze, me , game );
   // -- virtual nu  with 1 Z
   const HELFermion n09  (ne  , z1, glzn, grzn, mn , gamn );
   const HELFermion n10  (ne  , z2, glzn, grzn, mn , gamn );
   // -- virtual nu  with 2 Zs
   const HELFermion n11  (n09 , z2, glzn, grzn, mn , gamn );
   const HELFermion n12  (n10 , z1, glzn, grzn, mn , gamn );
   // -- virtual nub with 1 Z
   const HELFermion nb13 (neb , z1, glzn, grzn, mn , gamn );
   const HELFermion nb14 (neb , z2, glzn, grzn, mn , gamn );
   // -- virtual nub with 2 Zs
   const HELFermion nb15 (nb13, z2, glzn, grzn, mn , gamn );
   const HELFermion nb16 (nb14, z1, glzn, grzn, mn , gamn );
   // -- virtual W
   const HELVector  w17  (em  , ne, glwf, grwf, mw , gamw );
   const HELVector  w18  (neb , ep, glwf, grwf, mw , gamw );

   //-----------
   // Non-fusion
   //-----------
   // (t-channel e+e- -> nnb with 2 Zs attached to fermion legs)
   Complex_t amp01 = AmpEEtoNNT (em03, ep  ,  ne, neb );
   Complex_t amp02 = AmpEEtoNNT (em04, ep  ,  ne, neb );
 
   Complex_t amp03 = AmpEEtoNNT (em01, ep06,  ne, neb );
   Complex_t amp04 = AmpEEtoNNT (em02, ep05,  ne, neb );

   Complex_t amp05 = AmpEEtoNNT (em01, ep  , n10, neb );
   Complex_t amp06 = AmpEEtoNNT (em02, ep  , n09, neb );

   Complex_t amp07 = AmpEEtoNNT (em01, ep  , ne , nb14);
   Complex_t amp08 = AmpEEtoNNT (em02, ep  , ne , nb13);

   Complex_t amp09 = AmpEEtoNNT (em  , ep07, ne , neb );
   Complex_t amp10 = AmpEEtoNNT (em  , ep08, ne , neb );

   Complex_t amp11 = AmpEEtoNNT (em  , ep05, n10, neb );
   Complex_t amp12 = AmpEEtoNNT (em  , ep06, n09, neb );

   Complex_t amp13 = AmpEEtoNNT (em  , ep05, ne , nb14);
   Complex_t amp14 = AmpEEtoNNT (em  , ep06, ne , nb13);

   Complex_t amp15 = AmpEEtoNNT (em  , ep  , n11, neb );
   Complex_t amp16 = AmpEEtoNNT (em  , ep  , n12, neb );

   Complex_t amp17 = AmpEEtoNNT (em  , ep  , n09, nb14);
   Complex_t amp18 = AmpEEtoNNT (em  , ep  , n10, nb13);

   Complex_t amp19 = AmpEEtoNNT (em  , ep  , ne , nb15);
   Complex_t amp20 = AmpEEtoNNT (em  , ep  , ne , nb16);

   //-----------
   // Non-fusion
   //-----------
   // (1 Z from a fermion leg and another from internal W)
   HELVector wem1ne(em01, ne, glwf, grwf, mw, gamw);
   Complex_t amp21 = HELVertex (w18, wem1ne, z2, gwwz);
   HELVector wem2ne(em02, ne, glwf, grwf, mw, gamw);
   Complex_t amp22 = HELVertex (w18, wem2ne, z1, gwwz);

   HELVector wnebep5(neb, ep05, glwf, grwf, mw, gamw);
   Complex_t amp23 = HELVertex (wnebep5, w17, z2, gwwz);
   HELVector wnebep6(neb, ep06, glwf, grwf, mw, gamw);
   Complex_t amp24 = HELVertex (wnebep6, w17, z1, gwwz);

   HELVector wemn09(em, n09, glwf, grwf, mw, gamw);
   Complex_t amp25 = HELVertex (w18, wemn09, z2, gwwz);
   HELVector wemn10(em, n10, glwf, grwf, mw, gamw);
   Complex_t amp26 = HELVertex (w18, wemn10, z1, gwwz);

   HELVector wnb13ep(nb13, ep, glwf, grwf, mw, gamw);
   Complex_t amp27 = HELVertex (wnb13ep, w17, z2, gwwz);
   HELVector wnb14ep(nb14, ep, glwf, grwf, mw, gamw);
   Complex_t amp28 = HELVertex (wnb14ep, w17, z1, gwwz);

   //-----------
   // Fusion
   //-----------
   Complex_t amp29 = HELVertex(w18,z1,w17,z2,gwwz,gwwz,mw,gamw,kTRUE);

   //-----------
   // Higgs
   //-----------
   HELScalar hs(z1, z2, gzzh, mh, gamh);
   Complex_t amp30 = HELVertex(w18, w17, hs, gwwh);

   //--------------------------
   // Sum up all the amplitudes
   //--------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 + amp05
                 + amp06 + amp07 + amp08 + amp09 + amp10
                 + amp11 + amp12 + amp13 + amp14 + amp15
                 + amp16 + amp17 + amp18 + amp19 + amp20
		 + amp21 + amp22 + amp23 + amp24 + amp25
                 + amp26 + amp27 + amp28
		 + amp29
		 + amp30;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNT()
// --------------------------
//
Complex_t NNZZBases::AmpEEtoNNT(const HELFermion &em,
                                const HELFermion &ep,
                                const HELFermion &ne,
                                const HELFermion &neb)
{
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
   // t-channel W exchange
   //---------------------------
   Double_t  mw     = fWBosonPtr->GetMass();
   Double_t  gamw   = fWBosonPtr->GetWidth();

   Double_t  glwf  = -kGw*kSqh;
   Double_t  grwf  = 0.;

   HELVector wt(em, ne, glwf, grwf, mw, gamw);
   Complex_t amp = HELVertex(neb, ep, wt, glwf, grwf);
#endif

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void NNZZBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("NNZZBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("NNZZBases.BeamstrahlungFilename","trc500"));
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
  if (!fZ1BosonPtr) fZ1BosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZ1BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZ1BosonPtr->GetMode(m); 
     if (mp && (m<fZ1ModesLo || m>fZ1ModesHi)) {
        mp->Lock();
     }
  }
  fZ1BosonPtr->DebugPrint();

  if (!fZ2BosonPtr) fZ2BosonPtr = new GENPDTZBoson();
  for (Int_t m=1; m<=fZ2BosonPtr->GetEntries(); m++) {
     GENDecayMode *mp = fZ2BosonPtr->GetMode(m); 
     if (mp && (m<fZ2ModesLo || m>fZ2ModesHi)) {
        mp->Lock();
     }
  }
  fZ2BosonPtr->DebugPrint();

  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  fWBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs    = fEcmInit;
  Double_t s     = rs*rs;
  Double_t mw    = fWBosonPtr->GetMass();
  Double_t xilo  = TMath::Log(4.*mw*mw/s);
  Double_t xihi  = 0.;
  Double_t dm1   = mw*mw/s;
  Double_t dp1   = 1.;
  Double_t etlo  = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi  =  TMath::Log((1.+dp1)/dm1)/2.;
  Double_t qzzlo = 0.; //2*mz - 40.;
  Double_t qzzhi = rs;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,  qzzlo,  qzzhi,       50, "qzz"   );
  Xh_init( 9,    -1.,    +1.,       50, "cosz1" );
  Xh_init(10,     0.,   k2Pi,       50, "phiz1" );
  Xh_init(11,    70.,   110.,       50, "Mz1"   );
  Xh_init(12,    -1.,    +1.,       50, "CosF3" );
  Xh_init(13,     0.,   k2Pi,       50, "PhiF3" );
  Xh_init(14,    70.,   110.,       50, "Mz2"   );
  Xh_init(15,    -1.,    +1.,       50, "CosF5" );
  Xh_init(16,     0.,   k2Pi,       50, "PhiF5" );
  Xh_init(17,     0.,     2.,        2, "Helin ");
  Xh_init(18,     0.,     2.,        4, "Helot ");
  Xh_init(19,     0.,    12.,       12, "Z1 mode");
  Xh_init(20,     0.,    12.,       12, "Z2 mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void NNZZBases::Userout()
{
  cout << "End of NNZZBases----------------------------------- "  << endl
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
void NNZZBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][6] = {{-1, +1, -1, +1, -1, +1},
                                           {-1, +1, -1, +1, +1, -1},
                                           {-1, +1, +1, -1, -1, +1},
                                           {-1, +1, +1, -1, +1, -1}};
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
