//*****************************************************************************
//* =====================
//*  NNZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- --> NNZH generator
//*
//* (Update Record)
//*    2015/01/08  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "NNZHSpring.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__
//#define __PHASESPACE__

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(NNZHSpring)
ClassImp(NNZHSpringBuf)
ClassImp(NNZHBases)

//-----------------------------------------------------------------------------
// ==============================
//  class NNZHSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNZHSpring::NNZHSpring(const char     *name,
                       const char     *title,
                            NNZHBases *bases)
           : JSFSpring(name, title, bases)
{
  fEventBuf = new NNZHSpringBuf("NNZHSpringBuf",
                               "NNZHSpring event buffer",
                               this);
  if (!bases) { 
    SetBases(new NNZHBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
NNZHSpring::~NNZHSpring()
{
  //delete fEventBuf;
  //delete GetBases();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
Bool_t NNZHSpring::Initialize()
{
  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    NNZHBases *bs = static_cast<NNZHBases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> NNZHBases written to file" << endl;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ==============================
//  class NNZHSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t NNZHSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  NNZHBases     *bases   = static_cast<NNZHBases *>(
                          static_cast<NNZHSpring *>(Module())->GetBases());

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  fNparton = 6;

  static ANL4DVector *pv = 0;
  if (!pv) {
    pv = new ANL4DVector [fNparton];
  }
  pv[0] = bases->fP[0]; // ne
  pv[1] = bases->fP[1]; // neb
  pv[2] = bases->fP[2]; // f3
  pv[3] = bases->fP[3]; // f4b
  pv[4] = bases->fP[4]; // h
  pv[5] = pv[2] + pv[3];// Z

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
  fQ2ZH         = bases->GetQ2ZH();
  fCosZ         = bases->GetCosZ();
  fPhiZ         = bases->GetPhiZ();
  fQ2Z          = bases->GetQ2Z();
  fCosZF        = bases->GetCosZF ();
  fPhiZF        = bases->GetPhiZF ();
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

  Int_t    idne    = 12;                          // PDG code for nu_e

  Int_t    idh     = 25;                          // PDG code for h
  Double_t mass    = bases->GetMass();

  Int_t    idz     = 23;                          // PDG code for Z

  // F
  Int_t    idf3    = bases->f3Ptr->GetPID   ();   // PDG code for f3
  Double_t chrg3   = bases->f3Ptr->GetCharge();   // f3 charge
  Double_t m3      = bases->f3Ptr->GetMass  ();   // f3 mass
  Int_t    hel3    = bases->fHelFinal[2];         // f3 helicity
  Double_t color3  = bases->f3Ptr->GetColor();    // color factor for f3

  Int_t    idf4    = bases->f4Ptr->GetPID   ();   // PDG code for f4
  Double_t chrg4   = bases->f4Ptr->GetCharge();   // f4 charge
  Double_t m4      = bases->f4Ptr->GetMass  ();   // f4 mass
  Int_t    hel4    = bases->fHelFinal[3];         // f4 helicity

  Int_t    islevz  = color3 > 1. ? 101 : 0; 	  // shower level
  Int_t    icfz    = 2;                           // color flux id
  Double_t rq2z    = pv[4].Mag();

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
    qcm += pv[ip];
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
  new (partons[2]) JSFSpringParton(3,  idz, rq2z,    0., *qp[5], 2, 5,  0,    0,    0,      0);
  new (partons[3]) JSFSpringParton(4,  idh, mass,    0., *qp[4], 0, 0,  0,    0,    0,      0);
  new (partons[4]) JSFSpringParton(5, idf3,   m3, chrg3, *qp[2], 0, 0,  4, hel3, icfz, islevz);
  new (partons[5]) JSFSpringParton(6, idf4,   m4, chrg4, *qp[3], 0, 0,  4, hel4, icfz, islevz);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class NNZHBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
NNZHBases::NNZHBases(const char *name, const char *title)
         : JSFBases   (name, title), 
           fMass      ( 125.),
           fWidth     ( 0.004),
           fEcmInit   (1000.),
           fISR       ( 1),
           fBeamStr   ( 1),
           fBeamWidth (0.002),
           fPole      (0.),
           fZModesLo  ( 1),
           fZModesHi  (12),
	   fNCALL     (80000),
	   fACC1      (0.05),
	   fACC2      (0.05),
	   fITMX1     (20),
	   fITMX2     (40),
           fZBosonPtr ( 0),
           fWBosonPtr ( 0),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
	   fXi        (0.),
	   fEta1      (0.),
	   fEta2      (0.),
	   fPhi1      (0.),
	   fPhi2      (0.),
	   fQ2ZH      (0.),
	   fCosZ      (0.),
	   fPhiZ      (0.),
	   fQ2Z       (0.),
	   fCosZF     (0.),
	   fPhiZF     (0.),
           fZModePtr  (0),
           f3Ptr      (0),
           f4Ptr      (0),
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
  stringstream ins(gJSF->Env()->GetValue("NNZHBases.MassH","120.")); // M_x [GeV]
  ins >> fMass;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.WidthH","0.006")); // M_x [GeV]
  ins >> fWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.Ecm","500.")); // E_cm (1TeV)
  ins >> fEcmInit;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.BeamWidth","0.002")); // BmStr (on)
  ins >> fBeamWidth;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.Bremsstrahlung","1"));// ISR (on)
  ins >> fISR;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.Pole","0."));         // electron polarization
  ins >> fPole;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ZModesLo","1"));      // Z decay mode lo
  ins >> fZModesLo;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ZModesHi","12"));     // Z decay mode hi
  ins >> fZModesHi;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ACC1","0.05"));
  ins >> fACC1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ACC2","0.05"));
  ins >> fACC2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ITMX1","20"));
  ins >> fITMX1;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.ITMX2","40"));
  ins >> fITMX2;

  ins.str("");
  ins.clear();
  ins.str(gJSF->Env()->GetValue("NNZHBases.NCALL","80000"));
  ins >> fNCALL;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Initial state helicity combination
  //--
  DefineVariable(fHelCombInitial, 0., 1., 1, 1);
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fZDecayMode    , 0., 1., 0, 1);
  DefineVariable(fXQ2ZH         , 0., 1., 1, 1);
  DefineVariable(fXQ2Z          , 0., 1., 1, 1);
  //--
  //  xi, eta1, eata2, phi1, phi21:=phi2-phi1
  //--
  DefineVariable(fXXi    , 0., 1., 1, 1);
  DefineVariable(fXEta1  , 0., 1., 1, 1);
  DefineVariable(fXEta2  , 0., 1., 1, 1);
  DefineVariable(fXPhi1  , 0., 1., 1, 1);
  DefineVariable(fXPhi21 , 0., 1., 1, 1);
  DefineVariable(fXQ2ZH  , 0., 1., 1, 1);
  DefineVariable(fXCosZ  , 0., 1., 0, 1);
  DefineVariable(fXPhiZ  , 0., 1., 0, 1);
  DefineVariable(fXQ2Z   , 0., 1., 1, 1);
  DefineVariable(fXCosZF , 0., 1., 0, 1);
  DefineVariable(fXPhiZF , 0., 1., 0, 1);

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
NNZHBases::~NNZHBases()
{
  delete fZBosonPtr;
  delete fWBosonPtr;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t NNZHBases::Func()
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

  GENDecayMode *fZModePtr = fZBosonPtr->PickMode(fZDecayMode, weight, fZMode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Check if reduced Ecm enough for prod.
  // --------------------------------------------
  Double_t m1   = 0.; // Mne
  Double_t m2   = 0.; // Mneb
  Double_t m5   = fMass; // m_h
  if (fEcmIP < m1 + m2 + m3 + m4 + m5) {
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
  // Z
  Double_t qzmin = m3 + m4;
  Double_t qzmax = rs - (m1 + m2 + m5);
#ifndef __ZEROWIDTH__
  fQ2Z = fZBosonPtr->GetQ2BW(qzmin, qzmax, fXQ2Z, weight);
#else
  fQ2Z = TMath::Power(fZBosonPtr->GetMass(),2);
  weight = kPi*fZBosonPtr->GetMass()*fZBosonPtr->GetWidth();
#endif
  bsWeight *= weight;

  // --------------------------------------------
  //  Handle kinematics here
  // --------------------------------------------
  Double_t rq2z   = TMath::Sqrt(fQ2Z);
  Double_t qzh2mn = TMath::Power(rq2z + m5,2);
  Double_t qzh2mx = TMath::Power(rs - (m1+m2),2);
  fQ2ZH     = qzh2mn + (qzh2mx-qzh2mn)*fXQ2ZH;
  bsWeight *= qzh2mx - qzh2mn;

  Double_t qzh   = TMath::Sqrt(fQ2ZH);
  Double_t xilo  = TMath::Log(qzh*(qzh+m1+m2)/s);
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

  fCosZ  = -1. + 2.*fXCosZ;
  fPhiZ  = k2Pi * fXPhiZ;
  bsWeight *= 2.*k2Pi;

  fCosZF = -1. + 2.*fXCosZF;
  fPhiZF = k2Pi * fXPhiZF;
  bsWeight *= 2.*k2Pi;

  fM[0] = m1;
  fM[1] = m2;

  GENBranch zbranch(fQ2Z  , fCosZF, fPhiZF, m3*m3   , m4*m4);
  GENBranch zhbranch(fQ2ZH, fCosZ , fPhiZ , &zbranch, m5*m5);

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t sigma = DSigmaDX(zhbranch);

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
  Xh_fill( 8, TMath::Sqrt(fQ2ZH), (bsWeight*sigma));
  Xh_fill( 9, fCosZ             , (bsWeight*sigma));
  Xh_fill(10, fPhiZ             , (bsWeight*sigma));
  Xh_fill(11, TMath::Sqrt(fQ2Z) , (bsWeight*sigma));
  Xh_fill(12, fCosZF            , (bsWeight*sigma));
  Xh_fill(13, fPhiZF            , (bsWeight*sigma));
  Xh_fill(14, (Double_t)fJCombI , (bsWeight*sigma));
  Xh_fill(15, (Double_t)fJCombF , (bsWeight*sigma));
  Xh_fill(16, (Double_t)fZMode  , (bsWeight*sigma));

  return (bsWeight * sigma);
}

//_____________________________________________________________________________
// --------------------------
//  DSigmaDX
// --------------------------
Double_t NNZHBases::DSigmaDX(GENBranch &zhbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t weight = 1.;

  Double_t rs   = fEcmIP;
  Double_t s    = rs*rs; 
  Double_t eb   = rs/2.;

  Double_t q2zh = zhbranch.GetQ2();  // q^2_{ZH}
  Double_t qzh  = TMath::Sqrt(q2zh);
  Double_t q2z  = zhbranch.GetM12(); // Z
  Double_t q2h  = zhbranch.GetM22(); // h

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
    Double_t rximn = (qzh-m1+m2)*(qzh+m1+m2)/s;
    Double_t rximx = 1. - 2*m1/rs;
    x1    = 1. - rxi;
    x2    = (rxi - q2zh/s)/(1-x1*(1.-cs12)/2.);
    if (1.-x2 < rximn || 1.-x2 > rximx) return 0.;
    weight *= s*x1*x2*x2/(rxi - q2zh/s);
  } else {
    Double_t rximn = (qzh+m1-m2)*(qzh+m1+m2)/s;
    Double_t rximx = 1. - 2*m2/rs;
    x2    = 1. - rxi;
    x1    = (rxi - q2zh/s)/(1-x2*(1.-cs12)/2.);
    if (1.-x1 < rximn || 1.-x1 > rximx) return 0.;
    weight *= s*x2*x1*x1/(rxi - q2zh/s);
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
  ANL4DVector pzh = qcm - fP[0] - fP[1];
  GENFrame    cmframe;
  GENPhase2   phaseZH(pzh, q2z, q2h, cmframe, fCosZ, fPhiZ, 1);
  ANL4DVector pz    = phaseZH.GetFourMomentum(0);
  fP[4] = phaseZH.GetFourMomentum(1);
  fM[4] = TMath::Sqrt(q2h);
  Double_t    betaz = phaseZH.GetBetaBar();
  if (betaz <= 0.) return 0.;
  weight *= betaz;

  GENBranch &zbranch = *zhbranch.GetBranchPtr(0);
  Double_t coszf  = zbranch.GetCosTheta();
  Double_t phizf  = zbranch.GetPhi     ();
  Double_t m32    = zbranch.GetM12();
  Double_t m42    = zbranch.GetM22();
  //GENFrame zhframe = phaseZH.GetFrame(1);
  //GENPhase2 phaseZ(pz, m32, m42, zhframe, coszf, phizf, 1);
  GENPhase2 phaseZ(pz, m32, m42, cmframe, coszf, phizf, 1);
  fP[2] = phaseZ.GetFourMomentum(0);
  fP[3] = phaseZ.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betazf = phaseZ.GetBetaBar();
  if (betazf <= 0.) return 0.;
  weight *= betazf;

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------
#ifdef __DEBUG__
  for (Int_t i=0; i<2; i++) {
    cerr << " fK[" << i << "] = (" 
         << fK[i].E () << ","
         << fK[i].Px() << ","
         << fK[i].Py() << ","
         << fK[i].Pz() << ")" << endl;
  }
  for (Int_t i=0; i<4; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qvz = fP[2] + fP[3];
  cerr << " pz  = (" 
       << qvz.E () << ","
       << qvz.Px() << ","
       << qvz.Py() << ","
       << qvz.Pz() << ")" << endl;
  cerr << " ph  = (" 
       << fP[4].E () << ","
       << fP[4].Px() << ","
       << fP[4].Py() << ","
       << fP[4].Pz() << ")" << endl;
  ANL4DVector pcm = qvz + fP[0] + fP[1] + fP[4];
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
  cerr << " zmass = " << qvz.GetMass() << endl;
  cerr << " hmass = " << fP[4].GetMass() << endl;
#endif
  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();
  
  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 4;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                            // identical particle factor
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
Double_t NNZHBases::AmpSquared()
{
  Double_t  color = f3Ptr->GetColor();

  Complex_t amp   = FullAmplitude();
  Double_t  amp2  = TMath::Power(abs(amp),2) * color;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t NNZHBases::FullAmplitude()
{
   Double_t mz     = fZBosonPtr->GetMass();
   Double_t gamz   = fZBosonPtr->GetWidth();

   Double_t qf3    = f3Ptr->GetCharge();
   Double_t t3f3   = f3Ptr->GetISpin();
   Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
   Double_t grzf3  = -kGz*(     - qf3*kSin2W);

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
   HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);

   HELFermion ne (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing);
   HELFermion neb(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);

   HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // f3
   HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // f4b
   HELVector  zf(f4b, f3, glzf3, grzf3, mz, gamz);               // Z1

   HELScalar  hs(fP[4]);

   Complex_t amp = AmpEEtoNNZH(em, ep, ne, neb, zf, hs);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNZH()
// --------------------------
Complex_t NNZHBases::AmpEEtoNNZH(const HELFermion &em,
                                 const HELFermion &ep,
                                 const HELFermion &ne,
                                 const HELFermion &neb,
                                 const HELVector  &zf,
                                 const HELScalar  &hs)
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
   Double_t  mw     = fWBosonPtr->GetMass();
   Double_t  gamw   = fWBosonPtr->GetWidth();
   Double_t  mz     = fZBosonPtr->GetMass();
   Double_t  gamz   = fZBosonPtr->GetWidth();

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

   Double_t  gwwz  = kGw*kCosW;
   Double_t  gwwh  = kGw*mw;
   Double_t  gzzh  = kGz*mz;

   //------------------------------------------
   // Amplitude (1): nnZ* followed by Z* -> Zh
   //------------------------------------------
   const HELVector  zvzfhs(zf  , hs, gzzh, mz, gamz);
   Complex_t amp01 = AmpEEtoNNZ(em, ep, ne, neb, zvzfhs);

   //----------------------------------------------
   // Amplitude (2-5): e+(*) + e-(*) -> n(*)n(*)H
   //----------------------------------------------
   // (2)
   const HELFermion emvemzf(em, zf, glze, grze, me , game);
   Complex_t amp02 = AmpEEtoNNH(emvemzf, ep, ne, neb, hs);

   // (3)
   const HELFermion epvepzf(ep, zf, glze, grze, me , game);
   Complex_t amp03 = AmpEEtoNNH(em, epvepzf, ne, neb, hs);
   
   // (4)
   const HELFermion nevnezf(ne, zf, glzn, grzn, mn , gamn);
   Complex_t amp04 = AmpEEtoNNH(em, ep, nevnezf, neb, hs);

   // (5)
   const HELFermion nebvnebzf(neb, zf, glzn, grzn, mn , gamn);
   Complex_t amp05 = AmpEEtoNNH(em, ep, ne, nebvnebzf, hs);

   //-----------------
   // Amplitude (6-7)
   //-----------------
   HELVector wmvemne (em , ne, glwf, grwf, mw, gamw);
   HELVector wpvepneb(neb, ep, glwf, grwf, mw, gamw);

   HELVector wvwmvzf(wmvemne, zf, gwwz, mw, gamw);
   Complex_t amp06 = HELVertex(wpvepneb, wvwmvzf, hs, gwwh);

   HELVector wvzfwpv(zf, wpvepneb, gwwz, mw, gamw);
   Complex_t amp07 = HELVertex(wvzfwpv, wmvemne, hs, gwwh);

   //--------------------------
   // Sum up all the amplitudes
   //--------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 + amp05 + amp06 + amp07;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNZ()
// --------------------------
Complex_t NNZHBases::AmpEEtoNNZ(const HELFermion &em,
                                const HELFermion &ep,
                                const HELFermion &ne,
                                const HELFermion &neb,
                                const HELVector  &zf)
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

   Double_t  gwwz  = kGw*kCosW;

   //---------------
   // Internal lines
   //---------------
   // (1)
   const HELVector  w1     (em  , ne, glwf, grwf, mw , gamw );
   const HELFermion epnebw1(neb , w1, glwf, grwf, me , game );
   Complex_t amp01 = HELVertex(epnebw1, ep, zf, glze, grze);

   // (2)
   const HELVector  w2     (neb , ep, glwf, grwf, mw , gamw );
   const HELFermion ememzf (em  , zf, glze, grze, me , game );
   Complex_t amp02 = HELVertex(ememzf, ne, w2, glwf, grwf);
   
   // (3)
   const HELFermion nenezf (ne  , zf, glzn, grzn, mn , gamn );
   Complex_t amp03 = HELVertex(em, nenezf, w2, glwf, grwf);

   // (4)
   const HELFermion nbnbzf (neb , zf, glzn, grzn, mn , gamn );
   Complex_t amp04 = HELVertex(nbnbzf, ep, w1, glwf, grwf);

   //-----------
   // Fusion
   //-----------
   // (5)
   Complex_t amp05 = HELVertex(w2, w1, zf, gwwz);

   //--------------------------
   // Sum up all the amplitudes
   //--------------------------
   Complex_t amp = amp01 + amp02 + amp03 + amp04 + amp05;
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  AmpEEtoNNH()
// --------------------------
Complex_t NNZHBases::AmpEEtoNNH(const HELFermion &em,
                                const HELFermion &ep,
                                const HELFermion &ne,
                                const HELFermion &neb,
                                const HELScalar  &hs)
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
   Double_t mw     = fWBosonPtr->GetMass();
   Double_t gamw   = fWBosonPtr->GetWidth();
   Double_t glwf   = -kGw*kSqh;
   Double_t grwf   = 0.;

   HELVector w1(em , ne, glwf, grwf, mw, gamw);
   HELVector w2(neb, ep, glwf, grwf, mw, gamw);

   Double_t gwwh   = kGw*mw;
   Complex_t amp = HELVertex(w1, w2, hs, gwwh);
#endif /* end __PHASESPACE__ */

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void NNZHBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("NNZHBases.BeamstrahlungFilepath",
                                  "/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("NNZHBases.BeamstrahlungFilename","trc500"));
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

  if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
  fWBosonPtr->DebugPrint();

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Double_t rs    = fEcmInit;
  Double_t s     = rs*rs;
  Double_t mw    = fWBosonPtr->GetMass();
  Double_t mz    = fZBosonPtr->GetMass();
  Double_t xilo  = TMath::Log(mz*mz/s);
  Double_t xihi  = 0.;
  Double_t dm1   = mw*mw/s;
  Double_t dp1   = 1.;
  Double_t etlo  = -TMath::Log((1.+dm1)/dp1)/2.;
  Double_t ethi  =  TMath::Log((1.+dp1)/dm1)/2.;

  Xh_init( 1,     0., fEcmInit*1.1, 50, "Ecm"   );
  Xh_init( 2,   xilo,   xihi,       50, "xi"    );
  Xh_init( 3,   etlo,   ethi,       50, "eta1"  );
  Xh_init( 4,  -ethi,  -etlo,       50, "eta2"  );
  Xh_init( 5,     0.,   k2Pi,       50, "phi1"  );
  Xh_init( 6,     0.,   k2Pi,       50, "phi2"  );
  Xh_init( 7,     0.,   k2Pi,       50, "phi21" );
  Xh_init( 8,   190.,  1190.,       50, "mzh"   );
  Xh_init( 9,    -1.,    +1.,       50, "cosz"  );
  Xh_init(10,     0.,   k2Pi,       50, "phiz"  );
  Xh_init(11,    70.,   110.,       50, "mz"    );
  Xh_init(12,    -1.,    +1.,       50, "coszf" );
  Xh_init(13,     0.,   k2Pi,       50, "phizf" );
  Xh_init(14,     0.,     2.,        2, "Helin ");
  Xh_init(15,     0.,     4.,        4, "Helot ");
  Xh_init(16,     0.,    12.,       12, "Z mode");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void NNZHBases::Userout()
{
  cout << "End of NNZHBases----------------------------------- "  << endl
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
void NNZHBases::SelectHelicities(Double_t &weight)
{
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][5] = {{-1, +1, -1, +1, 0},
                                           {-1, +1, +1, -1, 0}};
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
   weight = kNf;
}
