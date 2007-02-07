//*****************************************************************************
//* =====================
//*  RSAXSpring
//* =====================
//*  
//* (Description)
//*    Kaluza-Klein  ee -> hh generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "RSAXSpring.h"

#include <sstream>

ClassImp(RSAXSpring)
ClassImp(RSAXSpringBuf)
ClassImp(RSAXBases)
ClassImp(GENFrame)
ClassImp(GENPhase2)

//-----------------------------------------------------------------------------
// ==============================
//  class RSAXSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSAXSpring::RSAXSpring(const char      *name,
                       const char      *title,
                             RSAXBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new RSAXSpringBuf("RSAXSpringBuf",
                                "RSAXSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new RSAXBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
RSAXSpring::~RSAXSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class RSAXSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t RSAXSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  RSAXBases    *bases   = (RSAXBases*)((RSAXSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------
  TLorentzVector p2 = bases->fPA;
  TLorentzVector p3 = bases->fP1;
  TLorentzVector p4 = bases->fP2;
  TLorentzVector p1 = p3 + p4;
#if 0
  cerr << " -------- " << endl;
  cerr << " px = (" << p1.E () << ", " << p1.Px() << ", " << p1.Py() << ", " << p1.Pz() << ") " << endl;
  cerr << " pa = (" << p2.E () << ", " << p2.Px() << ", " << p2.Py() << ", " << p2.Pz() << ") " << endl;
  cerr << " p1 = (" << p3.E () << ", " << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << ") " << endl;
  cerr << " p2 = (" << p4.E () << ", " << p4.Px() << ", " << p4.Py() << ", " << p4.Pz() << ") " << endl;
#endif

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  fEcmIP        = bases->GetEcmIP();
  fZBoost       = bases->GetZBoost();
  fCosTheta     = bases->GetCosTheta();
  fPhi          = bases->GetPhi();
  fCosThetaA    = bases->GetCosThetaA();
  fPhiA         = bases->GetPhiA();
  Double_t elab = TMath::Sqrt(fEcmIP*fEcmIP + fZBoost*fZBoost);
  p1.Boost(0.,0.,fZBoost/elab);
  p2.Boost(0.,0.,fZBoost/elab);
  p3.Boost(0.,0.,fZBoost/elab);
  p4.Boost(0.,0.,fZBoost/elab);

  // ----------------------------------------------
  //  Set final state parton infomation
  // ----------------------------------------------

  TVector q1(4);
  TVector q2(4);
  TVector q3(4);
  TVector q4(4);

  fNparton = 4;

  q1(0) = p1.E ();
  q1(1) = p1.Px();
  q1(2) = p1.Py();
  q1(3) = p1.Pz();

  q2(0) = p2.E ();
  q2(1) = p2.Px();
  q2(2) = p2.Py();
  q2(3) = p2.Pz();

  q3(0) = p3.E ();
  q3(1) = p3.Px();
  q3(2) = p3.Py();
  q3(3) = p3.Pz();

  q4(0) = p4.E ();
  q4(1) = p4.Px();
  q4(2) = p4.Py();
  q4(3) = p4.Pz();

#if 0
  Int_t    idx     = 25; 	// PDG code for higgs
#else
  Int_t    idx     = 20000000; 	// PDG code for higgs
#endif
  Int_t    ida     = 22; 	// PDG code for higgs
  Double_t charge  = 0;  	// X charge
  Int_t    islev   = 0;  	// shower level
  Int_t    icf     = 0;  	// color flux id
  Double_t mass    = bases->GetMass();

  new (partons[0]) JSFSpringParton(1,  idx, mass,  charge, q1,
                                   2, 3, 0, 0, icf, islev);
  new (partons[1]) JSFSpringParton(2,  ida, 0., -charge, q2,
                                   0, 0, 0, 0, icf, islev);
  new (partons[2]) JSFSpringParton(3,  ida, 0.,  charge, q3,
                                   0, 0, 1, 0, icf, islev);
  new (partons[3]) JSFSpringParton(4,  ida, 0., -charge, q4,
                                   0, 0, 1, 0, icf, islev);
  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class RSAXBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
RSAXBases::RSAXBases(const char *name, const char *title)
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
           fCosTheta  (0.),
           fPhi       (0.),
           fCosThetaA (0.),
           fPhiA      (0.),
           fZBoost    (0.),
           fEcmIP     (fEcmInit),
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

  cout << "Init rsaxbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("RSAXBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.CosthARange","-1.0 1.0"));
  ins >> fXL[2] >> fXU[2];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.PhiAOverPiRange","0.0 2.0"));
  ins >> fXL[3] >> fXU[3];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Lambda","1000."));    // Lambda [GeV]
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C0","0.")); 	  	 // C_0
  ins >> fC0;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C1","1.")); 		 // C_1
  ins >> fC1;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C2","1.")); 		 // C_2
  ins >> fC2;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.C3","1.")); 		 // C_3
  ins >> fC3;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.MassX","120.")); 	 // M_x [GeV]
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Ecm","1000."));       // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("RSAXBases.Bremstrahlung","1")); // ISR (on)
  ins >> fISR;

  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSAXBases.BeamstrahlungFilepath","/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("RSAXBases.BeamstrahlungFilename","trc500"));
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
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();
  fXL[3] = fXL[3]*TMath::Pi();
  fXU[3] = fXU[3]*TMath::Pi();

  DefineVariable(fCosTheta , fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi      , fXL[1], fXU[1], 0, 1);
  DefineVariable(fCosThetaA, fXL[2], fXU[2], 0, 1);
  DefineVariable(fPhiA     , fXL[3], fXU[3], 0, 1);

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
  SetIteration1(0.05, 100);
  SetIteration2(0.05, 100);

  // --------------------------------------------
  //  Define some plots
  // --------------------------------------------
  Xh_init(1, fXL[0], fXU[0],       50, "Costh" );
  Xh_init(2, fXL[1], fXU[1],       50, "Phi"   );
  Xh_init(3, fXL[2], fXU[2],       50, "CosthA");
  Xh_init(4, fXL[3], fXU[3],       50, "PhiA"  );
  Xh_init(5,     0., fEcmInit*1.1, 50, "Ecm"   );
}
// --------------------------
//  D-tor
// --------------------------
RSAXBases::~RSAXBases()
{
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t RSAXBases::Func()
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


  // --------------------------------------------
  //  Check if reduced Ecm enough for pair prod.
  // --------------------------------------------

  if (fEcmIP < fMass) {
    bsWeight = 0.;
    return 0.;
  }

  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t betax  = (1. - fMass/fEcmIP)*(1. + fMass/fEcmIP);
  Double_t beta_e = 1.;  // ignore electron mass
  Double_t rs     = fEcmIP;
  Double_t s      = rs * rs;

  // -------------------
  //  Amplitude squared
  // -------------------

  Double_t a1     = fC1*kCos2W + fC2*kSin2W;
  Double_t a2     = (fC1 - fC2) * kSinCosW;
  Double_t betax2 = betax * betax;


  Double_t part1  = ((k4Pi*kAlpha)*a1*a1/(4.*fLambda*fLambda)) * s
                    * betax2 * (1. + fCosTheta*fCosTheta);
  Double_t part2  = ((k4Pi*kAlpha/(4.*kSinCosW*kSinCosW)) *a2*a2 * (1. - 4.*kSin2W + 8.*kSin2W*kSin2W) /(8.*fLambda*fLambda)) * s
                    * TMath::Power(s/((rs - kM_z)*(rs + kM_z)), 2)
                    * betax2 * (1. + fCosTheta*fCosTheta);
  Double_t part3  = ((k4Pi*kAlpha/(2.*kSinCosW)) *a1*a2* (1. - 4.*kSin2W) /(4.*fLambda*fLambda)) * s
                    * s/((rs - kM_z)*(rs + kM_z))
                    * betax2 * (1. + fCosTheta*fCosTheta);

  Double_t amp2   = part1 + part2 - part3; // amplitude squared

  Double_t & beta_bar = betax;

  static const Double_t k8Pi = 8. * TMath::Pi();

  Double_t identp = 1.;                           // identical particle factor
  Double_t dPhase = (beta_bar/k8Pi) / k4Pi;       // 2-body phase space factor
  Double_t flux   = 1./(2.*fEcmIP*fEcmIP*beta_e); // beam flux factor
  Double_t spin   = 1./4.;                        // spin average

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  // --------------------------------------------
  //  Phase space
  // --------------------------------------------

  TLorentzVector qcm(0., 0., 0., rs);
  GENFrame       cmframe;

  GENPhase2 phaseCM(qcm, fMass*fMass, 0., cmframe, fCosTheta, fPhi, 0);
  GENFrame       xframe = phaseCM.GetFrame();
  TLorentzVector px     = phaseCM.GetFourMomentum(0);
                 fPA    = phaseCM.GetFourMomentum(1);

  GENPhase2 phaseX(px, 0., 0., cmframe, fCosThetaA, fPhiA, 1);
  fP1 = phaseX.GetFourMomentum(0);
  fP2 = phaseX.GetFourMomentum(1);

  bsWeight *= 1./k4Pi;

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill(1, fCosTheta , (bsWeight*sigma));
  Xh_fill(2, fPhi      , (bsWeight*sigma));
  Xh_fill(3, fCosThetaA, (bsWeight*sigma));
  Xh_fill(4, fPhiA     , (bsWeight*sigma));
  Xh_fill(5, fEcmIP    , (bsWeight*sigma));

  return (bsWeight * sigma) ;
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void RSAXBases::Userout()
{
  cout << "End of RSAXBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Brehmstrahlung       = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}

//-----------------------------------------------------------------------------
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

GENFrame::GENFrame(const TLorentzVector &q, const GENFrame &eb)
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
TLorentzVector GENFrame::Transform(const TLorentzVector &pb)
{
   return TLorentzVector(pb.X()*fEV[0]+pb.Y()*fEV[1]+pb.Z()*fEV[2], pb.E());
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPhase2
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPhase2::GENPhase2(const TLorentzVector &q,
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
TLorentzVector GENPhase2::GetFourMomentum(Int_t i)
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
