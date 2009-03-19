//*****************************************************************************
//* =====================
//*  KKhhSpring
//* =====================
//*  
//* (Description)
//*    Kaluza-Klein  ee -> hh generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2003/XX/XX  Nicolas Delerue	Original version.
//*    2005/03/01  K.Fujii		Cleaned up and fixed some bugs with 
//*                                     higgs momenta and Lorentz boost.
//*
//*****************************************************************************

#include "JSFSteer.h"
#include "KKhhSpring.h"

#include <sstream>

ClassImp(KKhhSpring)
ClassImp(KKhhSpringBuf)
ClassImp(KKhhBases)


//-----------------------------------------------------------------------------
// ==============================
//  class KKhhSpring
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
KKhhSpring::KKhhSpring(const char      *name,
                       const char      *title,
                             KKhhBases *bases)
          : JSFSpring(name, title, bases)
{
  fEventBuf = new KKhhSpringBuf("KKhhSpringBuf",
                                "KKhhSpring event buffer",
                                this);
  if (!bases) { 
    SetBases(new KKhhBases);
  }
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
KKhhSpring::~KKhhSpring()
{
  delete fEventBuf;
  delete GetBases();
}


//-----------------------------------------------------------------------------
// ==============================
//  class KKhhSpringBuf
// ==============================
//_____________________________________________________________________________
// --------------------------
//  SetPartons
// --------------------------
Bool_t KKhhSpringBuf::SetPartons()
{
  TClonesArray &partons = *fPartons;
  KKhhBases    *bases   = (KKhhBases*)((KKhhSpring*)Module())->GetBases();

  // ----------------------------------------------
  //  Set 4-momenta of final state particles in CM
  // ----------------------------------------------

  fNparton = 2;
  Double_t fZBoost = bases->GetZboost();
  Double_t ecm     = bases->GetEcmIP();
  Double_t eh      = ecm/2.;
  fCosth           = bases->GetCosth();
  fPhi             = bases->GetPhi();
  Double_t fMass   = bases->GetMass();

  Double_t sinth     = TMath::Sqrt((1. - fCosth)*(1. + fCosth)); 
  Double_t ph        = TMath::Sqrt((eh - fMass)*(eh + fMass));

  TLorentzVector lp1 = TLorentzVector(ph*sinth*TMath::Cos(fPhi), 
                                      ph*sinth*TMath::Sin(fPhi), 
                                      ph*fCosth, 
                                      eh);
  TLorentzVector lp2 = TLorentzVector(-lp1.Px(),
                                      -lp1.Py(),
                                      -lp1.Pz(),
                                      eh);

  // ----------------------------------------------
  //  Boost them to lab. frame
  // ----------------------------------------------

  Double_t elab = TMath::Sqrt(ecm*ecm + fZBoost*fZBoost);
  lp1.Boost(0.,0.,fZBoost/elab);
  lp2.Boost(0.,0.,fZBoost/elab);

  // ----------------------------------------------
  //  Set final state parton infomation
  // ----------------------------------------------
  //     we generate only higgs (0 charge)  
  //

  TVector p1(4);
  TVector p2(4);

  p1(0) = lp1.E();
  p1(1) = lp1.Px();
  p1(2) = lp1.Py();
  p1(3) = lp1.Pz();

  p2(0) = lp2.E();
  p2(1) = lp2.Px();
  p2(2) = lp2.Py();
  p2(3) = lp2.Pz();

  Int_t    id     = 25; // PDG code for higgs
  Double_t charge = 0;  // higgs charge
  Int_t    islev  = 0;  // shower level
  Int_t    icf    = 0;  // color flux id

  if (id < 10) { islev = 101; icf = 1; }

  new (partons[0]) JSFSpringParton(1,  id, fMass,  charge, p1,
                                   0, 0, 0, 0, icf, islev);
  new (partons[1]) JSFSpringParton(2,  id, fMass, -charge, p2,
                                   0, 0, 0, 0, icf, islev);

  return kTRUE ;
}

//-----------------------------------------------------------------------------
// ==============================
//  class KKhhBases
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
KKhhBases::KKhhBases(const char *name, const char *title)
         : JSFBases(name, title), 
           fLambda (   1.),
           fMScale (2000.),
           fMass   ( 120.),
           fEcmInit(1000.),
           fISR    ( 1),
           fBeamStr( 1),
           fBeamWidth(0.002),
           fCosth  (0.),
           fPhi    (0.),
           fZBoost (0.),
           fEcmIP  (fEcmInit),
           fR_BW_m (0),
           fR_BW_p (0),
           fR_BS_m (0),
           fR_BS_p (0),
           fR_ISR_var (0),
           fR_ISR_side(0)
{
  //  Constructor of bases.  Default parameter should be initialized here
  //
  // --------------------------------------------
  //  Get parameters from jsf.conf, if specified
  // --------------------------------------------

  cout << "Init kkbases " << endl;
  
  using namespace std;
  stringstream ins(gJSF->Env()->GetValue("KKhhBases.CosthRange","-1.0 1.0"));
  ins >> fXL[0] >> fXU[0];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.PhiOverPiRange","0.0 2.0"));
  ins >> fXL[1] >> fXU[1];

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.lambda","+1"));      // Lambda (+1)
  ins >> fLambda;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.MassScale","2000.")); // M_s (2TeV)
  ins >> fMScale;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.MassHiggs","120."));  // M_h (120GeV)
  ins >> fMass;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.Ecm","1000."));      // E_cm (1TeV)
  ins >> fEcmInit;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.beamstrahlung","1")); // BmStr (on)
  ins >> fBeamStr;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.BeamWidth","0.002")); // Beam energy spread
  ins >> fBeamWidth;

  ins.clear();
  ins.str(gJSF->Env()->GetValue("KKhhBases.bremstrahlung","1")); // ISR (on)
  ins >> fISR;

  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  cos(theta) and phi
  //--
  fXL[1] = fXL[1]*TMath::Pi();
  fXU[1] = fXU[1]*TMath::Pi();

  DefineVariable(fCosth, fXL[0], fXU[0], 1, 1);
  DefineVariable(fPhi,   fXL[1], fXU[1], 0, 1);

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
}

// --------------------------
//  D-tor
// --------------------------
KKhhBases::~KKhhBases()
{
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t KKhhBases::Func()
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
    eplus   *= 0.5*fEcmInit;
    eminus  *= 0.5*fEcmInit;
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

  if (fEcmIP < 2.*fMass) {
    bsWeight = 0.;
    return 0.;
  }

  fZBoost = eminus - eplus; // P_z of the cm system after ISR and beamstrahlung

  // --------------------------------------------
  //  Calcuate differential cross section
  // --------------------------------------------

  Double_t eh     = fEcmIP/2.;

  Double_t beta_h = TMath::Sqrt((eh - fMass)*(eh + fMass))/eh;
  Double_t beta_e = 1.;  // ignore electron mass


  // -------------------
  //  Amplitude squared
  // -------------------
  static const Double_t k4Pi = 4. * TMath::Pi();

  Double_t part1  = 1./2.;
  Double_t part2  = (k4Pi * fLambda) * TMath::Power(fEcmIP/fMScale, 4);
           part2 *= part2;

  Double_t part3  = (TMath::Power(beta_h, 4)/4.)
                    * (fCosth*fCosth) * (1. - fCosth) * (1. + fCosth);
           part3  = part3 > 0. ? part3 : 0.;

  Double_t amp2   = part1 * part2 * part3; // amplitude squared
#if 0
  Double_t X_1 = TMath::Power(fMass/fEcmIP,2);
  Double_t X_2 = TMath::Power(fMass/fEcmIP,2);
  Double_t sqrt_beta_bar = 1 - 2*(X_1 + X_2) + TMath::Power((X_1 - X_2), 2);
  Double_t beta_bar      = TMath::Sqrt(sqrt_beta_bar);   
#else
  Double_t & beta_bar = beta_h;
#endif

  static const Double_t k8Pi = 8. * TMath::Pi();

  Double_t identp = 1./2.;                        // identical particle factor
  Double_t dPhase = (beta_bar/k8Pi) / k4Pi;       // 2-body phase space factor
  Double_t flux   = 1./(2.*fEcmIP*fEcmIP*beta_e); // beam flux factor
  Double_t spin   = 1./4.;                        // spin average

  Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
           sigma *= kGeV2fb;                              // now in [fb]

  // --------------------------------------------
  //  Fill plots
  // --------------------------------------------

  Xh_fill(1, fCosth, (bsWeight*sigma));
  Xh_fill(2, fPhi,   (bsWeight*sigma));
  Xh_fill(3, fEcmIP, (bsWeight*sigma));

  return (bsWeight * sigma) ;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void KKhhBases::Userin()
{
  // --------------------------------------------
  //  Open beamstrahlung data
  // --------------------------------------------
  string bsfiledir;
  string bsfilename;

  if (fBeamStr) {
    stringstream ins(gJSF->Env()->GetValue("KKhhBases.beamstrahlungFilepath","/proj/soft/jsf/pro/data/bsdata/"));
    ins >> bsfiledir;

    ins.clear();
    ins.str(gJSF->Env()->GetValue("KKhhBases.beamstrahlungFilename","trc500"));
    ins >> bsfilename;
    
    TFile *fBeamFile = new TFile((bsfiledir+bsfilename+".root").data());
    if (!fBeamFile) {
      cerr << " Unable to open a file for beam strahlung" << endl;
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
  //  Define some plots
  // --------------------------------------------
  Xh_init(1, fXL[0], fXU[0],       50, "Costh");
  Xh_init(2, fXL[1], fXU[1],       50, "Phi");
  Xh_init(3, 0.,     fEcmInit*1.1, 50, "Ecm");
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void KKhhBases::Userout()
{
  cout << "End of KKhhBases----------------------------------- "  << endl
       << "Ecm                  = " << fEcmInit << " [GeV]   "    << endl
       << "Beamstrahlung        = " << (fBeamStr ? "on" : "off")  << endl
       << "Brehmstrahlung       = " << (fISR     ? "on" : "off")  << endl
       << "Total Cross section  = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [fb]"  << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}
