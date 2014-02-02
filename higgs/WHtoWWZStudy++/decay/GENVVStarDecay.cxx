//*****************************************************************************
//* =====================
//*  GENVVStarDecay
//* =====================
//*  
//* (Description)
//*    Calculate X -> V V(*) -> f1 f2 f3 f4
//*
//* (Update Record)
//*    2014/02/01  K.Fujii	Original version.
//*
//*****************************************************************************

#include "GENVVStarDecay.h"

//#define __DEBUG__
//#define __ZEROWIDTH__
//#define __PHASESPACE__
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

#ifdef __NODECAY__
#define __ZEROWIDTH__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(GENVVStarDecay)

//-----------------------------------------------------------------------------
// ==============================
//  class VVStarDecay
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor & D-tor
// --------------------------
GENVVStarDecay::GENVVStarDecay(GENPDTEntry *xp,
                               GENPDTEntry *v1p,
		               GENPDTEntry *v2p)
         : JSFBases("VVStarDecay", "VVStarDecay"), 
           fXBosonPtr  (xp),
           fV1BosonPtr (v1p),
           fV2BosonPtr (v2p),
           fAmpUtilsPtr(0),
           fQ2V1       (0.),
           fQ2V2       (0.),
           fV1ModePtr  (0),
           fV2ModePtr  (0),
           f1Ptr       (0),
           f2Ptr       (0),
           f3Ptr       (0),
           f4Ptr       (0),
           fCosThetaV1 (0.),
           fPhiV1      (0.),
           fXQ2V1      (0.),
           fCosThetaV1F(0.),
           fPhiV1F     (0.),
           fXQ2V2      (0.),
           fCosThetaV2F(0.),
           fPhiV2F     (0.)
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENVVStarDecay::Initialize()
{
  // --------------------------------------------
  //  Registor random numbers
  // --------------------------------------------
  //--
  //  Final state helicity combination, decay modes, and Q2
  //--
  DefineVariable(fHelCombFinal  , 0., 1., 0, 1);
  DefineVariable(fV1DecayMode   , 0., 1., 0, 1);
  DefineVariable(fV2DecayMode   , 0., 1., 0, 1);
  DefineVariable(fXQ2V1         , 0., 1., 1, 1);
  DefineVariable(fXQ2V2         , 0., 1., 1, 1);
  //--
  //  cos(theta) and phi
  //--
  DefineVariable(fCosThetaV1  , -1.,  +1., 1, 1);
  DefineVariable(fPhiV1       ,  0., k2Pi, 0, 1);
  DefineVariable(fCosThetaV1F , -1.,  +1., 0, 1);
  DefineVariable(fPhiV1F      ,  0., k2Pi, 0, 1);
  DefineVariable(fCosThetaV2F , -1.,  +1., 0, 1);
  DefineVariable(fPhiV2F      ,  0., k2Pi, 0, 1);

  // --------------------------------------------
  //  Set Bases integration parameters
  // --------------------------------------------
  SetNoOfSample(80000);

  SetTuneValue (1.5);
  SetIteration1(0.05, 10);
  SetIteration2(0.05, 20);
}

GENVVStarDecay::~GENVVStarDecay()
{
}

//_____________________________________________________________________________
// --------------------------
//  CalculateWidth
// --------------------------
Double_t GENVVStarDecay::CalculateWidth()
{
  // --------------------------------------------
  //  Calculate width.
  // --------------------------------------------
  Double_t gam = 0.;
  //--
  // X --> V(*)V(*)
  //--
  Bases();
  Userout();
  gam = GetEstimate();
  return gam;
}

//_____________________________________________________________________________
// --------------------------
//  Func
// --------------------------
Double_t GENVVStarDecay::Func()
{
  //  Bases Integrand
  //
  Double_t bsWeight = 1.; // Jacobian factor

  // --------------------------------------------
  //  Select final state
  // --------------------------------------------
  Double_t weight = 1.;
  GENDecayMode *fV1ModePtr = fV1BosonPtr->PickMode(fV1DecayMode, weight, fV1Mode);
  bsWeight *= weight;
  f1Ptr = static_cast<GENPDTEntry *>(fV1ModePtr->At(0));
  f2Ptr = static_cast<GENPDTEntry *>(fV1ModePtr->At(1));
  Double_t m1   = f1Ptr->GetMass();
  Double_t m2   = f2Ptr->GetMass();

  GENDecayMode *fV2ModePtr = fV2BosonPtr->PickMode(fV2DecayMode, weight, fV2Mode);
  bsWeight *= weight;
  f3Ptr = static_cast<GENPDTEntry *>(fV2ModePtr->At(0));
  f4Ptr = static_cast<GENPDTEntry *>(fV2ModePtr->At(1));
  Double_t m3   = f3Ptr->GetMass();
  Double_t m4   = f4Ptr->GetMass();

  // --------------------------------------------
  //  Select helicity combination
  // --------------------------------------------
  //  Notice that spin average is taken
  //  care of here
  fAmpUtilsPtr->SelectHelicities(fHelCombFinal, fHelFinal, weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Decide Q2
  // --------------------------------------------
  Double_t qmin = m1 + m2;
  Double_t qmax = fXBosonPtr->GetMass() - (m3 + m4);
  fQ2V1 = fV1BosonPtr->GetQ2BW(qmin, qmax, fXQ2V1, weight);
  bsWeight *= weight;
  Double_t qw = TMath::Sqrt(fQ2V1);

  qmin = m3 + m4;
  qmax = fXBosonPtr->GetMass() - qw;
  fQ2V2 = fV2BosonPtr->GetQ2BW(qmin, qmax, fXQ2V2, weight);
  bsWeight *= weight;

  // --------------------------------------------
  //  Calculate width
  // --------------------------------------------
  GENBranch v1branch (fQ2V1,  fCosThetaV1F, fPhiV1F, m1*m1, m2*m2);
  GENBranch v2branch (fQ2V2,  fCosThetaV2F, fPhiV2F, m3*m3, m4*m4);
  GENBranch xbranch (TMath::Power(fXBosonPtr->GetMass(),2), fCosThetaV1, fPhiV1, &v1branch, &v2branch);

  Double_t gamma = DGammaDX(xbranch);

  return (bsWeight * gamma);
}

//_____________________________________________________________________________
// --------------------------
//  DGammaDX
// --------------------------
Double_t GENVVStarDecay::DGammaDX(GENBranch &xbranch)
{
  // --------------------------------------------
  //  Phase space
  // --------------------------------------------
  Double_t mass = TMath::Sqrt(xbranch.GetQ2());
  ANL4DVector px(mass, 0., 0., 0.);
  GENFrame    xframe;

  Double_t q2v1  = xbranch.GetM12();
  Double_t q2v2  = xbranch.GetM22();
  Double_t cosv1 = xbranch.GetCosTheta();
  Double_t phiv1 = xbranch.GetPhi     ();
  GENPhase2 phaseX(px, q2v1, q2v2, xframe, cosv1, phiv1, 0);
  ANL4DVector pv1 = phaseX.GetFourMomentum(0);
  ANL4DVector pv2 = phaseX.GetFourMomentum(1);
  Double_t betav1 = phaseX.GetBetaBar();
  if (betav1 <= 0.) return 0.;

  GENBranch &v1branch = *xbranch.GetBranchPtr(0);
  Double_t cosv1f = v1branch.GetCosTheta();
  Double_t phiv1f = v1branch.GetPhi     ();
  Double_t m12    = v1branch.GetM12();
  Double_t m22    = v1branch.GetM22();
  GENPhase2 phaseV1(pv1, m12, m22, xframe, cosv1f, phiv1f, 1);
  fP[0] = phaseV1.GetFourMomentum(0);
  fP[1] = phaseV1.GetFourMomentum(1);
  fM[0] = TMath::Sqrt(m12);
  fM[1] = TMath::Sqrt(m22);
  Double_t betav1f = phaseV1.GetBetaBar();
  if (betav1f <= 0.) return 0.;

  GENBranch &v2branch = *xbranch.GetBranchPtr(1);
  Double_t cosv2f = v2branch.GetCosTheta();
  Double_t phiv2f = v2branch.GetPhi     ();
  Double_t m32    = v2branch.GetM12();
  Double_t m42    = v2branch.GetM22();
  GENPhase2 phaseV2(pv2, m32, m42, xframe, cosv2f, phiv2f, 1);
  fP[2] = phaseV2.GetFourMomentum(0);
  fP[3] = phaseV2.GetFourMomentum(1);
  fM[2] = TMath::Sqrt(m32);
  fM[3] = TMath::Sqrt(m42);
  Double_t betav2f = phaseV2.GetBetaBar();
  if (betav2f <= 0.) return 0.;

  fK.SetXYZT(0., 0., 0., mass);

#if 0
  cerr << " fK = (" 
       << fK.E () << ","
       << fK.Px() << ","
       << fK.Py() << ","
       << fK.Pz() << ")" << endl;
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
  ANL4DVector qw = fP[0] + fP[1];
  cerr << " qw= (" 
       << qw.E () << ","
       << qw.Px() << ","
       << qw.Py() << ","
       << qw.Pz() << ")" << endl;
  ANL4DVector qz = fP[2] + fP[3];
  cerr << " qz= (" 
       << qz.E () << ","
       << qz.Px() << ","
       << qz.Py() << ","
       << qz.Pz() << ")" << endl;
  ANL4DVector qx = qw + qz;
  cerr << " px= (" 
       << qx.E () << ","
       << qx.Px() << ","
       << qx.Py() << ","
       << qx.Pz() << ")" << endl;
#endif

  // -------------------
  //  Amplitude squared
  // -------------------
  Double_t amp2 = AmpSquared();

  // -------------------
  //  Put them together
  // -------------------
  static const Int_t    kNbr  = 3;
  static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));

  Double_t identp = 1.;                                 // identical particle factor
  Double_t dPhase = kFact * betav1 * betav1f * betav2f; // phase space factor
  Double_t flux   = 1./(2.* mass);                      // initial state factor
  Double_t spin   = 1.;                                 // spin factor

  Double_t gamma = identp * flux * spin * amp2 * dPhase; // in [GeV]

  return gamma;
}

//_____________________________________________________________________________
// --------------------------
//  AmpSquared()
// --------------------------
Double_t GENVVStarDecay::AmpSquared()
{
  Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor();
  Int_t     ig1   = f1Ptr->GetGenNo() - 1;
  Int_t     ig2   = f2Ptr->GetGenNo() - 1;
  Double_t  mix   = TMath::Power(kVkm[ig1][ig2],2);
  Complex_t amp   = fAmpUtilsPtr->FullAmplitude(fK,
		                              fP,
		                              fM,
			                      fHelFinal,
		                              fV1BosonPtr,
			            	      fV2BosonPtr,
				              f1Ptr,
				              f2Ptr,
				              f3Ptr,
				              f4Ptr);
  Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;

  return amp2;
}

//_____________________________________________________________________________
// --------------------------
//  Userin
// --------------------------
void GENVVStarDecay::Userin()
{
}

//_____________________________________________________________________________
// --------------------------
//  Userout
// --------------------------
void GENVVStarDecay::Userout()
{
  cout << "GENVVStarDecay Results -------------------------- "    << endl
       << "Width                = " << GetEstimate()  << " +/- "
                                    << GetError()     << " [GeV]" << endl
       << "Number of iterations = " << GetNoOfIterate()           << endl;
}
