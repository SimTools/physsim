#ifndef GENVVSTARDECAY_H
#define GENVVSTARDECAY_H
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

#include "TNamed.h"
#include "JSFBases.h"

#include "HELLib.h"
#include "GENLib.h"

class GENVVStarDecay;

//_______________________________________________________________________
// ==========================
//  class GENVVStarDecay
// ==========================
//-----------------------------------------------------------------------

class GENVVStarDecay: public JSFBases {
public:
  class AmpUtils {
  public:
     virtual Complex_t FullAmplitude(ANL4DVector &k,
                                     ANL4DVector  p[4],
                                     Double_t     m[4],
                                     Int_t        helFinal[4],
                                     GENPDTEntry *v1p,
                                     GENPDTEntry *v2p,
                                     GENPDTEntry *f1p,
                                     GENPDTEntry *f2p,
                                     GENPDTEntry *f3p,
                                     GENPDTEntry *f4p) = 0;

     virtual void   SelectHelicities(Double_t  helComb,
                                     Int_t     helFinal[4],
                                     Double_t &weight) = 0;
  };

public:
  GENVVStarDecay(GENPDTEntry *xp,
                 GENPDTEntry *v1,
                 GENPDTEntry *v2);
  virtual ~GENVVStarDecay();

  virtual Double_t CalculateWidth(); // calculates width
  inline  void     SetAmpUtilsPtr(GENVVStarDecay::AmpUtils *aup) { fAmpUtilsPtr =aup; }

  // ----------------------
  //   Base class methods
  // ----------------------
  virtual void     Userin ();  // Bases user initialization
  virtual void     Userout();  // Bases user output 

  Double_t Func();     // Bases integration function.

  // ----------------------
  //   Utility methods
  // ----------------------
private:
  Double_t  DGammaDX     (GENBranch &hbranch);
  Double_t  AmpSquared   ();

  void      Initialize();
  Double_t  BetaBar(Double_t x1, Double_t x2)
  {
     return TMath::Sqrt(1. - 2.*(x1+x2) + TMath::Power(x1-x2,2));
  }

private:
  // ----------------
  //  Particle Data
  // ----------------
  GENPDTEntry *fXBosonPtr;       //! PD table entry of "V1"
  GENPDTEntry *fV1BosonPtr;      //! PD table entry of "V1"
  GENPDTEntry *fV2BosonPtr;      //! PD table entry of "V2"
  AmpUtils    *fAmpUtilsPtr;     //! pointer to amplitude utility functions

  // ----------------
  //  Event info
  // ----------------
  Double_t       fQ2V1;           // q^2 of final state V1
  Double_t       fQ2V2;           // q^2 of final state V2
  GENDecayMode  *fV1ModePtr;      // pointer to V1 decay mode
  GENDecayMode  *fV2ModePtr;      // pointer to V2 decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st V1 daughter (fermion)
  GENPDTEntry   *f2Ptr;           // point to 2nd V1 daughter (anti-fermion)
  GENPDTEntry   *f3Ptr;           // point to 1st V2 daughter (fermion)
  GENPDTEntry   *f4Ptr;           // point to 2nd V2 daughter (anti-fermion)
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fV1Mode;         // W decay mode
  Int_t          fV2Mode;         // Z decay mode
  ANL4DVector    fK;              // h
  ANL4DVector    fP[4];           // [0,1,2,3] = [fv1, fbv1, fv2, fbv2]
  Double_t       fM[4];           // [0,1,2,3] = [ m1,   m2,  m3,   m4]
  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fV1DecayMode;    // decay mode selector for V1
  Double_t       fV2DecayMode;    // decay mode selector for V2
  Double_t       fCosThetaV1;     // cos(theta_V1) in X   frame
  Double_t       fPhiV1;          // phi_V1        in X   frame
  Double_t       fXQ2V1;          // q^2 of final state V1
  Double_t       fCosThetaV1F;    // cos(theta_f)  in V1  frame
  Double_t       fPhiV1F;         // phi_f         in V1  frame
  Double_t       fXQ2V2;          // q^2 of final state V2
  Double_t       fCosThetaV2F;    // cos(theta_f)  in V2  frame
  Double_t       fPhiV2F;         // phi_f         in V2  frame

  ClassDef(GENVVStarDecay, 1)  // GENVVStarDecay class
};
#endif
