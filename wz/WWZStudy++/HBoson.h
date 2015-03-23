#ifndef HBOSON_H
#define HBOSON_H
//*****************************************************************************
//* =====================
//*  HBoson class
//* =====================
//*  
//* (Description)
//*    A PDT entry representing X boson
//*
//* (Update Record)
//*    2014/11/21  K.Fujii	Original version.
//*    2015/03/23  K.Fujii	Improved Gamma_H calculations.
//*
//*****************************************************************************

#include "HELLib.h"
#include "GENLib.h"
#include "GENVVStarDecay.h"

// =====================
//  class HBoson
// =====================
//-----------------------------------------------------------------------

class HBoson: public GENPDTEntry {
friend class GENVVStarDecay;
public:
  HBoson(Double_t m     = 125.,
         Double_t l     = 1000., // Lambda [GeV]
         Double_t a     = 0.,    // VVh
         Double_t b     = 0.,    // FFh
         Double_t bt    = 0.);   // FFth
  virtual ~HBoson();

  // ----------------------
  //   Utility class
  // ----------------------
private:
  class AmpUtils: public GENVVStarDecay::AmpUtils {
  public:
     // pure virtual functions in GENVVStarDecay to be implemented
     AmpUtils(HBoson *xbp = 0) : fHBosonPtr(xbp) {}
     virtual ~AmpUtils() {};
     Complex_t FullAmplitude(ANL4DVector &k,
                             ANL4DVector *p,
                             Double_t    *m,
                             Int_t       *helFinal,
                             GENPDTEntry *v1p,
                             GENPDTEntry *v2p,
                             GENPDTEntry *f1p,
                             GENPDTEntry *f2p,
                             GENPDTEntry *f3p,
                             GENPDTEntry *f4p);

     void   SelectHelicities(Double_t  helComb,
                             Int_t     helFinal[4],
                             Double_t &weight);
  public:
     HBoson   *fHBosonPtr; // pointer to parent Xboson object
  };

  // ----------------------
  //   Utility methods
  // ----------------------
private:
  Double_t  BetaBar(Double_t x1, Double_t x2)
  {
     return TMath::Sqrt(1. - 2.*(x1+x2) + TMath::Power(x1-x2,2));
  }
  Double_t  GammaToFF(Double_t mf, Double_t cf);
  Double_t  GammaToGG();
  Complex_t F (Double_t t);
  Complex_t Fs(Double_t t); // J=0   loop
  Complex_t Ff(Double_t t); // J=1/2 loop
  Complex_t Fv(Double_t t); // J=1   loop
  Double_t  GetLambdaQCD(Double_t alphas, Double_t q);
  Double_t  GetAlphaS   (Double_t lambda, Double_t q);

private:
  Double_t      fLambda;
  Double_t      fA;
  Double_t      fB;
  Double_t      fBtilde;
  // ----------------
  //  Particle Data
  // ----------------
  HBoson::AmpUtils *fWWAmpUtilsPtr;     //! pointer to HBoson amplitude utility functions
  HBoson::AmpUtils *fZZAmpUtilsPtr;     //! pointer to HBoson amplitude utility functions
  GENVVStarDecay   *fWWStarDecayPtr;    //! PD table entry of "h"
  GENVVStarDecay   *fZZStarDecayPtr;    //! PD table entry of "h"
  GENPDTWBoson     *fWMBosonPtr;        //! PD table entry of "W-"
  GENPDTWBoson     *fWPBosonPtr;        //! PD table entry of "W+"
  GENPDTZBoson     *fZ1BosonPtr;        //! PD table entry of "Z1"
  GENPDTZBoson     *fZ2BosonPtr;        //! PD table entry of "Z2"

  ClassDef(HBoson, 1)  // X boson class
};
#endif
