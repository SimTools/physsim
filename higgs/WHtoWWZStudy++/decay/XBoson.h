#ifndef XBOSON_H
#define XBOSON_H
//*****************************************************************************
//* =====================
//*  XBoson class
//* =====================
//*  
//* (Description)
//*    A PDT entry representing X boson
//*
//* (Update Record)
//*    2014/02/01  K.Fujii	Original version.
//*
//*****************************************************************************

#include "HELLib.h"
#include "GENLib.h"
#include "GENVVStarDecay.h"

// =====================
//  class XBoson
// =====================
//-----------------------------------------------------------------------

class XBoson: public GENPDTEntry {
friend class GENVVStarDecay;
public:
  XBoson(Double_t m      = 150.,
         Double_t fhwz   = 1.,
         Double_t fhwa   = 0.);
  virtual ~XBoson();

  // ----------------------
  //   Utility class
  // ----------------------
private:
  class AmpUtils: public GENVVStarDecay::AmpUtils {
  public:
     // pure virtual functions in GENVVStarDecay to be implemented
     AmpUtils(XBoson *xbp = 0) : fXBosonPtr(xbp) {}
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
     XBoson *fXBosonPtr; // pointer to parent Xboson object
  };

  // ----------------------
  //   Utility methods
  // ----------------------
private:
  Double_t  BetaBar(Double_t x1, Double_t x2)
  {
     return TMath::Sqrt(1. - 2.*(x1+x2) + TMath::Power(x1-x2,2));
  }

private:
  Double_t      fFhwz;
  Double_t      fFhwa;
  // ----------------
  //  Particle Data
  // ----------------
  XBoson::AmpUtils *fAmpUtilsPtr;       //! pointer to XBoson amplitude utility functions
  GENVVStarDecay   *fVVStarDecayPtr;    //! PD table entry of "H+"
  GENPDTWBoson     *fWBosonPtr;         //! PD table entry of "W"
  GENPDTZBoson     *fZBosonPtr;         //! PD table entry of "Z"

  ClassDef(XBoson, 1)  // X boson class
};
#endif
