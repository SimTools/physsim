#ifndef __FLAVOURGETTER__
#define __FLAVOURGETTER__
//*************************************************************************
//* =====================
//*  FlavourGetter Class
//* =====================
//*
//* (Description)
//*    A very primitive flavour of jets getting class.
//* (Requires)
//*     class ANLJetFinder
//*     class ANLTrack
//*     class JSFSIMDST, etc.
//* (Provides)
//*     class FlavourGetter
//* (Update Recored)
//*    2000/02/13  K.Ikematsu   Original version.
//*    2000/04/28  K.Ikemtasu   Modified GetPrimaryHadronPID method.
//*
//*************************************************************************
//
#include "JSFSteer.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "JSFGeneratorParticle.h"
#include "JSFCDCTrack.h"
#include "Anlib.h"
//_____________________________________________________________________
//  -------------------
//  FlavourGetter Class
//  -------------------
//
class FlavourGetter {
public:
  FlavourGetter() {
    JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
    JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
    fGen   = evt->GetGeneratorParticles();
  }                                                // default constructor

  virtual ~FlavourGetter() {}                      // default destructor
  Int_t operator()(const ANLJet &jet);
  Int_t GetPrimaryHadronPID(const ANLTrack &t);

private:
  TClonesArray     *fGen;         //! Pointer to GeneratorParticle

  ClassDef(FlavourGetter, 1)      //  FlavourGetter Example
};

#endif
