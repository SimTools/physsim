#ifndef __FLAVOURGETTER__
#define __FLAVOURGETTER__
//*************************************************************************
//* =====================
//*  FlavourGetter Class
//* =====================
//*
//* (Description)
//*    A flavour of jets getting class
//* (Requires)
//*     class ANLJetFinder
//*     class ANLTrack
//*     class JSFSIMDST, etc
//* (Provides)
//*     class FlavourGetter
//* (Update Recored)
//*    2000/02/13  K.Ikematsu   Original version
//*    2000/04/28  K.Ikemtasu   Modified GetPrimaryHadronPID method
//*    2001/07/20  K.Ikemtasu   Added GetPrimaryHadronSN method
//*    2001/07/20  K.Ikemtasu   Added GetPartonID method
//*    2001/07/26  K.Ikematsu   Added TTL4JFlavourGetter class
//*
//* $Id$
//*************************************************************************
//
#include "JSFSteer.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "JSFGeneratorParticle.h"
#include "JSFCDCTrack.h"
#include "JSFSpring.h"
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
    fGpid  = 0; fGsn = 0; fGmsn = 0;
  }                                                // default constructor

  virtual ~FlavourGetter() {}                      // default destructor
  Int_t operator()(const ANLJet &jet);

  //private:
  Int_t GetPrimaryHadronPID(const ANLTrack &t);
  Int_t GetPrimaryHadronSN(const ANLTrack &t);
  Int_t GetPartonID(const ANLTrack &t);

private:
  void SearchPrimaryHadron(const ANLTrack &t);

private:
  TClonesArray     *fGen;         //! Pointer to GeneratorParticles Array
  Int_t             fGpid;        //  PID of GeneratorParticle
  Int_t             fGsn;         //  S.N of GeneratorParticle
  Int_t             fGmsn;        //  Mother S.N of GeneratorParticle

  ClassDef(FlavourGetter, 1)      //  FlavourGetter class
};

//_____________________________________________________________________
//  ------------------------
//  TTL4JFlavourGetter Class
//  ------------------------
//
class TTL4JFlavourGetter : public FlavourGetter {
public:
  TTL4JFlavourGetter(JSFSpring &sp) {
    JSFSpringBuf  *spptn   = (JSFSpringBuf *)sp.EventBuf();
    fSpgen = spptn->GetPartons();
  }                                                // default constructor

  Int_t operator()(const ANLJet &jet);

private:
  TClonesArray     *fSpgen;       //! Pointer to SpringParton

  ClassDef(TTL4JFlavourGetter, 1)    //  TTFlavourGetter class
};

#endif
