#ifndef __ANLVTXTAGGER__
#define __ANLVTXTAGGER__

//*************************************************************************
//* ====---=============
//*  ANLVTXTagger Class
//* ====================
//*
//* (Description)
//*    A very primitive vertex tagging class.
//* (Requires)
//*	class ANLJetFinder
//*	class ANLTrack
//*	class JSFSIMDST, etc.
//* (Provides)
//*     class ANLVTXTagger
//* (Update Recored)
//*    1999/10/09  K.Fujii	Original version.
//*    1999/10/18  K.Ikematsu   Impliment Ks, Lambda, Sigma removal.
//*    2000/03/18  K.Ikematsu   Added SetNsig and SetNoff method.
//*    2000/03/18  K.Ikematsu   Added GetNsig and GetNoff method.
//*
//*************************************************************************
//
#include "JSFSteer.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "JSFGeneratorParticle.h"
#include "JSFCDCTrack.h"
#include "ANLJetFinder.h"
#include "ANLTrack.h"
//_____________________________________________________________________
//  -------------------
//  ANLVTXTagger Class
//  -------------------
//
class ANLVTXTagger {
public:
  ANLVTXTagger(Double_t nsig=0., Int_t noff=1)
            : fNsigCut(nsig), fNoffVTracks(noff) {
    JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
    JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
    fParam = sds->Param();
    fGen   = evt->GetGeneratorParticles();
  }                                                // default constructor

  virtual ~ANLVTXTagger(){}                        // default destructor

  Bool_t   operator()(const ANLJet &jet);
  Double_t GetNsig() const { return fNsigCut; }
  Int_t    GetNoff() const { return fNoffVTracks; }
  void     SetNsig(Double_t nsig);                 // sets nsig
  void     SetNoff(Int_t noff);                    // sets noff

  //private:
  Double_t Getbnorm(const ANLTrack &t);

private:
  JSFQuickSimParam *fParam;	  //! Pointer to QuickSimParam
  TClonesArray     *fGen;         //! Pointer to GeneratorParticles
  Double_t          fNsigCut;     //  Cut value of Nsig
  Int_t             fNoffVTracks; //  Number of Off-vertex tracks

  ClassDef(ANLVTXTagger, 1) 	// ANLVTXTagger Example
};

#endif
