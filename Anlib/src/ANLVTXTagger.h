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
//* (To do)
//*     Ks, Lambda removal.
//* (Update Recored)
//*    1999/10/09  K.Fujii	Original version.
//*
//*************************************************************************
//
#include "JSFSteer.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "ANLJetFinder.h"
#include "ANLTrack.h"
//_____________________________________________________________________
//  -------------------
//  ANLVTXTagger Class
//  -------------------
//
class ANLVTXTagger {
public:
  ANLVTXTagger(Double_t nsig=3., Int_t noff=3)
            : fNsigCut(nsig), fNoffVTracks(noff) {
    JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
    JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
    fParam = sds->Param();
    fGen   = evt->GetGeneratorParticles();
  }                                                // default constructor

  virtual ~ANLVTXTagger(){}                        // default destructor
  Bool_t operator()(const ANLJet &jet);

private:
  Double_t Getbnorm(const ANLTrack &t);

private:
  JSFQuickSimParam *fParam;	  //! Pointer to QuickSimParam
  TClonesArray     *fGen;         //! Pointer to GeneratorParticles
  Double_t          fNsigCut;     //  Cut value of Nsig
  Int_t             fNoffVTracks; //  Number of Off-vertex tracks

  ClassDef(ANLVTXTagger, 1) 	// ANLVTXTagger Example
};

#endif

