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
#include "ANLVTXTagger.h"
//_____________________________________________________________________
//  ------------------
//  ANLVTXTagger Class
//  ------------------
//
ClassImp(ANLVTXTagger)

//_____________________________________________________________________
//*--
//*  Getters
//*--
Bool_t   ANLVTXTagger::operator()(const ANLJet &jet){

  Int_t noffv = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Double_t bnorm = Getbnorm(*tp);
    if ( bnorm > fNsigCut ) noffv++;
  }

  if ( noffv >= fNoffVTracks) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________
Double_t ANLVTXTagger::Getbnorm(const ANLTrack &t){
  JSFCDCTrack *cdctp = t.GetLTKCLTrack()->GetCDC();
  if (!cdctp) { 
    return 0.;
  } else {
    JSFCDCTrack cdct(*cdctp);
    cdct.MovePivotToIP(fParam);
    Float_t helix[5];
    cdct.GetHelix(helix);
    Double_t err[15];
    cdct.GetError(err);
    Double_t dr   = helix[0];
    Double_t dz   = helix[3];
    Double_t drdr = err[0];
    Double_t dzdz = err[9];
    return TMath::Sqrt(dr*dr/drdr + dz*dz/dzdz);
  }
}
