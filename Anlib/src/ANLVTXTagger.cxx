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
#include "ANLVTXTagger.h"
//_____________________________________________________________________
//  ------------------
//  ANLVTXTagger Class
//  ------------------
//
ClassImp(ANLVTXTagger)

//_____________________________________________________________________
//*--
//*  Setters
//*--
void ANLVTXTagger::SetNsig(Double_t nsig) {
  if (nsig != fNsigCut) {
    fNsigCut = nsig;
  }
}

void ANLVTXTagger::SetNoff(Int_t noff) {
  if (noff != fNoffVTracks) {
    fNoffVTracks = noff;
  }
}

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
    if (bnorm > fNsigCut) noffv++;
  }

  if (noffv >= fNoffVTracks) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________
Double_t ANLVTXTagger::Getbnorm(const ANLTrack &t){
  JSFCDCTrack *cdctp = t.GetLTKCLTrack()->GetCDC();
  if (!cdctp) { 
    return 0.;
  } else {

    Int_t gsn  = t.GetLTKCLTrack()->GetCDC()->GetGenID();

    // Warnig!! TClonesArray *fGen starts from i=0.
    JSFGeneratorParticle *g = (JSFGeneratorParticle *)fGen->UncheckedAt(gsn-1);
    //                                                                  ^^^^^
    Int_t gmsn = g->GetMother();

    // If this generator particle comes from SpringParton directly,
    // we cannot get pointer to JSFGeneratorParticle.
    if (gmsn < 0) {
      return 0.;
    }
    JSFGeneratorParticle *gm = (JSFGeneratorParticle *)fGen->UncheckedAt(gmsn-1);
    Int_t gmpid = gm->GetID();

    if (TMath::Abs(gmpid) == 310 || TMath::Abs(gmpid) == 3122
        || TMath::Abs(gmpid) == 3112 || TMath::Abs(gmpid) == 3222) {
      return -99999.;
    }

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
