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
//*
//* $Id$
//*************************************************************************
//
#include "FlavourGetter.h"
//_____________________________________________________________________
//  -------------------
//  FlavourGetter Class
//  -------------------
//
ClassImp(FlavourGetter)

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t   FlavourGetter::operator()(const ANLJet &jet){

  Int_t nbhad = 0;
  Int_t nchad = 0;
  Int_t nhad = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Int_t hadpid = TMath::Abs(GetPrimaryHadronPID(*tp));
    if ( hadpid != 0 ) ++nhad;
    Int_t partonid = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      partonid = hadpid/1000;
    } else {
      partonid = (hadpid/100)%10;
    }
    if ( partonid == 5 ) ++nbhad;
    if ( partonid == 4 ) ++nchad;
  }

  Double_t blike = (Double_t)nbhad/(Double_t)nhad;
  Double_t clike = (Double_t)nchad/(Double_t)nhad;
  if ( blike >= 0.2 && blike >= clike ) return 3;
  else if ( clike >= 0.2 && clike > blike ) return 2;
  else return 1;
}

//_____________________________________________________________________
Int_t FlavourGetter::GetPrimaryHadronPID(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGpid;
}

Int_t FlavourGetter::GetPrimaryHadronSN(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGsn;
}

Int_t FlavourGetter::GetPartonID(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGmsn;
}

//_____________________________________________________________________
void FlavourGetter::SearchPrimaryHadron(const ANLTrack &t){

  JSFCDCTrack *cdctp = t.GetLTKCLTrack()->GetCDC();
  if (!cdctp) {
    //cerr << "No pointer to JSFCDCTrack !" << endl;
    fGsn  = 0;
    fGpid = 0;
    fGmsn = 0;
    return;
  }

  fGsn = 0;
  fGsn = t.GetLTKCLTrack()->GetCDC()->GetGenID();
  fGpid = 0;
  while ( fGsn > 0 ) {
    JSFGeneratorParticle *g = (JSFGeneratorParticle *)fGen->UncheckedAt(fGsn-1);
    fGpid = g->GetID();
    fGmsn = g->GetMother();

    /*
    cerr << "(PID, S.N, M.S.N) = ("
         << fGpid << ","
         << fGsn  << ","
         << fGmsn << ")" << endl;
    */
    fGsn = fGmsn;
  }
  //cerr << "-- Search ended --" << endl;
}
