//*************************************************************************
//* ==================
//*  ANLTrack Classes
//* ==================
//*
//* (Description)
//*    Track Class for JLC analyses
//* (Requires)
//*     class TLorentzVector
//*     class Lockable
//*     class ANL4DVector
//* (Provides)
//*     class ANLTrack
//* (Update Recored)
//*    1999/08/17  K.Ikematsu   Original version.
//*    1999/09/05  K.Ikematsu   Replaced LockableLVector with ANL4DVector.
//*    1999/09/13  K.Ikematsu   Added GetLTKCLTrack Method.
//*
//*************************************************************************
#include "ANLTrack.h"
//_____________________________________________________________________
//  --------------
//  ANLTrack Class
//  --------------

ClassImp(ANLTrack)

//*--
//*  Constructors
//*--
ANLTrack::ANLTrack() : ANL4DVector(0.), fTrackPntr(0) {}
#if 0
ANLTrack::ANLTrack(const TObject *track) :
  ANL4DVector(((JSFLTKCLTrack *)track)->GetPV()), fTrackPntr(track) {}
#else
ANLTrack::ANLTrack(const TObject *track) :
  ANL4DVector(((JSFLTKCLTrack *)track)->GetE(),
  	     ((JSFLTKCLTrack *)track)->GetPx(),
  	     ((JSFLTKCLTrack *)track)->GetPy(),
  	     ((JSFLTKCLTrack *)track)->GetPz()), 
  fTrackPntr(track) {}
#endif

//*--
//*  Destructor
//*--
ANLTrack::~ANLTrack() {}

//*--
//*  Getters
//*--
Bool_t   ANLTrack::IsElectron() const {
  return ( ((JSFLTKCLTrack *)fTrackPntr)->GetType() == 11 );
}
Bool_t   ANLTrack::IsMuon() const {
  return ( ((JSFLTKCLTrack *)fTrackPntr)->GetType() == 13 );
}
Bool_t   ANLTrack::IsLepton() const {
  return ( IsElectron() || IsMuon() );
}
Double_t ANLTrack::GetCharge() const {
  return ((JSFLTKCLTrack *)fTrackPntr)->GetCharge();
}
Double_t ANLTrack::GetConeEnergy(const Double_t cth,
				 const TObjArray *tracks) const {
  Double_t e = 0.;
  TIter next(tracks);
  ANLTrack *t;
  while ((t = (ANLTrack *)next())) {
    if ( t == this ) continue;
    if ( CosTheta(*t) > cth ) e += t->E();
  }
  return e;
}
JSFLTKCLTrack *ANLTrack::GetLTKCLTrack() const {
  return (JSFLTKCLTrack *)fTrackPntr;
}
