//*************************************************************************
//* ==================
//*  ANLTrack Classes
//* ==================
//*
//* (Description)
//*    Track Class for JLC analyses
//* (Requires)
//*     class TVector
//*     class JSFLTKCLTrack
//*     class ANL4DVector
//* (Provides)
//*     class ANLTrack
//* (Update Recored)
//*    1999/08/17  K.Ikematsu   Original version.
//*    1999/09/05  K.Ikematsu   Replaced LockableLVector with ANL4DVector.
//*    1999/09/13  K.Ikematsu   Added GetLTKCLTrack Method.
//*    1999/09/13  K.Ikematsu   Added GetTrackPtr Method.
//*    2001/07/09  K.Ikematsu   Changed IsElectron, IsMuon & GetCharge
//*                             method to virtual.
//*
//* $Id$
//*************************************************************************
//
#include "ANLTrack.h"
//_____________________________________________________________________
//  --------------
//  ANLTrack Class
//  --------------
//
ClassImp(ANLTrack)

//*--
//*  Constructors
//*--
ANLTrack::ANLTrack() : ANL4DVector(0.), fTrackPtr(0) {}
ANLTrack::ANLTrack(const TObject *track) :
  ANL4DVector(((JSFLTKCLTrack *)track)->GetPV()), fTrackPtr(track) {}
ANLTrack::ANLTrack(const TVector &pv, const TObject *ptr) :
  ANL4DVector(pv), fTrackPtr(ptr) {}

//*--
//*  Destructor
//*--
ANLTrack::~ANLTrack() {}

//*--
//*  Getters
//*--
Bool_t   ANLTrack::IsElectron() const {
  return ( ((JSFLTKCLTrack *)fTrackPtr)->GetType() == 11 );
}
Bool_t   ANLTrack::IsMuon() const {
  return ( ((JSFLTKCLTrack *)fTrackPtr)->GetType() == 13 );
}
Bool_t   ANLTrack::IsLepton() const {
  return ( IsElectron() || IsMuon() );
}
Double_t ANLTrack::GetCharge() const {
  return ((JSFLTKCLTrack *)fTrackPtr)->GetCharge();
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
TObject *ANLTrack::GetTrackPtr() const {
  return (TObject *)fTrackPtr;
}
JSFLTKCLTrack *ANLTrack::GetLTKCLTrack() const {
  return (JSFLTKCLTrack *)fTrackPtr;
}
