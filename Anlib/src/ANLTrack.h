#ifndef __ANLTRACK__
#define __ANLTRACK__
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
#include <iostream.h>
#include "TObjArray.h"
#include "JSFLTKCLTrack.h"
#include "ANL4DVector.h"
//_____________________________________________________________________
//  --------------
//  ANLTrack Class
//  --------------

class ANLTrack : public ANL4DVector {
public:
  ANLTrack();
  ANLTrack(const TObject *track);
  virtual ~ANLTrack();

  Bool_t   IsElectron() const;
  Bool_t   IsMuon() const;
  Bool_t   IsLepton() const;
  Double_t GetCharge() const;
  Double_t GetConeEnergy(const Double_t cth, const TObjArray *tracks) const;
  JSFLTKCLTrack *GetLTKCLTrack() const;

private:
  const TObject * fTrackPntr;   // Pointer to JSFLTKCLTrack

  ClassDef(ANLTrack,1)   // ANLTrack class
};

#endif
