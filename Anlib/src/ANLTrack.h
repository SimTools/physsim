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
#include <iostream.h>
#include "TVector.h"
#include "TObjArray.h"
#include "JSFLTKCLTrack.h"
#include "ANL4DVector.h"
//_____________________________________________________________________
//  --------------
//  ANLTrack Class
//  --------------
//
class ANLTrack : public ANL4DVector {
public:
  ANLTrack();
  ANLTrack(const TObject *track);
  ANLTrack(const TVector &pv, const TObject *ptr);
  virtual ~ANLTrack();

  virtual Bool_t   IsElectron() const;
  virtual Bool_t   IsMuon() const;
  Bool_t   IsLepton() const;
  virtual Double_t GetCharge() const;
  Double_t GetConeEnergy(const Double_t cth, const TObjArray *tracks) const;
  TObject *GetTrackPtr() const;
  JSFLTKCLTrack *GetLTKCLTrack() const;

private:
  const TObject * fTrackPtr;   // Pointer to JSFLTKCLTrack

  ClassDef(ANLTrack,1)   // ANLTrack class
};

#endif
