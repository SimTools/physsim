#ifndef __ANL2DVECTOR__
#define __ANL2DVECTOR__
//*************************************************************************
//* ===================
//*  ANL2DVector Class
//* ===================
//*
//* (Description)
//*    A very primitive lockable 2D vector class.
//* (Requires)
//*	class TVector2
//* 	class Lockable
//* (Provides)
//* 	class ANL2DVector
//* (Update Recored)
//*    2001/02/16  K.Ikematsu	Original version.
//*    2001/07/12  K.Ikematsu   Added DebugPrint() function.
//*
//* $Id$
//*************************************************************************
//
#include <iostream>
#include "TVector.h"
#include "TVector2.h"
#include "Lockable.h"

using namespace std;

//_____________________________________________________________________
//  -----------------
//  Lockable TVector2
//  -----------------
//
class ANL2DVector : public TVector2, public Lockable {
public:
   ANL2DVector(Double_t px=0., Double_t py=0.) : TVector2(px,py) {}
   ANL2DVector(Float_t px, Float_t py=0.) {
     fX = px;
     fY = py;
   }
   ANL2DVector(const TVector &q) {
     fX = q(0);
     fY = q(1);
   }
   ANL2DVector(const TVector2 &q) : TVector2(q) {}
   ANL2DVector(const ANL2DVector &q) : TVector2(q), Lockable(q) {}

   virtual ~ANL2DVector() {}

   inline Double_t & operator()(Int_t i) {
     switch (i) {
     case 1: return fX;
     case 2: return fY;
     default: Error("ANL2DVector::operator()(int)", "bad index (%d) returning &fX",i); return fX;
     }
   }

   inline Double_t operator()(Int_t i) const {
     switch (i) {
     case 1: return fX;
     case 2: return fY;
     default: Error("ANL2DVector::operator()(int)", "bad index (%d) returning 0",i); return 0.;
     }
   }

   inline friend ANL2DVector operator+ (const ANL2DVector &q1) { return q1; }
   inline friend ANL2DVector operator- (const ANL2DVector &q1) {
       ANL2DVector ans(-q1.Px(),-q1.Py()); return ans;
   }

   inline virtual void DebugPrint(const Char_t *opt = "Brief") const {
     cerr << "p    = " << operator()(1) << " " << operator()(2) << endl;
     if (opt == "Detailed") {
       cerr << "ap   = " << Mod()  << endl;
     }
   }

   ClassDef(ANL2DVector,1)  // Lockable 2D vector class
};

#endif
