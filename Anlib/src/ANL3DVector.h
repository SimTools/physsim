#ifndef __ANL3DVECTOR__
#define __ANL3DVECTOR__
//*************************************************************************
//* ===================
//*  ANL3DVector Class
//* ===================
//*
//* (Description)
//*    A very primitive lockable 3D vector class.
//* (Requires)
//* 	class Lockable
//*	class TVector3
//* (Provides)
//* 	class ANL3DVector
//* (Update Recored)
//*    1999/09/13  K.Ikematsu	Original version.
//*    2000/03/23  K.Ikematsu	Added GetNorm method.
//*    2000/03/23  K.Ikematsu	Added GetTheta method.
//*    2000/03/28  K.Ikematsu	Added Acol method.
//*
//*************************************************************************
//
#include <iostream.h>
#include "Lockable.h"
#include "TVector3.h"
//_____________________________________________________________________
//  -----------------
//  Lockable TVector3
//  -----------------
//
class ANL3DVector : public TVector3, public Lockable {
public:
   ANL3DVector(Double_t px=0., Double_t py=0., Double_t pz=0.) 
   			               : TVector3(px,py,pz) {}
   ANL3DVector(Float_t px, Float_t py=0., Float_t pz=0.) {
     TVector3::operator[](0) = px;
     TVector3::operator[](1) = py;
     TVector3::operator[](2) = pz;
   }
   ANL3DVector(const TVector &q) {
     TVector3::operator[](0) = q(1);
     TVector3::operator[](1) = q(2);
     TVector3::operator[](2) = q(3);
   }
   ANL3DVector(const TVector3 &q) : TVector3(q) {}
   ANL3DVector(const ANL3DVector &q) : TVector3(q), Lockable(q) {}

   virtual ~ANL3DVector() {}

   inline Double_t & Px() { return TVector3::operator[](0); }
   inline Double_t & Py() { return TVector3::operator[](1); }
   inline Double_t & Pz() { return TVector3::operator[](2); }
   inline Double_t & X()  { return TVector3::operator[](0); }
   inline Double_t & Y()  { return TVector3::operator[](1); }
   inline Double_t & Z()  { return TVector3::operator[](2); }
   inline Double_t & operator()(Int_t i) {
     return TVector3::operator[](i-1);
   }

   inline Double_t Px() const { return TVector3::operator[](0); }
   inline Double_t Py() const { return TVector3::operator[](1); }
   inline Double_t Pz() const { return TVector3::operator[](2); }
   inline Double_t X()  const { return TVector3::operator[](0); }
   inline Double_t Y()  const { return TVector3::operator[](1); }
   inline Double_t Z()  const { return TVector3::operator[](2); }
   inline Double_t operator()(Int_t i) const {
     return TVector3::operator[](i-1);
   }

   inline friend ANL3DVector operator+ (const ANL3DVector &q1,
				        const ANL3DVector &q2) {
     ANL3DVector ans = q1; ans += q2; return ans;
   }
   inline friend ANL3DVector operator- (const ANL3DVector &q1,
 				        const ANL3DVector &q2) {
     ANL3DVector ans = q1; ans -= q2; return ans;
   }
   inline friend ANL3DVector operator+ (const ANL3DVector &q1) { return q1; }
   inline friend ANL3DVector operator- (const ANL3DVector &q1) {
       ANL3DVector ans(-q1(1),-q1(2),-q1(3)); return ans;
   }
   inline friend Double_t operator* (const ANL3DVector &q1,
				     const ANL3DVector &q2) {
     return ( q1(1)*q2(1) + q1(2)*q2(2) + q1(3)*q2(3) );
   }


   inline friend ANL3DVector operator^ (const ANL3DVector &q1,
				        const ANL3DVector &q2) {
     ANL3DVector ans( q1.Cross(q2) );
     return ans;
   }

   inline Double_t GetPt2()   const { return ( operator()(1)*operator()(1) 
                                       + operator()(2)*operator()(2) ); }
   inline Double_t GetMag2()  const { return ( GetPt2() 
                                       + operator()(3)*operator()(3) ); }
   inline Double_t GetPt()    const { return TMath::Sqrt( GetPt2() );       }
   inline Double_t GetMag()   const { return TMath::Sqrt( GetMag2() );      }
   inline Double_t GetTheta() const { return (180.*(TMath::ACos(CosTheta()))/TMath::Pi()); }
   inline Double_t GetTheta(const ANL3DVector &q) const {
     return (180.*(TMath::ACos(CosTheta(q)))/TMath::Pi());
   }

   inline ANL3DVector GetNorm() const {
     ANL3DVector n(operator()(1)/GetMag(),
		   operator()(2)/GetMag(),
		   operator()(3)/GetMag());
     return n;
   }

   inline Double_t Acol(const ANL3DVector &q) const {
     return (180.*(TMath::Pi()-TMath::ACos(CosTheta(q)))/TMath::Pi());
   }
   inline Double_t Acop(const ANL3DVector &q) const {
     Double_t c = (operator()(1)*q(1)+operator()(2)*q(2))
                  /(GetPt()*q.GetPt());
     return (180.*(TMath::Pi()-TMath::ACos(c))/TMath::Pi());
   }

   inline Double_t CosTheta() const { return TVector3::CosTheta(); }
   inline Double_t CosTheta(const ANL3DVector &q) const {
     return (operator()(1)*q(1)+operator()(2)*q(2)+operator()(3)*q(3))
                  /(GetMag()*q.GetMag());
   }

   inline virtual void DebugPrint(const Char_t *opt = "Brief") const { 
     cerr << "p    = " << operator()(1) << " "
                       << operator()(2) << " " << operator()(3) << endl;
     if (opt == "Detailed") {
       cerr << "pt   = " << GetPt()   << endl;
       cerr << "ap   = " << GetMag()  << endl;
     }
   }

   ClassDef(ANL3DVector,1)  // Lockable 3D vector class
};

#endif
