//*****************************************************************************
//* =====================
//*  HELLib 
//* =====================
//*  
//* (Description)
//*    Class library for helicity amplitude calculations.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version based on HELAS.
//*
//*****************************************************************************

#include "HELLib.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__

  
ClassImp(HELFermion)
ClassImp(HELVector)
ClassImp(HELScalar)
ClassImp(HELVertex)

//-----------------------------------------------------------------------------
// ==============================
//  class HELFermion
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------------------
// IXXXXX() and OXXXXX()
//----------------------
HELFermion::HELFermion(const ANL4DVector &p,
                             Double_t     m,
                             Int_t        hel,
                             Int_t        nsf,
                             Bool_t       isincom)
          : TVectorC(4),
            fP(nsf*p.E(), nsf*p.Px(), nsf*p.Py(), nsf*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSF(nsf),
            fIsIncoming(isincom)
{
   Double_t  sf[2], sfomeg[2], omega[2], pp, pp3, sqpop3, sqm;
   Complex_t chi[2];
   Int_t nh = hel*nsf;
   Int_t in = isincom ? 1 : -1;
   if (m == 0.) {
      sqpop3 = nsf*TMath::Sqrt(TMath::Max(p(0)+p(3), 0.));
      chi[0] = sqpop3;
      if (sqpop3 == 0.) chi[1] = -hel * TMath::Sqrt(2.*p(0));
      else              chi[1] = Complex_t(nh*p(1), in*p(2))/sqpop3;
      Int_t iu = (1-in)/2;
      Int_t id = (1+in)/2;
      if (nh == 1*in) {
         (*this)[0] = 0.;
         (*this)[1] = 0.;
         (*this)[2] = chi[iu];
         (*this)[3] = chi[id];
      } else {
         (*this)[0] = chi[id];
         (*this)[1] = chi[iu];
         (*this)[2] = 0.;
         (*this)[3] = 0.;
      }
   } else {
      pp = TMath::Min(p(0), p.Vect().Mag());
      if (pp == 0.) {
         sqm = TMath::Sqrt(m);
         Int_t ip =    (1+in*nh)/2;
         Int_t im = in*(1-in*nh)/2;
         (*this)[0] = ip       * sqm;
         (*this)[1] = im * nsf * sqm;
         (*this)[2] = ip * nsf * sqm;
         (*this)[3] = im       * sqm;
      } else {
         sf[0] = (1+nsf+(1-nsf)*nh)*0.5;
         sf[1] = (1+nsf-(1-nsf)*nh)*0.5;
         omega[0] = TMath::Sqrt(p(0)+pp);
         omega[1] = m/omega[0];
         Int_t ip = (1+nh)/2;
         Int_t im = (1-nh)/2;
         sfomeg[0] = sf[0] * omega[ip];
         sfomeg[1] = sf[1] * omega[im];
         pp3    = TMath::Max(pp+p(3), 0.);
         chi[0] = TMath::Sqrt(pp3*0.5/pp);
         if (pp3 == 0.) chi[1] = -nh;
         else           chi[1] = Complex_t(nh*p(1), in*p(2))/TMath::Sqrt(2.*pp*pp3);
         Int_t iu = (1-in)/2;
         Int_t id = (1+in)/2;
         (*this)[0] = sfomeg[iu] * chi[im];
         (*this)[1] = sfomeg[iu] * chi[ip];
         (*this)[2] = sfomeg[id] * chi[im];
         (*this)[3] = sfomeg[id] * chi[ip];
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (isincom ? "fin" : "fot")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " nsf*fP = (" << nsf*fP.E () << ", "
                                 << nsf*fP.Px() << ", "
                                 << nsf*fP.Py() << ", "
                                 << nsf*fP.Pz() << ") " << endl;
#endif
}

//----------------------
// FVIXXX() or FVOXXX()
//----------------------
HELFermion::HELFermion(const HELFermion &f,
                       const HELVector  &v,
                             Double_t    gl,
                             Double_t    gr,
                             Double_t    m,
                             Double_t    gam)
          : TVectorC(4),
            fP((f.fIsIncoming ? f.fP - v.fP : f.fP + v.fP)),
	    fM(m),
            fHel(0),
            fNSF(0),
            fIsIncoming(f.fIsIncoming)
{
   static const Complex_t kI(0., 1.);
   
   Double_t  pf2 = fP.Mag2();
   Complex_t d   = -1./Complex_t(pf2-m*m, TMath::Max(TMath::Sign(m*gam, pf2), 0.));
   if (fIsIncoming) {
      Complex_t sl1 = (v[0] +    v[3])*f[0]
                    + (v[1] - kI*v[2])*f[1];
      Complex_t sl2 = (v[1] + kI*v[2])*f[0]
                    + (v[0] -    v[3])*f[1];
      if (gr == 0.) {
         (*this)[0] = gl * ((fP(0) - fP(3))*sl1 - (fP(1) - fP(2)*kI)*sl2)*d;
         (*this)[1] = gl * ((fP(0) + fP(3))*sl2 - (fP(1) + fP(2)*kI)*sl1)*d;
         (*this)[2] = gl * m * sl1 * d;
         (*this)[3] = gl * m * sl2 * d;
      } else {
         Complex_t sr1 =   (v[0] -    v[3])*f[2]
                         - (v[1] - kI*v[2])*f[3];
         Complex_t sr2 = - (v[1] + kI*v[2])*f[2]
                         + (v[0] +    v[3])*f[3];
         (*this)[0] = ( gl * ((fP(0) - fP(3))*sl1 - (fP(1) - fP(2)*kI)*sl2)
                      + gr * m * sr1) * d;
         (*this)[1] = ( gl * ((fP(0) + fP(3))*sl2 - (fP(1) + fP(2)*kI)*sl1)
                      + gr * m * sr2) * d;
         (*this)[2] = ( gr * ((fP(0) + fP(3))*sr1 + (fP(1) - fP(2)*kI)*sr2)
                      + gl * m * sl1) * d;
         (*this)[3] = ( gr * ((fP(0) - fP(3))*sr2 + (fP(1) + fP(2)*kI)*sr1)
                      + gl * m * sl2) * d;
      }
   } else {
      Complex_t sl1 = (v[0] +    v[3])*f[2]
                    + (v[1] + kI*v[2])*f[3];
      Complex_t sl2 = (v[1] - kI*v[2])*f[2]
                    + (v[0] -    v[3])*f[3];
      if (gr == 0.) {
         (*this)[0] = gl * m * sl1 * d;
         (*this)[1] = gl * m * sl2 * d;
         (*this)[2] = gl * ((fP(0) - fP(3))*sl1 - (fP(1) + fP(2)*kI)*sl2)*d;
         (*this)[3] = gl * ((fP(0) + fP(3))*sl2 - (fP(1) - fP(2)*kI)*sl1)*d;
      } else {
         Complex_t sr1 =   (v[0] -    v[3])*f[0]
                         - (v[1] - kI*v[2])*f[1];
         Complex_t sr2 = - (v[1] + kI*v[2])*f[0]
                         + (v[0] +    v[3])*f[1];
         (*this)[0] = ( gr * ((fP(0) + fP(3))*sr1 + (fP(1) + fP(2)*kI)*sr2)
                      + gl * m * sl1) * d;
         (*this)[1] = ( gr * ((fP(0) - fP(3))*sr2 + (fP(1) - fP(2)*kI)*sr1)
                      + gl * m * sl2) * d;
         (*this)[2] = ( gl * ((fP(0) - fP(3))*sl1 - (fP(1) + fP(2)*kI)*sl2)
                      + gr * m * sr1) * d;
         (*this)[3] = ( gl * ((fP(0) + fP(3))*sl2 - (fP(1) - fP(2)*kI)*sl1)
                      + gr * m * sr2) * d;
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (fIsIncoming ? "fvi" : "fvo")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELVector
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// VXXXXX()
//----------
HELVector::HELVector(const ANL4DVector &p,
                           Double_t     m,
                           Int_t        hel,
                           Int_t        nsv)
          : TVectorC(4),
            fP(nsv*p.E(), nsv*p.Px(), nsv*p.Py(), nsv*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSV(nsv)
{
   static const Double_t kSqh = TMath::Sqrt(0.5);
   Double_t  hel0, pt, pt2, pp, pzpt, emp;
   Int_t     nsvahl;
   nsvahl = nsv*TMath::Abs(hel);
   pt2 = p.GetPt2();
   pp  = TMath::Min(p(0), p.Vect().Mag());
   pt  = TMath::Min(pp, TMath::Sqrt(pt2));
   if (m == 0.) {
      pp = p(0);
      pt = p.GetPt();
      (*this)[0] = 0.;
      (*this)[3] = hel*pt/pp*kSqh;
      if (pt != 0.) {
         pzpt = p(3)/(pp*pt)*kSqh*hel;
         (*this)[1] = Complex_t(-p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
         (*this)[2] = Complex_t(-p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
      } else {
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
      }
   } else {
      hel0 = 1. - TMath::Abs(hel);
      if (pp == 0.) {
         (*this)[0] = 0.;
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*kSqh);
         (*this)[3] = hel0;
      } else {
         emp = p(0)/(m*pp);
	 (*this)[0] = hel0*pp/m;
         (*this)[3] = hel0*p(3)*emp + hel*pt/pp*kSqh;
	 if (pt != 0.) {
            pzpt = p(3)/(pp*pt)*kSqh*hel;
	    (*this)[1] = Complex_t(hel0*p(1)*emp - p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
	    (*this)[2] = Complex_t(hel0*p(2)*emp - p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
	 } else {
	    (*this)[1] = -hel*kSqh;
	    (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
	 }
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (nsv < 0 ? "vin" : "vot")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " nsv*fP = (" << nsv*fP.E () << ", "
                                 << nsv*fP.Px() << ", "
                                 << nsv*fP.Py() << ", "
                                 << nsv*fP.Pz() << ") " << endl;
#endif
}

//----------
// JIOXXX()
//----------
HELVector::HELVector(const HELFermion &fin,
                     const HELFermion &fout,
                           Double_t    glv,
                           Double_t    grv,
                           Double_t    mv,
                           Double_t    gamv)
          : TVectorC(4),
            fP(fout.fP - fin.fP),
	    fM(mv),
            fHel(0),
            fNSV(0)
{
   Complex_t c0, c1, c2, c3, cs, d;
   Double_t  q2, vm2, dd;
   q2  = fP.Mag2();
   vm2 = mv*mv;
   if (mv == 0.) {
      dd = 1./q2;
      if (grv == 0.) { // purely left-handed
         dd *= glv;
	 (*this)[0] = ( fout[2]*fin[0] + fout[3]*fin[1]) * dd;
	 (*this)[1] = (-fout[2]*fin[1] - fout[3]*fin[0]) * dd;
	 (*this)[2] = ( fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., dd);
	 (*this)[3] = (-fout[2]*fin[0] + fout[3]*fin[1]) * dd;
      } else {
         (*this)[0] = (  glv * ( fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] + fout[1]*fin[3])) * dd;
         (*this)[1] = (- glv * ( fout[2]*fin[1] + fout[3]*fin[0])
                       + grv * ( fout[0]*fin[3] + fout[1]*fin[2])) * dd;
         (*this)[2] = (  glv * ( fout[2]*fin[1] - fout[3]*fin[0])
                       + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., dd);
         (*this)[3] = (  glv * (-fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] - fout[1]*fin[3])) * dd;
      }
   } else {
      d = 1./Complex_t(q2 - vm2, TMath::Max(TMath::Sign(mv*gamv, q2), 0.));
      if (grv == 0.) { // purely left-handed
         d *= glv;
	 c0 =  fout[2]*fin[0] + fout[3]*fin[1];
	 c1 = -fout[2]*fin[1] - fout[3]*fin[0];
	 c2 = (fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., 1.);
	 c3 = -fout[2]*fin[0] + fout[3]*fin[1];
      } else {
         c0 =   glv * ( fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] + fout[1]*fin[3]);
         c1 = - glv * ( fout[2]*fin[1] + fout[3]*fin[0])
              + grv * ( fout[0]*fin[3] + fout[1]*fin[2]);
         c2 =  (glv * ( fout[2]*fin[1] - fout[3]*fin[0])
              + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., 1.);
         c3 =   glv * (-fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] - fout[1]*fin[3]);
      }
      cs = (fP(0)*c0 - fP(1)*c1 - fP(2)*c2 - fP(3)*c3)/vm2;
      (*this)[0] = (c0 - cs*fP(0))*d;
      (*this)[1] = (c1 - cs*fP(1))*d;
      (*this)[2] = (c2 - cs*fP(2))*d;
      (*this)[3] = (c3 - cs*fP(3))*d;
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " jio =(" << (*this)[0] << ", "
                             << (*this)[1] << ", "
                             << (*this)[2] << ", "
                             << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}

//----------
// J3XXXX()
//----------
HELVector::HELVector(const HELFermion &fin,
                     const HELFermion &fout,
                           Double_t    gla,
                           Double_t    gra,
                           Double_t    glz,
                           Double_t    grz,
                           Double_t    mz,
                           Double_t    gamz)
          : TVectorC(4),
            fP(fout.fP - fin.fP),
	    fM(mz),
            fHel(0),
            fNSV(0)
{
   Double_t q2  = fP.Mag2();
   Double_t zm2 = mz*mz;
   Double_t zmw = mz*gamz;

   Double_t  da   = 1./q2;
   Double_t  ww   = TMath::Max(TMath::Sign(zmw, q2), 0.);
   Complex_t dz   = 1./Complex_t(q2 - zm2, ww);
   Complex_t ddif = Complex_t(-zm2, ww) * da * dz;
   Double_t  cw   = 1./TMath::Sqrt(1.+ TMath::Power(grz/gra,2));
   Double_t  sw   = TMath::Sqrt((1.-cw)*(1+cw));
   Double_t  gn   = gra * sw;
   Double_t  gz3l = glz * cw;
   Double_t  ga3l = gla * sw;

   static const Complex_t kI(0., 1.);
   Complex_t c0l  =   fout[2] * fin[0] + fout[3] * fin[1];
   Complex_t c0r  =   fout[0] * fin[2] + fout[1] * fin[3];
   Complex_t c1l  = -(fout[2] * fin[1] + fout[3] * fin[0]);
   Complex_t c1r  =   fout[0] * fin[3] + fout[1] * fin[2];
   Complex_t c2l  =  (fout[2] * fin[1] - fout[3] * fin[0])*kI;
   Complex_t c2r  = (-fout[0] * fin[3] + fout[1] * fin[2])*kI;
   Complex_t c3l  =  -fout[2] * fin[0] + fout[3] * fin[1];
   Complex_t c3r  =   fout[0] * fin[2] - fout[1] * fin[3];
   Complex_t csl  = (-fP(0)*c0l+fP(1)*c1l+fP(2)*c2l+fP(3)*c3l)/zm2; 
   Complex_t csr  = (-fP(0)*c0r+fP(1)*c1r+fP(2)*c2r+fP(3)*c3r)/zm2; 

   (*this)[0] = gz3l * dz * (c0l  + csl * fP(0)) + ga3l * c0l *da
               + gn * (c0r * ddif - csr * fP(0)*dz);
   (*this)[1] = gz3l * dz * (c1l  + csl * fP(1)) + ga3l * c1l *da
               + gn * (c1r * ddif - csr * fP(1)*dz);
   (*this)[2] = gz3l * dz * (c2l  + csl * fP(2)) + ga3l * c2l *da
               + gn * (c2r * ddif - csr * fP(2)*dz);
   (*this)[3] = gz3l * dz * (c3l  + csl * fP(3)) + ga3l * c3l *da
               + gn * (c3r * ddif - csr * fP(3)*dz);
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " j3 =(" << (*this)[0] << ", "
                            << (*this)[1] << ", "
                            << (*this)[2] << ", "
                            << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELScalar
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// SXXXXX()
//----------
HELScalar::HELScalar(const ANL4DVector &p,
                           Int_t        nss)
          : Complex_t(1.,0.),
            fP(nss*p.E(), nss*p.Px(), nss*p.Py(), nss*p.Pz()),
            fNSS(nss)
{
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELVertex
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// IOVXXX()
//----------
HELVertex::HELVertex(const HELFermion &in,
                     const HELFermion &out,
                     const HELVector  &v,
                           Double_t    gl,
                           Double_t    gr)
{
   
   *this    =  gl * ( (out[2]*in[0] + out[3]*in[1]) * v[0]
                     +(out[2]*in[1] + out[3]*in[0]) * v[1]
                     -(out[2]*in[1] - out[3]*in[0]) * v[2] * Complex_t(0., 1.)
                     +(out[2]*in[0] - out[3]*in[1]) * v[3] );
   if (gr != 0.) {
      *this += gr * ( (out[0]*in[2] + out[1]*in[3]) * v[0]
                     -(out[0]*in[3] + out[1]*in[2]) * v[1]
                     +(out[0]*in[3] - out[1]*in[2]) * v[2] * Complex_t(0., 1.)
                     -(out[0]*in[2] - out[1]*in[3]) * v[3] );
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " iov =(" << (*this) << endl;
#endif
}

//----------
// VVSXXX()
//----------
HELVertex::HELVertex(const HELVector &v1,
                     const HELVector &v2,
                     const HELScalar &sc,
                           Double_t    g)
{
   *this = g * sc * (v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3]);
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " vvs =(" << (*this) << endl;
#endif
}

//----------
// VVVXXX()
//----------
HELVertex::HELVertex(const HELVector &wm,
                     const HELVector &wp,
                     const HELVector &w3,
                           Double_t    g)
{
   Complex_t v12 = wm[0]*wp[0] - wm[1]*wp[1] - wm[2]*wp[2] - wm[3]*wp[3];
   Complex_t v23 = wp[0]*w3[0] - wp[1]*w3[1] - wp[2]*w3[2] - wp[3]*w3[3];
   Complex_t v31 = w3[0]*wm[0] - w3[1]*wm[1] - w3[2]*wm[2] - w3[3]*wm[3];
   Complex_t xv1 = 0.;
   Complex_t xv2 = 0.;
   Complex_t xv3 = 0.;
   if (abs(wm[0]) != 0. 
       && abs(wm[0]) >= TMath::Max(TMath::Max(abs(wm[1]),
                                              abs(wm[2])),
                                              abs(wm[3]))*1.e-1) {
      xv1 = wm.fP(0)/wm[0];
   }
   if (abs(wp[0]) != 0. 
       && abs(wp[0]) >= TMath::Max(TMath::Max(abs(wp[1]),
                                              abs(wp[2])),
                                              abs(wp[3]))*1.e-1) {
      xv2 = wp.fP(0)/wp[0];
   }
   if (abs(w3[0]) != 0. 
       && abs(w3[0]) >= TMath::Max(TMath::Max(abs(w3[1]),
                                              abs(w3[2])),
                                              abs(w3[3]))*1.e-1) {
      xv3 = w3.fP(0)/w3[0];
   }
   Complex_t p12 = (wm.fP(0) - xv1*wm[0])*wp[0] - (wm.fP(1) - xv1*wm[1])*wp[1]
                 - (wm.fP(2) - xv1*wm[2])*wp[2] - (wm.fP(3) - xv1*wm[3])*wp[3];
   Complex_t p13 = (wm.fP(0) - xv1*wm[0])*w3[0] - (wm.fP(1) - xv1*wm[1])*w3[1]
                 - (wm.fP(2) - xv1*wm[2])*w3[2] - (wm.fP(3) - xv1*wm[3])*w3[3];
   Complex_t p21 = (wp.fP(0) - xv2*wp[0])*wm[0] - (wp.fP(1) - xv2*wp[1])*wm[1]
                 - (wp.fP(2) - xv2*wp[2])*wm[2] - (wp.fP(3) - xv2*wp[3])*wm[3];
   Complex_t p23 = (wp.fP(0) - xv2*wp[0])*w3[0] - (wp.fP(1) - xv2*wp[1])*w3[1]
                 - (wp.fP(2) - xv2*wp[2])*w3[2] - (wp.fP(3) - xv2*wp[3])*w3[3];
   Complex_t p31 = (w3.fP(0) - xv3*w3[0])*wm[0] - (w3.fP(1) - xv3*w3[1])*wm[1]
                 - (w3.fP(2) - xv3*w3[2])*wm[2] - (w3.fP(3) - xv3*w3[3])*wm[3];
   Complex_t p32 = (w3.fP(0) - xv3*w3[0])*wp[0] - (w3.fP(1) - xv3*w3[1])*wp[1]
                 - (w3.fP(2) - xv3*w3[2])*wp[2] - (w3.fP(3) - xv3*w3[3])*wp[3];
   *this = - (v12*(p13-p23) + v23*(p21-p31) + v31*(p32-p12)) * g;
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " vvv =(" << (*this) << endl;
#endif
}