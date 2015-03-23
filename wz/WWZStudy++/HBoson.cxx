//*****************************************************************************
//* =====================
//*  HBoson class
//* =====================
//*  
//* (Description)
//*    A PDT entry representing the Higgs boson
//*
//* (Update Record)
//*    2014/11/21  K.Fujii	Original version.
//*    2015/03/23  K.Fujii	Improved Gamma_H calculations.
//*
//*****************************************************************************

#include "HBoson.h"
#include "GENVVStarDecay.h"

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(HBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class HBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor & D-tor
// --------------------------
HBoson::HBoson(Double_t m,
               Double_t l,
               Double_t a,
               Double_t b,
               Double_t bt)
         : GENPDTEntry("H", // name
                        25, // PID
                        0., // charge
                        0., // spin
                         m),// mass
           fLambda(l),
           fA     (a),
           fB     (b),
           fBtilde(bt)
{
   // --------------------------------------------
   //  Initialize mass
   // --------------------------------------------
   fMass = m;
   // --------------------------------------------
   //  Initialize W/Z decay table
   // --------------------------------------------
   fWMBosonPtr     = new GENPDTWBoson();
   fWPBosonPtr     = new GENPDTWBoson();
   fZ1BosonPtr     = new GENPDTZBoson();
   fZ2BosonPtr     = new GENPDTZBoson();
   //--
   // Full amplitude and helicity selector should be implemented
   // and passed to VV* decay calculator class.
   //--
   fWWAmpUtilsPtr  = new HBoson::AmpUtils(this);
   fWWStarDecayPtr = new GENVVStarDecay(this, fWMBosonPtr, fWPBosonPtr);
   fWWStarDecayPtr->SetAmpUtilsPtr(fWWAmpUtilsPtr);
   fZZAmpUtilsPtr  = new HBoson::AmpUtils(this);
   fZZStarDecayPtr = new GENVVStarDecay(this, fZ1BosonPtr, fZ2BosonPtr);
   fZZStarDecayPtr->SetAmpUtilsPtr(fZZAmpUtilsPtr);
   // --------------------------------------------
   //  Calculate width.
   // --------------------------------------------
   Double_t gam = 0.;
   //--
   // h --> WW*
   //--
   gam = fWWStarDecayPtr->CalculateWidth();
   if (gam > 0.) {
      GENDecayMode *dmp = new GENDecayMode(gam);
      GENPDTEntry  *d1p = new GENPDTWBoson();
      GENPDTEntry  *d2p = new GENPDTWBoson();
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
#if 1
   //--
   // h --> ZZ*
   //--
   gam = fZZStarDecayPtr->CalculateWidth();
   if (gam > 0.) {
      gam *= 0.5; // identical particle factor
      GENDecayMode *dmp = new GENDecayMode(gam);
      GENPDTEntry  *d1p = new GENPDTZBoson();
      GENPDTEntry  *d2p = new GENPDTZBoson();
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
   //--
   // h --> ff
   //--
   Double_t lamb = GetLambdaQCD(kAlphaS, kM_z);
   Double_t alfs = GetAlphaS(lamb, fMass);
   for (Int_t ic=0; ic<2; ic++) {
      for (Int_t it=0; it<2; it++) { 
         for (Int_t ig=0; ig<3; ig++) {
            const Char_t   *name = kName[ic][it][ig];
                  Int_t     pid  = kPID [ic][it][ig];  
                  Double_t  t3   = (1 - 2*it)/2.;
                  Double_t  qf   = kChrg[ic][it][ig];
                  Double_t  spin = 0.5;
                  Double_t  cf   = 2*ic + 1;
                  Double_t  mass = kMass[ic][it][ig];
            GENDecayMode *dmp;
            GENPDTEntry  *d1p, *d2p;
            d1p  = new GENPDTEntry(name, pid, qf, spin, mass, ig+1,  t3, cf);
            d2p  = new GENPDTEntry(name,-pid,-qf, spin, mass, ig+1, -t3, cf);
#if 1
	    // Running mass
	    if (ig == 2 && it == 1 && ic == 1) {
		mass = 2.7645; // mb(MSbar)(mh)
	    } else if (ig == 1 && it == 0 && ic == 1) {
		mass = 0.61614; // mc(MSbar)(mh)
	    } 
	    if (cf > 1.) {
		cf *= (1. + (17./3)*alfs/kPi); // QCD correction.
	    }
#endif
            gam = GammaToFF(mass, cf);
            if (gam == 0.) continue;
            dmp = new GENDecayMode(gam);
            dmp->Add(d1p);
            dmp->Add(d2p);
            Add(dmp); 
         }
      }
   }
   //--
   // h --> gg
   //--
   gam = GammaToGG();
   if (gam > 0.) {
      GENDecayMode *dmp = new GENDecayMode(gam);
      GENPDTEntry  *d1p = new GENPDTGluon();
      GENPDTEntry  *d2p = new GENPDTGluon();
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
   //--
   // h --> gamma gamma
   //--
#endif
}

HBoson::~HBoson()
{
   delete fWWStarDecayPtr;
   delete fZZStarDecayPtr;
   delete fWMBosonPtr;
   delete fWPBosonPtr;
   delete fZ1BosonPtr;
   delete fZ2BosonPtr;
   delete fWWAmpUtilsPtr;
   delete fZZAmpUtilsPtr;
}

Double_t HBoson::GammaToFF(Double_t mf,
                           Double_t cf)
{
   Double_t gf   = kSqh*kPi*kAlpha/(kSin2W*kCos2W*kM_z*kM_z);
   Double_t fact = kSqh*gf/k4Pi;
   Double_t ef   = fMass/2.;
   Double_t beta = (ef-mf)*(ef+mf);
   if (beta <= 0.) return 0.;
   else            beta = TMath::Sqrt(beta)/ef;
   return fact*mf*mf*fMass*TMath::Power(beta,3)*cf;
}

Double_t HBoson::GammaToGG()
{
   Double_t lamb  = GetLambdaQCD(kAlphaS, kM_z);
   Double_t alfs  = GetAlphaS(lamb, fMass/2);
#if 0
   Double_t mt    = kMass[1][0][2];
#else
   Double_t mt    = 169.22; // MSbar [GeV] at mh = 125GeV
#endif
   Double_t fact  = TMath::Power(alfs*kGw/kM_w,2)
	           *TMath::Power(fMass/k4Pi,3)/8.; 
   Double_t t     = 4.*TMath::Power(mt/fMass,2);
   Double_t gam   = fact * TMath::Power(abs(Ff(t)),2);
   return gam;
}

Complex_t HBoson::F(Double_t t)
{
   Complex_t ans;
   static const Complex_t kI(0.,1.); 
   if (t >= 1.) {
      ans = TMath::Power(TMath::ASin(TMath::Sqrt(1/t)),2);
   } else {
      ans = TMath::Log((1+TMath::Sqrt(1-t))/(1-TMath::Sqrt(1-t))) - kI*kPi;
      ans = -(1./4.)*ans*ans;
   }
   return ans;
}

Complex_t HBoson::Fs(Double_t t) // J=0 loop
{
   return t*(1.-t*F(t));
}

Complex_t HBoson::Ff(Double_t t) // J=1/2 loop
{
   return -2.*t*(1.+(1.-t)*F(t));
}

Complex_t HBoson::Fv(Double_t t) // J=1 loop
{
   return 2. + 3.*t + 3.*t*(2.-t)*F(t);
}

Double_t HBoson::GetLambdaQCD(Double_t alfs, Double_t q)
{
   Int_t    nf    = 5; // 5 flavors
   Double_t beta0 = (33.-2.*nf)/6;
   Double_t beta1 = (153.-19.*nf)/12;
   Double_t a     = alfs/kPi;
   Double_t ff1   = -1./beta0/a;
   Double_t ff2   = beta1/(beta0*beta0)*TMath::Log(2./beta0*(1/a+beta1/beta0));
   Double_t ans   = q*TMath::Exp(ff1*ff2);
   return ans;
}

Double_t HBoson::GetAlphaS(Double_t lambda, Double_t q)
{
   Int_t    nf    = 5; // 5 flavors
   Double_t beta0 = 11.-2.*nf/3.;
   Double_t beta1 = 2*(51.-19.*nf/3.);
   Double_t QovL = q/lambda;
   Double_t ans  = 2.*TMath::Log(QovL) 
	           + (beta1/(beta0*beta0))*TMath::Log(2*TMath::Log(QovL));
            ans  = k4Pi/(beta0*ans);
   cerr << "Lambda QCD = " << lambda << endl;
   cerr << " Q = " << q << " alfs = " << ans << endl;
   return ans;
}
//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t HBoson::AmpUtils::FullAmplitude(ANL4DVector &k,
                                          ANL4DVector  p[4],
                                          Double_t     m[4],
                                          Int_t        helFinal[4],
                                          GENPDTEntry *v1p,
                                          GENPDTEntry *v2p,
                                          GENPDTEntry *f1p,
                                          GENPDTEntry *f2p,
                                          GENPDTEntry *f3p,
                                          GENPDTEntry *f4p)
{
   Double_t mv1    = v1p->GetMass();
   Double_t gamv1  = v1p->GetWidth();
   Double_t mv2    = v2p->GetMass();
   Double_t gamv2  = v2p->GetWidth();

   Double_t ghvv   = kGw*kM_w;
   Double_t glv1   = -kGw*kSqh;
   Double_t grv1   = 0;
   Double_t glv2   = -kGw*kSqh;
   Double_t grv2   = 0;

   Double_t g1     = ghvv + 2 * mv1 * mv2 * (fHBosonPtr->fA/fHBosonPtr->fLambda);
   Double_t g2     = -2 * (fHBosonPtr->fB/fHBosonPtr->fLambda);
   Double_t g3     = -4 * (fHBosonPtr->fBtilde/fHBosonPtr->fLambda);

   if (v1p->GetName().Contains("Z")) {
      Double_t qf     = f1p->GetCharge();
      Double_t t3f    = f1p->GetISpin();
               glv1   = -kGz*(t3f - qf*kSin2W);
               grv1   = -kGz*(    - qf*kSin2W);
               qf     = f3p->GetCharge();
               t3f    = f3p->GetISpin();
               glv2   = -kGz*(t3f - qf*kSin2W);
               grv2   = -kGz*(    - qf*kSin2W);
	       ghvv   = kGz*kM_z;
#if 1
               g1     = ghvv + 2 * mv1 * mv2 * (fHBosonPtr->fA/fHBosonPtr->fLambda);
               g2     = -2 * (fHBosonPtr->fB/fHBosonPtr->fLambda);
               g3     = -4 * (fHBosonPtr->fBtilde/fHBosonPtr->fLambda);
#else
               g1     = ghvv;
               g2     = 0.;
               g3     = 0.;
#endif
   }

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion f1 (p[0], m[0], helFinal[0], +1, kIsOutgoing);
   HELFermion f2b(p[1], m[1], helFinal[1], -1, kIsIncoming);
   HELVector  v1f(f2b, f1, glv1, grv1, mv1, gamv1);

   HELFermion f3 (p[2], m[2], helFinal[2], +1, kIsOutgoing);
   HELFermion f4b(p[3], m[3], helFinal[3], -1, kIsIncoming);
   HELVector  v2f(f4b, f3, glv2, grv2, mv2, gamv2);

   HELScalar  h(k);

#ifndef ANOM_WWH
   Complex_t amp = HELVertex(v1f, v2f, h, ghvv);
#else
   Complex_t amp = HELVertex(v1f, v2f, h, g1, g2, g3);
#endif

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void HBoson::AmpUtils::SelectHelicities(Double_t  helcombf,
                                        Int_t     helFinal[4],
                                        Double_t &weight)
{
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1}};
   weight = 1.;
   Int_t combf = (Int_t)(helcombf*kNf);
         combf = TMath::Min(combf, kNf-1);
   helFinal  [0] = kFHelComb[combf][0];
   helFinal  [1] = kFHelComb[combf][1];
   helFinal  [2] = kFHelComb[combf][2];
   helFinal  [3] = kFHelComb[combf][3];
   weight *= kNf;
}
