//*****************************************************************************
//* =====================
//*  XBoson class
//* =====================
//*  
//* (Description)
//*    A PDT entry representing X boson
//*
//* (Update Record)
//*    2014/02/01  K.Fujii	Original version.
//*
//*****************************************************************************

#include "XBoson.h"
#include "GENVVStarDecay.h"

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(XBoson)

//-----------------------------------------------------------------------------
// ==============================
//  class XBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor & D-tor
// --------------------------
XBoson::XBoson(Double_t m,
               Double_t fhwz,
               Double_t fhwa)
         : fFhwz(fhwz),
           fFhwa(fhwa)
{
   // --------------------------------------------
   //  Initialize mass
   // --------------------------------------------
   fMass = m;
   // --------------------------------------------
   //  Initialize W/Z decay table
   // --------------------------------------------
   fWBosonPtr      = new GENPDTWBoson();
   fZBosonPtr      = new GENPDTZBoson();
   //--
   // Full amplitude and helicity selector should be implemented
   // and passed to VV* decay calculator class.
   //--
   fAmpUtilsPtr    = new XBoson::AmpUtils(this);
   fVVStarDecayPtr = new GENVVStarDecay(this, fWBosonPtr, fZBosonPtr);
   fVVStarDecayPtr->SetAmpUtilsPtr(fAmpUtilsPtr);
   // --------------------------------------------
   //  Calculate width.
   // --------------------------------------------
   Double_t gam = 0.;
   //--
   // H+ --> W+Z
   //--
   if (fMass > kM_w + kM_z) { // 2-body H+ -> W+ Z decay is open.
      Double_t amp2  = 3.* TMath::Power(kGw * kM_w * fFhwz,2);
      Double_t betab = BetaBar(TMath::Power(kM_w/fMass,2),TMath::Power(kM_z/fMass,2));
      gam = (1./(2.*fMass)) * amp2 * (betab/(2.*k4Pi));
   } else { // 2-body H+ -> W+ Z decay is not open, need to do numerial integration.
      gam = fVVStarDecayPtr->CalculateWidth();
   }
   if (gam > 0.) {
      GENDecayMode *dmp = new GENDecayMode(gam);
      GENPDTEntry  *d1p = new GENPDTWBoson();
      GENPDTEntry  *d2p = new GENPDTZBoson();
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
   }
   //--
   // H+ --> W+ gamma
   //--
   if (fMass > kM_w) {
      Double_t amp2  = 3.* TMath::Power(kGw * kM_w * fFhwa,2);
      Double_t betab = BetaBar(TMath::Power(kM_w/fMass,2),0.);
      gam = (1./(2.*fMass)) * amp2 * (betab/(2.*k4Pi));
   } else {
      std::cerr << ">>>> ERROR in XBoson::C-tor >>>>>>>>>>>>>>>>" << std::endl
                << "   error type = case not supported."          << std::endl
                << "   fMass = " << fMass << " < kM_z = " << kM_z << std::endl
                << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  }
  if (gam > 0.) {
      GENDecayMode *dmp = new GENDecayMode(gam);
      GENPDTEntry  *d1p = new GENPDTWBoson();
      GENPDTEntry  *d2p = new GENPDTPhoton();
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp);
  }
}

XBoson::~XBoson()
{
   delete fVVStarDecayPtr;
   delete fWBosonPtr;
   delete fZBosonPtr;
   delete fAmpUtilsPtr;
}

//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
Complex_t XBoson::AmpUtils::FullAmplitude(ANL4DVector &k,
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
   Double_t gamw   = v1p->GetWidth();

   Double_t glw    = -kGw*kSqh;
   Double_t grw    = 0;

   Double_t gamz   = v2p->GetWidth();

   Double_t qf     = f3p->GetCharge();
   Double_t t3f    = f3p->GetISpin();
   Double_t glz    = -kGz*(t3f - qf*kSin2W);
   Double_t grz    = -kGz*(    - qf*kSin2W);

   Double_t gzwh   = kGw*kM_w*fXBosonPtr->fFhwz;

   static const Bool_t kIsIncoming = kTRUE;
   static const Bool_t kIsOutgoing = kFALSE;

   HELFermion f1 (p[0], m[0], helFinal[0], +1, kIsOutgoing);
   HELFermion f2b(p[1], m[1], helFinal[1], -1, kIsIncoming);
   HELVector  wf(f2b, f1, glw, grw, kM_w, gamw);

   HELFermion f3 (p[2], m[2], helFinal[2], +1, kIsOutgoing);
   HELFermion f4b(p[3], m[3], helFinal[3], -1, kIsIncoming);
   HELVector  zf (f4b, f3, glz, grz, kM_z, gamz);

   HELScalar  h(k);

   Complex_t amp = HELVertex(wf, zf, h, gzwh);

   return amp;
}

//_____________________________________________________________________________
// --------------------------
//  SelectHelicities
// --------------------------
void XBoson::AmpUtils::SelectHelicities(Double_t  helcombf,
                                        Int_t     helFinal[4],
                                        Double_t &weight)
{
   static const Int_t kNf = 2;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, -1, +1},
                                           {-1, +1, +1, -1}};
   weight = 1.;
   Int_t combf = (Int_t)(helcombf*kNf);
         combf = TMath::Min(combf, kNf-1);
   helFinal  [0] = kFHelComb[combf][0];
   helFinal  [1] = kFHelComb[combf][1];
   helFinal  [2] = kFHelComb[combf][2];
   helFinal  [3] = kFHelComb[combf][3];
   weight *= kNf;
}
