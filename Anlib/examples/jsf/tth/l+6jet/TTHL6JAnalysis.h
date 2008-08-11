#ifndef TTHL6JANALYSIS_H
#define TTHL6JANALYSIS_H
//*************************************************************************
//* =======================
//*  TTHL6JAnalysis Classes
//* =======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ttbar data to select lepton + 6-jet events.
//* (Requires)
//* 	library Anlib
//* 	library TTHStudy
//* (Provides)
//* 	class TTHL6JAnalysis
//* 	class TTHL6JAnalysisBuf
//* (Usage)
//*   Take a look at anlL6J.C.
//* (Update Recored)
//*   2002/08/15  K.Fujii	Derived from TTL4JAnalysis.h.
//*   2008/08/08  K.Fujii	Clean up.
//*
//*************************************************************************

#include "JSFModule.h"
#include "TMath.h"

//_____________________________________________________________________
//  -------------------
//  TTHL6JAnalysis Class
//  -------------------
//
//
class TTHL6JAnalysis : public JSFModule {
public:
  TTHL6JAnalysis() : JSFModule("TTHL6JAnalysis", "TTHL6JAnalysis Example") {}
  TTHL6JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~TTHL6JAnalysis();

  inline void SetEcm       (Double_t x) { fEcm        = x; }

  inline void SetNtrackCut (Int_t    x) { fNtracksCut = x; }
  inline void SetEtrackCut (Double_t x) { fEtrackCut  = x; }
  inline void SetEvisCut   (Double_t x) { fEvisCut    = x; }
  inline void SetPtCut     (Double_t x) { fPtCut      = x; }
  inline void SetPlCut     (Double_t x) { fPlCut      = x; }
  inline void SetElpCut    (Double_t x) { fEleptonCut = x; }
  inline void SetConeAngle (Double_t x)
                           { fCosConeCut = TMath::Cos(TMath::Pi()*x/180.); }
  inline void SetEconeCut  (Double_t x) { fEconeCut   = x; }
  inline void SetMinYcut   (Double_t x) { fYcutCut    = x; }
  inline void SetNjetCut   (Int_t    x) { fNjetsCut   = x; }
  inline void SetEjetCut   (Double_t x) { fEjetCut    = x; }
  inline void SetCosjetCut (Double_t x) { fCosjetCut  = x; }
  inline void SetCosbwCut  (Double_t x) { fCosbwCut   = x; }
  inline void SetM2jCut    (Double_t x) { fM2jCut     = x; }
  inline void SetM3jCut    (Double_t x) { fM3jCut     = x; }
  inline void SetThrustCut (Double_t x) { fThrustCut  = x; }
  inline void SetBtagNsig  (Double_t x) { fBtagNsig   = x; }
  inline void SetBtagNoffv (Int_t    x) { fBtagNoffv  = x; }
  inline void SetBTtagNsig (Double_t x) { fBTtagNsig  = x; }
  inline void SetBTtagNoffv(Int_t    x) { fBTtagNoffv = x; }

  inline Double_t GetEcm() const { return fEcm; }

  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();

private:
  Double_t fEcm;           // Ecm
  Int_t    fNtracksCut;    // No. of tracks
  Double_t fEtrackCut;     // track energy
  Double_t fEvisCut;       // Minimum visible energy
  Double_t fPtCut;         // Pt maximum
  Double_t fPlCut;         // Pl maximum
  Double_t fEleptonCut;    // Elepton mimimum
  Double_t fCosConeCut;    // cos(theta_cone)
  Double_t fEconeCut;      // Ecome maximum
  Double_t fYcutCut;       // y_cut to force the event to 4 jets
  Int_t    fNjetsCut;      // No. of jets
  Double_t fEjetCut;       // E_jet minimum
  Double_t fCosjetCut;     // |cos(theta_j)| maximum
  Double_t fCosbwCut;      // cos(theta_bw) maximum
  Double_t fM2jCut;        // |m_jj-m_W| maximum
  Double_t fM3jCut;        // |m_3j-m_t| maximum
  Double_t fThrustCut;     // Thrust maximum
  Double_t fBtagNsig;      // Nsig  for b-tag (loose)
  Int_t    fBtagNoffv;     // Noffv for b-tag (loose)
  Double_t fBTtagNsig;     // Nsig  for b-tag (tight)
  Int_t    fBTtagNoffv;    // Noffv for b-tag (tight)

  ClassDef(TTHL6JAnalysis, 1) // TTHL6JAnalysis Example
};

#endif
