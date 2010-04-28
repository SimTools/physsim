#ifndef NNH2JAnalysis_H
#define NNH2JAnalysis_H

//**********************************************************************
//* ======================
//*  NNH2JAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> NNH -> 2J.
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library NNHStudy (also in physsim libraries)
//* (Provides)
//*	class NNH2JAnalysis
//* (Usage)
//*	...
//* (Update Record)
//*   2009/04/28 K.Fujii	Original version.
//**********************************************************************

#include "JSFModule.h"
#include "TMath.h"

// ----------------------
//  NNH2JAnalysis Class
// ----------------------

class NNH2JAnalysis : public JSFModule {
public:
  NNH2JAnalysis() : JSFModule("NNH2JAnalysis","NNH2JAnalysis Example") {}
  NNH2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~NNH2JAnalysis();


  inline void SetEcm       (Double_t x) { fEcm        = x; }
  inline void SetNtrackCut (Int_t    x) { fNtracksCut = x; }
  inline void SetEtrackCut (Double_t x) { fEtrackCut  = x; }
  inline void SetEvisLoCut (Double_t x) { fEvisLoCut  = x; }
  inline void SetEvisHiCut (Double_t x) { fEvisHiCut  = x; }
  inline void SetPtCut     (Double_t x) { fPtCut      = x; }
  inline void SetPlCut     (Double_t x) { fPlCut      = x; }
  inline void SetEleptonCut(Double_t x) { fEleptonCut = x; }
  inline void SetConeAngle (Double_t x) 
           { fCosConeCut = TMath::Cos(TMath::Pi()*x/180.); }
  inline void SetEconeCut  (Double_t x) { fEconeCut   = x; }
  inline void SetNleptonCut(Int_t    x) { fNleptonCut = x; }
  inline void SetMinYcut   (Double_t x) { fYcutCut    = x; }
  inline void SetNjetCut   (Int_t    x) { fNjetsCut   = x; }
  inline void SetCosjetCut (Double_t x)	{ fCosjetCut  = x; }
  inline void SetCos2jCut  (Double_t x)	{ fCos2jCut   = x; }
  inline void SetM2jCut    (Double_t x) { fM2jCut     = x; }

  inline Double_t GetEcm() const { return fEcm; }

  Bool_t Initialize();
  Bool_t Process(Int_t evt);
  Bool_t Terminate();

private:
  Double_t fEcm;        // Ecm
  Int_t    fNtracksCut;	// No. of tracks
  Double_t fEtrackCut;	// track energy
  Double_t fEvisLoCut;	// minimum visible energy
  Double_t fEvisHiCut;	// maximum visible energy
  Double_t fPtCut;	// Pt minimum
  Double_t fPlCut;	// Pl maximum
  Double_t fEleptonCut; // Elepton minimum
  Double_t fCosConeCut; // cos(theta_cone)
  Double_t fEconeCut;   // Ecome maximum
  Int_t    fNleptonCut;	// maximum No. of isolated leptons
  Double_t fYcutCut;    // y_cut to force the event to 2 jets
  Int_t    fNjetsCut;   // No. of jets
  Double_t fEjetCut;    // E_jet minimum
  Double_t fCosjetCut;  // |cos(theta_jet)| maximum
  Double_t fCos2jCut;   // |cos(theta_qq)| maximum
  Double_t fM2jCut;	// |m_jj - m_h| maximum

  ClassDef(NNH2JAnalysis, 1)	// NNH2JAnalysis Example
};
#endif
