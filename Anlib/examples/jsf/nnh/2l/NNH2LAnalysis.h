#ifndef NNH2LAnalysis_H
#define NNH2LAnalysis_H

//**********************************************************************
//* ======================
//*  NNH2LAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> NNH -> 2L.
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library NNHStudy (also in physsim libraries)
//* (Provides)
//*	class NNH2LAnalysis
//* (Usage)
//*	...
//* (Update Record)
//*   2009/04/28 K.Fujii	Original version.
//**********************************************************************

#include "JSFModule.h"
#include "TMath.h"

// ----------------------
//  NNH2LAnalysis Class
// ----------------------

class NNH2LAnalysis : public JSFModule {
public:
  NNH2LAnalysis() : JSFModule("NNH2LAnalysis","NNH2LAnalysis Example") {}
  NNH2LAnalysis(const Char_t *name, const Char_t *title);
  virtual ~NNH2LAnalysis();


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
  inline void SetCosLepCut (Double_t x)	{ fCosLepCut  = x; }
  inline void SetM2lCut    (Double_t x) { fM2lCut     = x; }

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
  Double_t fCosLepCut;  // cos(theta_lepton)
  Double_t fM2lCut;	// |m_jj - m_w| maximum

  ClassDef(NNH2LAnalysis, 1)	// NNH2LAnalysis Example
};
#endif
