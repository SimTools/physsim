#ifndef XN1XN12SL2JANALYSIS_H
#define XN1XN12SL2JANALYSIS_H

//**********************************************************************
//* ============================
//*  XN1XN12SL2JAnalysis Classes
//* ============================
//*
//* (Description)
//*    Analysis class for e+e- -> XN1 XN1 -> (stau tau) (stau tau)
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library XN1XN1Study (also in physsim libraries)
//* (Provides)
//*	class XN1XN12SL2JAnalysis
//* (Usage)
//*	...
//* (Update Record)
//*   2010/10/10  K.Fujii	Original version
//**********************************************************************

#include "JSFModule.h"

class XN1XN12SL2JAnalysis : public JSFModule {
public:
  XN1XN12SL2JAnalysis() : JSFModule("XN1XN12SL2JAnalysis","XN1XN12SL2JAnalysis Example") {}
  XN1XN12SL2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~XN1XN12SL2JAnalysis();

  inline void SetEcm       (Double_t x) { fEcm        = x; }

  inline void SetNtrackCut (Int_t    x) { fNtracksCut = x; }
  inline void SetEtrackCut (Double_t x) { fEtrackCut  = x; }
  inline void SetEvisLoCut (Double_t x) { fEvisLoCut  = x; }
  inline void SetEvisHiCut (Double_t x) { fEvisHiCut  = x; }
  inline void SetPtCut     (Double_t x) { fPtCut      = x; }
  inline void SetPlCut     (Double_t x) { fPlCut      = x; }
  inline void SetEleptonCut(Double_t x) { fEleptonCut = x; }
  inline void SetCosConeCut(Double_t x) { fCosConeCut = x; }
  inline void SetEcone     (Double_t x) { fEconeCut   = x; }
  inline void SetMinYcut   (Double_t x) { fYcutCut    = x; }
  inline void SetNjetCut   (Int_t    x) { fNjetsCut   = x; }
  inline void SetEjetCut   (Double_t x) { fEjetCut    = x; }
  inline void SetCosjetCut (Double_t x) { fCosjetCut  = x; }
  inline void SetM2jCut    (Double_t x) { fM2jCut     = x; }
  inline void SetCosxCut   (Double_t x) { fCosxCut    = x; }

  inline void SetMM1Cut    (Double_t x) { fMM1Cut     = x; }
  inline void SetMM2Cut    (Double_t x) { fMM2Cut     = x; }
  inline void SetAcopCut   (Double_t x) { fAcopCut    = x; }

  inline Double_t GetEcm() const { return fEcm; }

  Bool_t Initialize();
  Bool_t Process(Int_t evt);
  Bool_t Terminate();

private:
  Double_t fEcm;           // Ecm
  Int_t    fNtracksCut;    // No. of tracks
  Double_t fEtrackCut;     // track energy
  Double_t fEvisLoCut;     // Minimum visible energy
  Double_t fEvisHiCut;     // Maximum visible energy
  Double_t fPtCut;         // Pt minimum
  Double_t fPlCut;         // Pl maximum
  Double_t fEleptonCut;    // El minimum
  Double_t fCosConeCut;    // cos(cone) cut
  Double_t fEconeCut;      // Econe cut
  Double_t fYcutCut;       // y_cut to force the event to 4 jets
  Int_t    fNjetsCut;      // No. of jets
  Double_t fEjetCut;       // E_jet minimum
  Double_t fCosjetCut;     // |cos(theta_j)| maximum
  Double_t fM2jCut;        // |m_jj-m_x| maximum
  Double_t fCosxCut;       // |cos(theta_x)| maximum
  Double_t fMM1Cut;        // mm_WW < fMM1Cut or
  Double_t fMM2Cut;        // mm_WW > fMM2Cut
  Double_t fAcopCut;       // Minimum Acoplanarity

  ClassDef(XN1XN12SL2JAnalysis, 2)	// XN1XN12SL2JAnalysis Example
};

#endif
