#ifndef ZHAH2JANALYSIS_H
#define ZHAH2JANALYSIS_H
//*************************************************************************
//* ========================
//*  ZHAH2JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC eta_R eta_I data. 
//* (Requires)
//* 	library Anlib
//* 	library ZHAHStudy
//* (Provides)
//* 	class ZHAH2JAnalysis
//* 	class ZHAH2JAnalysisBuf
//* (Usage)
//*   Take a look at anl2J.C.  
//* (Update Recored)
//*   2008/12/30  K.Fujii	Original version.
//*
//*************************************************************************

#include "JSFModule.h"
#include "TMath.h"

//_____________________________________________________________________
//  --------------------
//  ZHAH2JAnalysis Class
//  --------------------
//
//
class ZHAH2JAnalysis : public JSFModule {
public:
  ZHAH2JAnalysis() : JSFModule("ZHAH2JAnalysis", "ZHAH2JAnalysis") {}
  ZHAH2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~ZHAH2JAnalysis();

  inline void SetEcm       (Double_t x) { fEcm        = x; }

  inline void SetNtrackCut (Int_t    x) { fNtracksCut = x; }
  inline void SetEtrackCut (Double_t x) { fEtrackCut  = x; }
  inline void SetEvisLoCut (Double_t x) { fEvisLoCut  = x; }
  inline void SetEvisHiCut (Double_t x) { fEvisHiCut  = x; }
  inline void SetPtCut     (Double_t x) { fPtCut      = x; }
  inline void SetPlCut     (Double_t x) { fPlCut      = x; }
  inline void SetElCut     (Double_t x) { fElCut      = x; }
  inline void SetMinYcut   (Double_t x) { fYcutCut    = x; }
  inline void SetNjetCut   (Int_t    x) { fNjetsCut   = x; }
  inline void SetEjetCut   (Double_t x) { fEjetCut    = x; }
  inline void SetCosjetCut (Double_t x) { fCosjetCut  = x; }
  inline void SetM2jLoCut  (Double_t x) { fM2jLoCut   = x; }
  inline void SetM2jHiCut  (Double_t x) { fM2jHiCut   = x; }
  inline void SetCoshCut   (Double_t x) { fCoshCut    = x; }
  inline void SetBtagNsig  (Double_t x) { fBtagNsig   = x; }
  inline void SetBtagNoffv (Int_t    x) { fBtagNoffv  = x; }
  inline void SetBTtagNsig (Double_t x) { fBTtagNsig  = x; }
  inline void SetBTtagNoffv(Int_t    x) { fBTtagNoffv = x; }

  inline void SetMM1Cut    (Double_t x) { fMM1Cut     = x; }
  inline void SetMM2Cut    (Double_t x) { fMM2Cut     = x; }

  inline Double_t GetEcm() const { return fEcm; }
  
  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();   

private:
  Double_t fEcm;           // Ecm
  Int_t    fNtracksCut;    // No. of tracks
  Double_t fEtrackCut;     // track energy
  Double_t fEvisLoCut;     // Minimum visible energy
  Double_t fEvisHiCut;     // Maximum visible energy
  Double_t fPtCut;         // Pt minimum
  Double_t fPlCut;         // Pl maximum
  Double_t fElCut;         // El maximum
  Double_t fYcutCut;       // y_cut to force the event to 4 jets
  Int_t    fNjetsCut;      // No. of jets
  Double_t fEjetCut;       // E_jet minimum
  Double_t fCosjetCut;     // |cos(theta_j)| maximum
  Double_t fM2jLoCut;      // m_jj minimum
  Double_t fM2jHiCut;      // m_jj maximum
  Double_t fCoshCut;       // |cos(theta_h)| maximum
  Double_t fMM1Cut;        // mm_WW < fMM1Cut or
  Double_t fMM2Cut;        // mm_WW > fMM2Cut
  Double_t fBtagNsig;      // Nsig  for b-tag (loose)
  Int_t    fBtagNoffv;     // Noffv for b-tag (loose)
  Double_t fBTtagNsig;     // Nsig  for b-tag (tight)
  Int_t    fBTtagNoffv;    // Noffv for b-tag (tight)

  ClassDef(ZHAH2JAnalysis, 1) // ZHAH2JAnalysis Example
};

#endif
