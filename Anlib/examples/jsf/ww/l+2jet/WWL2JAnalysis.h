#ifndef __WWL2JANALYSIS__
#define __WWL2JANALYSIS__
//*************************************************************************
//* =======================
//*  WWL2JAnalysis Classes
//* =======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC W boson pair data. 
//* (Requires)
//* 	library Anlib
//* 	library WWStudy
//* (Provides)
//* 	class WWL2JAnalysis
//* 	class WWL2JAnalysisBuf
//* (Usage)
//*   Take a look at anlL2J.C.  
//* (Update Recored)
//*   2000/04/28  K.Ikematsu 	Derived from WW4JAnalysis.h.
//*
//*************************************************************************
//
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "JSFSteer.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "Anlib.h"

#define MAXCUT 50
//_____________________________________________________________________
//  ----------------------
//  WWL2JAnalysisBuf Class
//  ----------------------
//
//  This class is to store data summary of a selected event.
//  Add more data members as needed.
//
class WWL2JAnalysis;

class WWL2JAnalysisBuf : public JSFEventBuf {
friend class WWL2JAnalysis;
public:
  WWL2JAnalysisBuf(const char *name="WWL2JAnalysisBuf",
                   const char *title="WW Lepton+2-Jets Data", WWL2JAnalysis *mod=0);
  WWL2JAnalysisBuf(WWL2JAnalysis *mod, const char *name="WWL2JAnalysisBuf",
                   const char *title="WW Lepton+2-Jets Data");
  virtual ~WWL2JAnalysisBuf() {}

  inline Int_t    GetNtracks()   const { return fNtracks; }
  inline Double_t GetEvis()      const { return fEvis; }
  inline Double_t GetPt()        const { return fPt; }
  inline Double_t GetPl()        const { return fPl; }
  inline Int_t    GetNlptracks() const { return fNlptracks; }
  inline Int_t    GetNjets()     const { return fNjets; }
  inline Double_t GetYcut()      const { return fYcut; }
  inline Double_t GetAcop()      const { return fAcop; }
  inline Double_t GetEcm()       const { Double_t ecm = 339.4; return ecm; }

private:
  Int_t         fNtracks;        // track multiplicity
  Double_t      fEvis;           // visible energy
  Double_t      fPt;             // Pt
  Double_t      fPl;             // Pl
  Int_t         fNlptracks;      // isolated lepton multiplicity
  Double_t      fYcut;           // y_cut to force the event to 2 jets
  Int_t         fNjets;          // jet multiplicity
  Double_t      fAcop;           // Acoplanarity

  ClassDef(WWL2JAnalysisBuf, 1)  // WWL2JAnalysis Buffer Example
};

//_____________________________________________________________________
//  -------------------
//  WWL2JAnalysis Class
//  -------------------
//
//
class WWL2JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;    // No. of tracks
  Double_t xEtrack;     // track energy
  Double_t xEvis;       // Maximum visible energy
  Double_t xPt;         // Pt maximum
  Double_t xPl;         // Pl maximum
  Double_t xElepton;    // Elepton mimimum
  Double_t xCosCone;    // cos(theta_cone)
  Double_t xEcone;      // Ecome maximum
  Double_t xYcut;       // y_cut to force the event to 2 jets
  Int_t    xNjets;      // No. of jets
  Double_t xEjet;	// E_jet minimum
  Double_t xCosjet;	// |cos(theta_j)| maximum
  Double_t xCosw;	// |cos(theta_w)| maximum
  Double_t xM2j;	// |m_jj-m_W| maximum
  Double_t xAcop;	// Acoplanarity maximum

  static Int_t Ngoods;	// Number of good events
  Char_t cutName[MAXCUT][100]; // Cut names
public:
  TCanvas *cHist;
  TH1F *hStat;
  TH1F *hNtracks;
  TH1F *hEvis;
  TH1F *hPt;
  TH1F *hNlptracks;
  TH1F *hYcut;
  TH1F *hNjets;
  TH1F *hEjet;
  TH1F *hCosjet;
  TH2F *hEw1Ew2;
  TH2F *hCosw1Cosw2;
  TH2F *hMw1Mw2;
  TH2F *hEvisPl;
  TH1F *hAcop;
public:
  WWL2JAnalysis() : JSFModule("WWL2JAnalysis", "WWL2JAnalysis Example") {}
  WWL2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~WWL2JAnalysis();

  void CleanUp(TObjArray *objs);

  inline void SetNtrackCut(Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut(Double_t x) { xEtrack = x; }
  inline void SetEvisCut  (Double_t x) { xEvis = x; }
  inline void SetPtCut    (Double_t x) { xPt = x; }
  inline void SetPlCut    (Double_t x) { xPl = x; }
  inline void SetElpCut   (Double_t x) { xElepton = x; }
  inline void SetConeAngle(Double_t x) 
                          { xCosCone = TMath::Cos(TMath::Pi()*x/180.); }
  inline void SetEconeCut (Double_t x) { xEcone = x; }
  inline void SetMinYcut  (Double_t x) { xYcut  = x; }
  inline void SetNjetCut  (Int_t    x) { xNjets = x; }
  inline void SetEjetCut  (Double_t x) { xEjet = x; }
  inline void SetCosjetCut(Double_t x) { xCosjet = x; }
  inline void SetCoswCut  (Double_t x) { xCosw = x; }
  inline void SetM2jCut   (Double_t x) { xM2j = x; }
  inline void SetAcopCut  (Double_t x) { xAcop = x; }
  
  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();   
  void DrawHist();

  ClassDef(WWL2JAnalysis, 1) // WWL2JAnalysis Example
};

#endif
