#ifndef __USERANALYSIS__
#define __USERANALYSIS__
//*************************************************************************
//* ========================
//*  XCXC4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC chargino pair data. 
//* (Requires)
//* 	library Anlib
//* 	library XCXCStudy
//* (Provides)
//* 	class XCXC4JAnalysis
//* 	class XCXC4JAnalysisBuf
//* (Usage)
//*   Take a look at anl.C.  
//* (Update Recored)
//*   1999/08/01  K.Fujii	Original version.
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
//  -----------------------
//  XCXC4JAnalysisBuf Class
//  -----------------------
//
//
class XCXC4JAnalysis;

class XCXC4JAnalysisBuf : public JSFEventBuf {
friend class XCXC4JAnalysis;
private:
  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity

public:
  XCXC4JAnalysisBuf(const char *name="XCXC4JAnalysisBuf",
                 const char *title="XCXC 4-Jet Data", XCXC4JAnalysis *mod=0);
  XCXC4JAnalysisBuf(XCXC4JAnalysis *mod, const char *name="XCXC4JAnalysisBuf",
                 const char *title="XCXC 4-Jet Data");
  virtual ~XCXC4JAnalysisBuf() {}
  
  inline Int_t    GetNtracks()	const { return fNtracks; }
  inline Double_t GetEvis()	const { return fEvis; }
  inline Double_t GetPt()	const { return fPt; }
  inline Double_t GetPl()	const { return fPl; }
  inline Int_t    GetNjets()	const { return fNjets; }
  inline Double_t GetYcut()	const { return fYcut; }

  ClassDef(XCXC4JAnalysisBuf, 1) // XCXC4JAnalysis Buffer Example
};

//_____________________________________________________________________
//  --------------------
//  XCXC4JAnalysis Class
//  --------------------
//
//
class XCXC4JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;    // No. of tracks
  Double_t xEtrack;     // track energy
  Double_t xEvisLo;     // Minimum visible energy
  Double_t xEvisHi;     // Maximum visible energy
  Double_t xPt;         // Pt minimum
  Double_t xPl;         // Pl maximum
  Double_t xEl;         // Elepton maximum
  Double_t xYcut;       // y_cut to force the event to 4 jets
  Int_t    xNjets;      // No. of jets
  Double_t xEjet;	// E_jet minimum
  Double_t xCosjet1;	// |cos(theta_j)| maximum
  Double_t xCosjet2;	// |cos(theta_j)| maximum
  Double_t xCosw;	// |cos(theta_w)| maximum
  Double_t xM2jLo;	// m_W - xM2jLo < m_jj
  Double_t xM2jHi;	// m_jj < m_W + xM2jHi
  Double_t xMM1;        // mm_WW < xMM1 or
  Double_t xMM2;        // mm_WW > xMM2
  Double_t xAcop;	// Acoplanarity minimum

  static Int_t Ngoods;	// Number of good events
  Char_t cutName[MAXCUT][100]; // Cut names
public:
  TCanvas *cHist;
  TH1F *hStat;
  TH1F *hNtracks;
  TH1F *hEvis;
  TH1F *hPt;
  TH1F *hNjets;
  TH1F *hEjet;
  TH1F *hElep;
  TH1F *hCosjet;
  TH1F *hNsols;
  TH1F *hChi2;
  TH2F *hEw1Ew2;
  TH2F *hCosw1Cosw2;
  TH2F *hMw1Mw2;
  TH2F *hEvisPl;
  TH1F *hMM;
  TH1F *hEw;
  TH1F *hAcop;
public:
  XCXC4JAnalysis() : JSFModule("XCXC4JAnalysis", "XCXC4JAnalysis Example") {}
  XCXC4JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~XCXC4JAnalysis();

  inline void SetNtrackCut(Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut(Double_t x) { xEtrack = x; }
  inline void SetEvisCut  (Double_t x) { xEvisLo = x; }
  inline void SetEvisLoCut(Double_t x) { xEvisLo = x; }
  inline void SetEvisHiCut(Double_t x) { xEvisHi = x; }
  inline void SetPtCut    (Double_t x) { xPt = x; }
  inline void SetPlCut    (Double_t x) { xPl = x; }
  inline void SetElCut    (Double_t x) { xEl = x; }
  inline void SetMinYcut  (Double_t x) { xYcut  = x; }
  inline void SetNjetCut  (Int_t    x) { xNjets = x; }
  inline void SetEjetCut  (Double_t x) { xEjet = x; }
  inline void SetCosjetCut(Double_t x) { xCosjet2 = x; }
  inline void SetCosjetCut(Double_t x1, Double_t x2) 
                                       { xCosjet1 = x1; xCosjet2 = x2; }
  inline void SetCoswCut  (Double_t x) { xCosw = x; }
  inline void SetM2jCut   (Double_t x) { xM2jLo = x; }
  inline void SetM2jLoCut (Double_t x) { xM2jLo = x; }
  inline void SetM2jHiCut (Double_t x) { xM2jHi = x; }
  inline void SetMM1Cut   (Double_t x) { xMM1 = x; }
  inline void SetMM2Cut   (Double_t x) { xMM2 = x; }
  inline void SetAcopCut  (Double_t x) { xAcop = x; }

  void CleanUp(TObjArray *objs);
  
  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();   
  void DrawHist();

  ClassDef(XCXC4JAnalysis, 1) // XCXC4JAnalysis Example
};

#endif

