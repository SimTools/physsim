#ifndef __USERANALYSIS__
#define __USERANALYSIS__
//*************************************************************************
//* ======================
//*  WW4JAnalysis Classes
//* ======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC chargino pair data. 
//* (Requires)
//* 	library Anlib
//* 	library WWStudy
//* (Provides)
//* 	class WW4JAnalysis
//* 	class WW4JAnalysisBuf
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
//  ---------------------
//  WW4JAnalysisBuf Class
//  ---------------------
//
//
class WW4JAnalysis;

class WW4JAnalysisBuf : public JSFEventBuf {
friend class WW4JAnalysis;
private:
  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity

public:
  WW4JAnalysisBuf(const char *name="WW4JAnalysisBuf",
                  const char *title="WW 4-Jet Data", WW4JAnalysis *mod=0);
  WW4JAnalysisBuf(WW4JAnalysis *mod, const char *name="WW4JAnalysisBuf",
                  const char *title="WW 4-Jet Data");
  virtual ~WW4JAnalysisBuf() {}
  
  inline Int_t    GetNtracks()	const { return fNtracks; }
  inline Double_t GetEvis()	const { return fEvis; }
  inline Double_t GetPt()	const { return fPt; }
  inline Double_t GetPl()	const { return fPl; }
  inline Int_t    GetNjets()	const { return fNjets; }
  inline Double_t GetYcut()	const { return fYcut; }

  ClassDef(WW4JAnalysisBuf, 1) // WW4JAnalysis Buffer Example
};

//_____________________________________________________________________
//  ------------------
//  WW4JAnalysis Class
//  ------------------
//
//
class WW4JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;    // No. of tracks
  Double_t xEtrack;     // track energy
  Double_t xEvis;       // Minimum visible energy
  Double_t xPt;         // Pt maximum
  Double_t xPl;         // Pl maximum
  Double_t xYcut;       // y_cut to force the event to 4 jets
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
  TH1F *hNjets;
  TH1F *hEjet;
  TH1F *hCosjet;
  TH1F *hNsols;
  TH1F *hChi2;
  TH2F *hEw1Ew2;
  TH2F *hCosw1Cosw2;
  TH2F *hMw1Mw2;
  TH2F *hEvisPl;
  TH1F *hAcop;
public:
  WW4JAnalysis() : JSFModule("WW4JAnalysis", "WW4JAnalysis Example") {}
  WW4JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~WW4JAnalysis();

  void CleanUp(TObjArray *objs);

  inline void SetNtrackCut(Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut(Double_t x) { xEtrack = x; }
  inline void SetEvisCut  (Double_t x) { xEvis = x; }
  inline void SetPtCut    (Double_t x) { xPt = x; }
  inline void SetPlCut    (Double_t x) { xPl = x; }
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

  ClassDef(WW4JAnalysis, 1) // WW4JAnalysis Example
};

#endif

