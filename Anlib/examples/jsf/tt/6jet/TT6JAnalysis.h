#ifndef __TT6JAnalysis__
#define __TT6JAnalysis__
//*************************************************************************
//* ======================
//*  TT6JAnalysis Classes
//* ======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC chargino pair data. 
//* (Requires)
//* 	library Anlib
//* 	library TTStudy
//* (Provides)
//* 	class TT6JAnalysis
//* 	class TT6JAnalysisBuf
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
//  TT6JAnalysisBuf Class
//  ---------------------
//  
//  This class is to store data summary of a selected event.
//  Add more data members as needed.
//
class TT6JAnalysis;

class TT6JAnalysisBuf : public JSFEventBuf {
friend class TT6JAnalysis;
public:
  TT6JAnalysisBuf(const char *name="TT6JAnalysisBuf",
                  const char *title="TT 6-Jet Data", TT6JAnalysis *mod=0);
  TT6JAnalysisBuf(TT6JAnalysis *mod, const char *name="TT6JAnalysisBuf",
                  const char *title="TT 6-Jet Data");
  virtual ~TT6JAnalysisBuf() {}
  
  inline Int_t    GetNtracks()	const { return fNtracks; }
  inline Double_t GetEvis()	const { return fEvis; }
  inline Double_t GetPt()	const { return fPt; }
  inline Double_t GetPl()	const { return fPl; }
  inline Int_t    GetNjets()	const { return fNjets; }
  inline Double_t GetYcut()	const { return fYcut; }
  inline Double_t GetThrust()	const { return fThrust; }

private:
  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t  	fThrust;	// thrust

  ClassDef(TT6JAnalysisBuf, 1) // TT6JAnalysis Buffer Example
};

//_____________________________________________________________________
//  ------------------
//  TT6JAnalysis Class
//  ------------------
//
//
class TT6JAnalysis : public JSFModule {
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
  Double_t xCosbw;	// cos(theta_bw) maximum
  Double_t xM2j;	// |m_jj-m_W| maximum
  Double_t xM3j;	// |m_3j-m_t| maximum
  Double_t xThrust;	// Thrust maximum

  static Int_t Ngoods;	// Number of good events
  Char_t cutName[MAXCUT][256]; // Cut names
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
  TH2F *hCosbw1Cosbw2;
  TH2F *hMw1Mw2;
  TH2F *hMt1Mt2;
  TH2F *hEvisPl;
  TH1F *hThrust;
  TH1F *hYcut;
public:
  TT6JAnalysis() : JSFModule("TT6JAnalysis", "TT6JAnalysis Example") {}
  TT6JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~TT6JAnalysis();

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
  inline void SetCosbwCut (Double_t x) { xCosbw = x; }
  inline void SetM2jCut   (Double_t x) { xM2j = x; }
  inline void SetM3jCut   (Double_t x) { xM3j = x; }
  inline void SetThrustCut(Double_t x) { xThrust = x; }
  
  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();   
  void DrawHist();

  ClassDef(TT6JAnalysis, 1) // TT6JAnalysis Example
};

#endif



