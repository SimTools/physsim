#ifndef TT2L2JAnalysis_H
#define TT2L2JAnalysis_H
//*************************************************************************
//* ========================
//*  TT2L2JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ttbar data to select 2-lepton + 2-jet events.
//* (Requires)
//* 	library Anlib
//* 	library TTStudy
//* (Provides)
//* 	class TT2L2JAnalysis
//* 	class TT2L2JAnalysisBuf
//* (Usage)
//*   Take a look at anl2L2J.C.
//* (Update Recored)
//*   1999/08/19  K.Ikematsu    Derived from TTL4JAnalysis.h.
//*   2001/07/07  K.Ikematsu    Modified for MacOS X.
//*
//* $Id$
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
//  TT2L2JAnalysisBuf Class
//  -----------------------
//
//  This class is to store data summary of a selected event.
//  Add more data members as needed.
//
class TT2L2JAnalysis;

class TT2L2JAnalysisBuf : public JSFEventBuf {
friend class TT2L2JAnalysis;
public:
  TT2L2JAnalysisBuf(const char *name="TT2L2JAnalysisBuf",
		   const char *title="TT Lepton+4-Jet Data", TT2L2JAnalysis *mod=0);
  TT2L2JAnalysisBuf(TT2L2JAnalysis *mod, const char *name="TT2L2JAnalysisBuf",
		   const char *title="TT Lepton+4-Jet Data");
  virtual ~TT2L2JAnalysisBuf() {}

  inline Int_t    GetNtracks()	const { return fNtracks; }
  inline Double_t GetEvis()	const { return fEvis; }
  inline Double_t GetPt()	const { return fPt; }
  inline Double_t GetPl()	const { return fPl; }
  inline Int_t    GetNjets()	const { return fNjets; }
  inline Double_t GetYcut()	const { return fYcut; }
  inline Double_t GetThrust()	const { return fThrust; }
  inline Double_t GetEcm()      const { Double_t ecm = 349.4; return ecm; }

private:
  Int_t     	fNtracks;	// track multiplicity
  Int_t         fNlptracks;     // lepton multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t      ecm;            // Center of Energy
  Double_t  	fYcut;		// y_cut to force the event to 2 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t  	fThrust;	// thrust

  ClassDef(TT2L2JAnalysisBuf, 1) // TT2L2JAnalysis Buffer Example
};

//_____________________________________________________________________
//  --------------------
//  TT2L2JAnalysis Class
//  --------------------
//
//
class TT2L2JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;    // No. of tracks
  Double_t xEtrack;     // track energy
  Double_t xEvis;       // Minimum visible energy
  Double_t xPt;         // Pt maximum
  Double_t xPl;         // Pl maximum
  Double_t xElepton;    // Elepton mimimum
  Double_t xCosCone;    // cos(theta_cone)
  Double_t xEcone;      // Ecome maximum
  Double_t xYcut;       // y_cut to force the event to 2 jets
  Int_t    xNjets;      // No. of jets
  Double_t xEjet;	// E_jet minimum
  Double_t xCosjet;	// |cos(theta_j)| maximum
  Double_t xThrust;	// Thrust maximum

  static Int_t Ngoods;	// Number of good events
  Char_t cutName[MAXCUT][256]; // Cut names
public:
  TCanvas *cHist;
  TH1F *hStat;
  TH1F *hNtracks;
  TH1F *hEvis;
  TH1F *hPt;
  TH1F *hNlptracks;
  TH1F *hNlpmtracks;
  TH1F *hNlpptracks;
  TH2F *hElpmlpp;
  TH2F *hCoslpmlpp;
  TH1F *hNjets;
  TH1F *hEjet;
  TH1F *hCosjet;
  TH2F *hEvisPl;
  TH1F *hThrust;
  TH1F *hYcut;
public:
  TT2L2JAnalysis() : JSFModule("TT2L2JAnalysis", "TT2L2JAnalysis Example") {}
  TT2L2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~TT2L2JAnalysis();

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
  inline void SetThrustCut(Double_t x) { xThrust = x; }

  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();
  void DrawHist();

  ClassDef(TT2L2JAnalysis, 1) // TT2L2JAnalysis Example
};

#endif
