#ifndef __ZH4JAnalysis__
#define __ZH4JAnalysis__

//**********************************************************************
//* ======================
//*  ZH4JAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> ZH -> 4 jets.
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library ZHStudy (also in physsim libraries)
//* (Provides)
//*	class ZH4JAnalysis
//*	class ZH4JAnalysisBuf
//* (Usage)
//*	...
//* (Update Record)
//*	18 Nov 1999	A.L.C.Sanchez	Written based on K. Fujii's 
//*					sample physsim-99a-1 classes
//*					(still largely to be corrected!)
//*	22 Nov 1999	A.L.C.Sanchez	Additional modifications to suit ZH.
//**********************************************************************

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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ---------------------
//  ZH4JAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ZH4JAnalysis;

class ZH4JAnalysisBuf : public JSFEventBuf{
friend class ZH4JAnalysis;
private:
  Int_t		fNtracks;	// track multiplicity
  Double_t	fEvis;		// visible energy (GeV)
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

public:
  ZH4JAnalysisBuf(const char *name="ZH4JAnalysisBuf",
                  const char *title="ZH 4-Jet Data", ZH4JAnalysis *mod=0);
  ZH4JAnalysisBuf(ZH4JAnalysis *mod, const char *name="ZH4JAnalysisBuf",
                  const char *title="ZH 4-Jet Data");
  virtual ~ZH4JAnalysisBuf();

  inline Int_t    GetNtracks()	const {return fNtracks; }
  inline Double_t GetEvis()	const {return fEvis; }
  inline Double_t GetPt()	const {return fPt; }
  inline Double_t GetPl()	const {return fPl; }
  inline Int_t    GetNjets()	const {return fNjets; }
  inline Double_t GetYcut()	const {return fYcut; }

  ClassDef(ZH4JAnalysisBuf, 2) // ZH4JAnalysis Buffer Example
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  ZH4JAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ZH4JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;	// No. of tracks
  Double_t xEtrack;	// track energy
  Double_t xEvis;	// minimum visible energy
  Double_t xPt;		// Pt minimum
  Double_t xPl;		// Pl maximum
  Double_t xYcut;	// y_cut to force the event to 4 jets
  Int_t    xNjets;	// No. of Jets
  Double_t xEjet;	// E_jet minimum
  Double_t xCosjet;	// |cos(theta_j)| maximum
  Double_t xCoszh;	// |cos(theta_z)| and |cos(theta_h)| maximum
  Double_t xM2j;	// |m_jj - m_w| maximum
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
  TH2F *hEvisPl;
  TH2F *hEzEh;
  TH2F *hCoszCosh;
  TH1F *hAcop;
  TH1F *hMassh;
  TH1F *hMassz;

public:
  ZH4JAnalysis() : JSFModule("ZH4JAnalysis","ZH4JAnalysis Example") {}
  ZH4JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~ZH4JAnalysis();
  
  void CleanUp(TObjArray *objs);

  inline void SetNtrackCut (Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut (Double_t x) { xEtrack = x; }
  inline void SetEvisCut   (Double_t x) { xEvis = x; }
  inline void SetPtCut     (Double_t x) { xPt = x; }
  inline void SetPlCut     (Double_t x) { xPl = x; }
  inline void SetMinYcut   (Double_t x) { xYcut = x; }
  inline void SetNjetCut   (Int_t    x) { xNjets = x; }
  inline void SetEjetCut   (Double_t x) { xEjet = x; }
  inline void SetCosjetCut (Double_t x) { xCosjet = x; }
  inline void SetCoszhCut  (Double_t x)	{ xCoszh = x; }
  inline void SetM2jCut    (Double_t x) { xM2j = x; }
  inline void SetAcopCut   (Double_t x) { xAcop = x; }

  Bool_t Initialize();
  Bool_t Process(Int_t evt);
  Bool_t Terminate();
  void DrawHist();

  ClassDef(ZH4JAnalysis, 2)	// ZH4JAnalysis Example
};

#endif
