#ifndef __KKHHAnalysis__
#define __KKHHAnalysis__

//**********************************************************************
//* ======================
//*  KKHHAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> hh -> 4 jets.
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//* (Provides)
//*	class KKHHAnalysis
//*	class KKHHAnalysisBuf
//* (Usage)
//*	...
//* (Update Record)
//*	09 Dec 2003	Nicolas Delerue Adapted from ZH analysis
//*
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
//  KKHHAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class KKHHAnalysis;

class KKHHAnalysisBuf : public JSFEventBuf{
friend class KKHHAnalysis;
private:
  Int_t		fNtracks;	// track multiplicity
  Double_t	fEvis;		// visible energy (GeV)
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

public:
  KKHHAnalysisBuf(const char *name="KKHHAnalysisBuf",
                  const char *title="hh 4-Jet Data", KKHHAnalysis *mod=0);
  KKHHAnalysisBuf(KKHHAnalysis *mod, const char *name="KKHHAnalysisBuf",
                  const char *title="hh 4-Jet Data");
  virtual ~KKHHAnalysisBuf();

  inline Int_t    GetNtracks()	const {return fNtracks; }
  inline Double_t GetEvis()	const {return fEvis; }
  inline Double_t GetPt()	const {return fPt; }
  inline Double_t GetPl()	const {return fPl; }
  inline Int_t    GetNjets()	const {return fNjets; }
  inline Double_t GetYcut()	const {return fYcut; }

  ClassDef(KKHHAnalysisBuf, 2) // KKHHAnalysis Buffer Example
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  KKHHAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class KKHHAnalysis : public JSFModule {
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
  Double_t xMassTotDist; //Total distance to the higgs mass

  Double_t crossSection; // The events cross section

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
  TH1F *hAcopZoom;
  TH1F *hMassh1;
  TH1F *hMassh2;
  TH2F *hMasshh;

  TH1F *hMasshTot;

  TH1F *hMassh1raw;
  TH1F *hMassh2raw;
  TH2F *hMasshhraw;
  
  TH1F *hMassh1ac;
  TH1F *hMassh2ac;
  TH2F *hMasshhac;
  
  TH1F *hCosTheta;


public:
  KKHHAnalysis() : JSFModule("KKHHAnalysis","KKHHAnalysis Example") {}
  KKHHAnalysis(const Char_t *name, const Char_t *title);
  virtual ~KKHHAnalysis();
  
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
  inline void SetMassTotDist    (Double_t x) { xMassTotDist = x; }
  inline void SetAcopCut   (Double_t x) { xAcop = x; }
  inline void SetCrossSection   (Double_t x) { crossSection = x; }

  Bool_t Initialize();
  Bool_t Process(Int_t evt);
  Bool_t Terminate();
  void DrawHist();

  ClassDef(KKHHAnalysis, 2)	// KKHHAnalysis Example
};

#endif
