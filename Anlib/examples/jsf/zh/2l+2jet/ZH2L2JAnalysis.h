#ifndef __ZH2L2JAnalysis__
#define __ZH2L2JAnalysis__

//**********************************************************************
//* ======================
//*  ZH2L2JAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> ZH -> 4 jets.
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library ZHStudy (also in physsim libraries)
//* (Provides)
//*	class ZH2L2JAnalysis
//*	class ZH2L2JAnalysisBuf
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
//  ZH2L2JAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ZH2L2JAnalysis;

class ZH2L2JAnalysisBuf : public JSFEventBuf{
friend class ZH2L2JAnalysis;
private:
  Int_t		fNtracks;	// track multiplicity
  Int_t         fNlptracks;     // isolated lepton multiplicity
  Double_t	fEvis;		// visible energy (GeV)
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

public:
  ZH2L2JAnalysisBuf(const char *name="ZH2L2JAnalysisBuf",
                  const char *title="ZH 4-Jet Data", ZH2L2JAnalysis *mod=0);
  ZH2L2JAnalysisBuf(ZH2L2JAnalysis *mod, const char *name="ZH2L2JAnalysisBuf",
                  const char *title="ZH 4-Jet Data");
  virtual ~ZH2L2JAnalysisBuf();

  inline Int_t    GetNtracks()	const {return fNtracks; }
  inline Double_t GetEvis()	const {return fEvis; }
  inline Double_t GetPt()	const {return fPt; }
  inline Double_t GetPl()	const {return fPl; }
  inline Int_t    GetNjets()	const {return fNjets; }
  inline Double_t GetYcut()	const {return fYcut; }

  ClassDef(ZH2L2JAnalysisBuf, 2) // ZH2L2JAnalysis Buffer Example
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  ZH2L2JAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ZH2L2JAnalysis : public JSFModule {
private:
  Int_t    xNtracks;	// No. of tracks
  Double_t xEtrack;	// track energy
  Double_t xEvis;	// minimum visible energy
  Double_t xPt;		// Pt minimum
  Double_t xPl;		// Pl maximum
  Double_t xElepton;    // Elepton mimimum
  Double_t xCosCone;    // cos(theta_cone)
  Double_t xEcone;      // Ecome maximum
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
  TH1F *hNlptracks;
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
  ZH2L2JAnalysis() : JSFModule("ZH2L2JAnalysis","ZH2L2JAnalysis Example") {}
  ZH2L2JAnalysis(const Char_t *name, const Char_t *title);
  virtual ~ZH2L2JAnalysis();
  
  void CleanUp(TObjArray *objs);

  inline void SetNtrackCut (Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut (Double_t x) { xEtrack = x; }
  inline void SetEvisCut   (Double_t x) { xEvis = x; }
  inline void SetPtCut     (Double_t x) { xPt = x; }
  inline void SetPlCut     (Double_t x) { xPl = x; }
  inline void SetElpCut    (Double_t x) { xElepton = x; }
  inline void SetConeAngle (Double_t x) 
                          { xCosCone = TMath::Cos(TMath::Pi()*x/180.); }
  inline void SetEconeCut  (Double_t x) { xEcone = x; }
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

  ClassDef(ZH2L2JAnalysis, 2)	// ZH2L2JAnalysis Example
};

#endif
