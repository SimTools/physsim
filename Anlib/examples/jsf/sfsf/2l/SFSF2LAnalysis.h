#ifndef __SFSF2LANALYSIS__
#define __SFSF2LANALYSIS__
//* $Id$
//*************************************************************************
//* ========================
//*  SFSF2LAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyzes sfsf data to select acoplanar 2-lepton events.
//* (Requires)
//* 	library Anlib
//* 	library SFSFStudy
//* (Provides)
//* 	class SFSF2LAnalysis
//* 	class SFSF2LAnalysisBuf
//* (Usage)
//*   Take a look at anl2L2J.C.
//* (Update Recored)
//*   2002/11/01  K.Ikematsu    Original version.
//*
//*************************************************************************

#include <string.h>
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
//  SFSF2LAnalysisBuf Class
//  -----------------------
//
//  This class is to store data summary of a selected event.
//  Add more data members as needed.
//
class SFSF2LAnalysis;

class SFSF2LAnalysisBuf : public JSFEventBuf 
{
friend class SFSF2LAnalysis;
public:
  SFSF2LAnalysisBuf(const char *name="SFSF2LAnalysisBuf",
		   const char *title="SFSF 2L Data", SFSF2LAnalysis *mod=0);
  SFSF2LAnalysisBuf(SFSF2LAnalysis *mod, const char *name="SFSF2LAnalysisBuf",
		   const char *title="SFSF 2L Data");
  virtual ~SFSF2LAnalysisBuf()
  {
      if (fSpringName) delete fSpringName;
  }

  inline Double_t GetEcm()           const { return fEcm; }
  inline Double_t GetMsf()           const { return fMsf; }
  inline Double_t GetMlsp()          const { return fMlsp; }
  inline Char_t * GetSpringName()    const { return fSpringName; }

  inline void     SetEcm(Double_t rs=350.) { fEcm = rs; }
  inline void     SetSpringName(const Char_t *n="SFSFSpring") 
  {
      if (fSpringName) delete fSpringName;
      fSpringName = new Char_t [strlen(n)]; 
      strcpy(fSpringName,n);
  }

  inline Bool_t   SolveKinematics(Double_t rs, 
                                  Double_t mp, 
                                  Double_t md,
                                  const ANL4DVector &p1, 
                                  const ANL4DVector &p2,
                                  ANL3DVector esol[2]) const 
  {
     ANL3DVector e1(p1.Get3D().Unit());
     ANL3DVector e2(p2.Get3D().Unit());
     Double_t beta  = TMath::Sqrt((1-2*mp/rs)*(1+2*mp/rs));
     Double_t ee1   =  (rs*p1.E() - (mp-md)*(mp+md))/p1.E()/rs/beta;
     Double_t ee2   = -(rs*p2.E() - (mp-md)*(mp+md))/p2.E()/rs/beta;
     Double_t e1e2  = e1*e2;
     Double_t aee2  = (1-e1e2)*(1+e1e2);

     if (aee2 < 0.) return kFALSE;

     aee2 = 1/TMath::Sqrt(aee2);
     ANL3DVector ee[4];
     ee[1] = e1;
     ee[2] = aee2*(e2 - e1e2*e1);
     ee[3] = (ee[1].Cross(ee[2])).Unit();
     ANL3DVector eee[2];
     eee[0](1) = ee1;
     eee[0](2) = aee2*(ee2 - e1e2*ee1);
     eee[0](3) = TMath::Sqrt(TMath::Max(1-eee[0](1)*eee[0](1)
                                         -eee[0](2)*eee[0](2), 0.));
     eee[1](1) =  eee[0](1);
     eee[1](2) =  eee[0](2);
     eee[1](3) = -eee[0](3);

     for (Int_t j=0; j<2; j++) {
        ANL3DVector e;
        for (Int_t i=1; i<4; i++) {
           e += eee[j](i)*ee[i];
        }
        esol[j] = e;
     }
     return kTRUE;
  }

private:
   Char_t  *fSpringName;   //! Name of Spring module
   Double_t fEcm;          // CM energy (GeV)
   Int_t    fNtracks;      // Multiplicity
   Int_t    fNleptons;     // Lepton multiplicity
   Double_t fElepm;        // l^- energy
   Double_t fElepp;        // l^+ energy
   Double_t fEvis;         // Visible energy
   Double_t fPt;           // Pt
   Double_t fPl;           // Pl
   Double_t fCoslm;        // cos(theta_lepton^-)
   Double_t fCoslp;        // cos(theta_lepton^+)
   Double_t fCosSFm1;      // cos(theta_SF^-)_1
   Double_t fCosSFm2;      // cos(theta_SF^-)_2
   Double_t fCosSFmG;      // cos(theta_SF^-)_G
   Double_t fMll;          // m(l^+l^-)
   Double_t fAcop;         // Acoplanarity
   Double_t fMsf;	   // Msf
   Double_t fMlsp;	   // Mlsp


  ClassDef(SFSF2LAnalysisBuf, 1) // SFSF2LAnalysis Buffer Example
};

//_____________________________________________________________________
//  --------------------
//  SFSF2LAnalysis Class
//  --------------------
//
//
class SFSF2LAnalysis : public JSFModule 
{
private:
  Int_t    xNtracks;      // Number of tracks
  Int_t    xNleptons;     // Number of leptons
  Double_t xElLo;         // Lepton E minimum
  Double_t xElHi;         // Lepton E maximum
  Double_t xEvisLo;       // Minimum visible energy
  Double_t xEvisHi;       // Minimum visible energy
  Double_t xPt;           // Pt minimum
  Double_t xPl;           // Pl minimum
  Double_t xCosl;         // |cos(theta_l)| maximum
  Double_t xCoslw;        // -Q_l*cos(theta_l) maximum
  Double_t xMll;          // |m_ll-m_Z| minimum
  Double_t xAcop;         // Acoplanarity minimum

  static Int_t Ngoods;	// Number of good events
  Char_t cutName[MAXCUT][256]; // Cut names
public:
  TCanvas *cHist;	//!
  TH1F *hStat;		//!
  TH1F *hNtracks;	//!
  TH1F *hEvis;		//!
  TH1F *hPt;		//!
  TH1F *hNleptons;	//!
  TH1F *hElepton;	//!
  TH1F *hCoslm;		//!
  TH1F *hCoslp;		//!
  TH1F *hMll;		//!
  TH1F *hAcop;		//!

  TH1F *hMllFinal;      //!
  TH1F *hElmFinal;      //!
  TH1F *hElpFinal;      //!
  TH1F *hElFinal;       //!
  TH1F *hCoslmFinal;    //!
  TH1F *hCoslpFinal;    //!
  TH1F *hPtFinal;       //!
  TH1F *hAcopFinal;     //!
  TH1F *hCosSFmFinal;   //!
  TH2F *hEvisPlFinal;   //! 
  TH2F *hCosSFm1G;      //! (sol1,gen)
  TH2F *hCosSFm2G;      //! (sol2,gen)
  TH2F *hCosSFm12;      //! (sol1,sol2)
  TH2F *hCosSFmGB;      //! (good,bad)

public:
  SFSF2LAnalysis() : JSFModule("SFSF2LAnalysis", "SFSF2LAnalysis Example") {}
  SFSF2LAnalysis(const Char_t *name, const Char_t *title);
  virtual ~SFSF2LAnalysis();

  inline void SetNtracksCut(Int_t    x) { xNtracks  = x; }
  inline void SetNleptonCut(Int_t    x) { xNleptons = x; }
  inline void SetEleptonCut(Double_t x) { xElLo     = x; }
  inline void SetElepLoCut (Double_t x) { xElLo     = x; }
  inline void SetElepHiCut (Double_t x) { xElHi     = x; }
  inline void SetEvisCut   (Double_t x) { xEvisLo   = x; }
  inline void SetEvisLoCut (Double_t x) { xEvisLo   = x; }
  inline void SetEvisHiCut (Double_t x) { xEvisHi   = x; }
  inline void SetPtCut     (Double_t x) { xPt       = x; }
  inline void SetPlCut     (Double_t x) { xPl       = x; }
  inline void SetCosLepCut (Double_t x) { xCosl     = x; }
  inline void SetCosLepwCut(Double_t x) { xCoslw    = x; }
  inline void SetMllCut    (Double_t x) { xMll      = x; }
  inline void SetAcopCut   (Double_t x) { xAcop     = x; }

  inline void SetSpringName(const Char_t *n="SFSFSpring") 
  {
     ((SFSF2LAnalysisBuf *)EventBuf())->SetSpringName(n);
  }

  Bool_t Initialize();
  Bool_t Process(Int_t ev);
  Bool_t Terminate();
  void DrawHist();

  ClassDef(SFSF2LAnalysis, 1) // SFSF2LAnalysis Example
};

#endif
