#ifndef AXANALYSIS_H
#define AXANALYSIS_H

//**********************************************************************
//* ======================
//*  AXAnalysis Classes
//* ======================
//*
//* (Description)
//*    Analysis class for e+e- -> gamma X -> 3 gammas
//* (Requires)
//*	library Anlib (from physsim of K. Fujii)
//*	library AXStudy (also in physsim libraries)
//* (Provides)
//*	class AXAnalysis
//* (Usage)
//*	...
//* (Update Record)
//*    2007/01/17  K.Fujii	Original version.
//**********************************************************************

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "JSFModule.h"
#include "JSFSIMDST.h"
#include "Anlib.h"

#define MAXCUT 50

// -------------------------------------------------------------------
//  ------------------
//  AXAnalysis Class
//  ------------------

class AXAnalysis : public JSFModule {
public:
  AXAnalysis() : JSFModule("AXAnalysis","AXAnalysis Example") {}
  AXAnalysis(const Char_t *name, const Char_t *title);
  virtual ~AXAnalysis();
  
  inline void SetNtrackCut (Int_t    x) { xNtracks = x; }
  inline void SetEtrackCut (Double_t x) { xEtrack = x;  }
  inline void SetEvisCut   (Double_t x) { xEvis = x;    }
  inline void SetPtCut     (Double_t x) { xPt = x;      }
  inline void SetPlCut     (Double_t x) { xPl = x;      }
  inline void SetMinYcut   (Double_t x) { xYcut = x;    }
  inline void SetNjetCut   (Int_t    x) { xNjets = x;   }
  inline void SetEjetCut   (Double_t x) { xEjet = x;    }
  inline void SetCosjetCut (Double_t x) { xCosjet = x;  }
  inline void SetCosrsaxCut(Double_t x) { xCosrsax = x; }
  inline void SetM2jCut    (Double_t x) { xM2j = x;     }

  Bool_t Initialize();
  Bool_t Process(Int_t evt);
  Bool_t Terminate();

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
  Double_t xCosrsax;	// |cos(theta_z)| and |cos(theta_h)| maximum
  Double_t xM2j;	// |m_jj - m_w| maximum

  ClassDef(AXAnalysis, 2)	// AXAnalysis Example
};

#endif
