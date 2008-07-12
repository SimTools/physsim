#ifndef ZX4JANALYSIS_H
#define ZX4JANALYSIS_H

//**********************************************************************
//* ========================
//*  ZX4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   Analysis class for e+e- -> ZX -> 4 jets.
//* (Requires)
//*   library Anlib     (in LEDA and JSF)
//*   library RSZXStudy (in physsim libraries)
//* (Provides)
//*   class ZX4JAnalysis
//*   class ZX4JAnalysisBuf
//* (Usage)
//*     ...
//* (Update Record)
//*   2007/02/22 K.Fujii         Original version.
//**********************************************************************

#include "JSFModule.h"

//______________________________________________________________________
// ======================= 
//  ZX4JAnalysis Class
// ======================= 

class ZX4JAnalysis : public JSFModule {
public:
  ZX4JAnalysis() : JSFModule("ZX4JAnalysis","ZX4JAnalysis Example") {}
  ZX4JAnalysis(const Char_t *name,
                 const Char_t *title);
  virtual ~ZX4JAnalysis();
  
  void CleanUp(TObjArray *objs);

  Bool_t Process   (Int_t evt);
  Bool_t Terminate ();

  Double_t GetEtrackCut() const     { return fEtrackCut; }
  Double_t GetYcutCut  () const     { return fYcutCut;   }
  Double_t GetM2jCut   () const     { return fM2jCut;    }

  void     SetEtrackCut(Double_t x) { fEtrackCut = x;    }
  void     SetYcutCut  (Double_t x) { fYcutCut   = x;    }
  void     SetM2jCut   (Double_t x) { fM2jCut    = x;    }

private:
  Double_t fEtrackCut;    // cut on track energy
  Double_t fYcutCut;      // cut on Ycut
  Double_t fM2jCut;       // cut on Invariant mass for X

  ClassDef(ZX4JAnalysis, 2)	// ZX4JAnalysis Example
};
#endif
