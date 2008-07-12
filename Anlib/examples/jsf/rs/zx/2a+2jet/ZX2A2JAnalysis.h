#ifndef ZX2A2JANALYSIS_H
#define ZX2A2JANALYSIS_H

//**********************************************************************
//* ========================
//*  ZX2A2JAnalysis Classes
//* ========================
//*
//* (Description)
//*   Analysis class for e+e- -> ZX -> 4 jets.
//* (Requires)
//*   library Anlib     (in LEDA and JSF)
//*   library RSZXStudy (in physsim libraries)
//* (Provides)
//*   class ZX2A2JAnalysis
//*   class ZX2A2JAnalysisBuf
//* (Usage)
//*     ...
//* (Update Record)
//*   2007/02/03 K.Fujii         Original version.
//**********************************************************************

#include "JSFModule.h"

//______________________________________________________________________
// ======================= 
//  ZX2A2JAnalysis Class
// ======================= 

class ZX2A2JAnalysis : public JSFModule {
public:
  ZX2A2JAnalysis() : JSFModule("ZX2A2JAnalysis","ZX2A2JAnalysis Example") {}
  ZX2A2JAnalysis(const Char_t *name,
                 const Char_t *title);
  virtual ~ZX2A2JAnalysis();
  
  Bool_t Process   (Int_t evt);
  Bool_t Terminate ();

  Double_t GetEtrackCut() const     { return fEtrackCut; }
  Double_t GetEgammaCut() const     { return fEgammaCut; }
  Double_t GetCosCone  () const     { return fCosCone;   }
  Double_t GetEconeCut () const     { return fEconeCut;  }
  Double_t GetYcutCut  () const     { return fYcutCut;   }
  Double_t GetM2aCut   () const     { return fM2aCut;    }

  void     SetEtrackCut(Double_t x) { fEtrackCut = x;    }
  void     SetEgammaCut(Double_t x) { fEgammaCut = x;    }
  void     SetCosCone  (Double_t x) { fCosCone   = x;    }
  void     SetEconeCut (Double_t x) { fEconeCut  = x;    }
  void     SetYcutCut  (Double_t x) { fYcutCut   = x;    }
  void     SetM2aCut   (Double_t x) { fM2aCut    = x;    }

private:
  Double_t fEtrackCut;    // cut on track energy
  Double_t fEgammaCut;    // cut on energy of isolated gammas
  Double_t fCosCone;      // cone angle
  Double_t fEconeCut;     // cut on cone energy for isolated gammas
  Double_t fYcutCut;      // cut on Ycut
  Double_t fM2aCut;       // cut on Invariant mass for X

  ClassDef(ZX2A2JAnalysis, 2)	// ZX2A2JAnalysis Example
};
#endif
