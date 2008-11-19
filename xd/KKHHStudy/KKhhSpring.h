#ifndef KKHHSPRING_H
#define KKHHSPRING_H
//*****************************************************************************
//* =====================
//*  KKhhSpring
//* =====================
//*  
//* (Description)
//*    Kaluza-Klein  ee -> hh generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2003/XX/XX  Nicolas Delerue      Original version.
//*    2005/03/01  K.Fujii              Cleaned up and fixed some bugs with 
//*                                     higgs momenta and Lorentz boost.
//*
//*****************************************************************************

#include "TNamed.h"
#include "TH1.h"

#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
const Double_t   kGeV2fb = 0.389379292e12;   // GeV to fb
const Double_t   kAlpha  = 0.0078125;        // alpha(mz)  = 1/128.
const Double_t   kAlpha0 = 1./137.0359895;   // alpha(q=0) = 1/137.
const Double_t   kM_e    = 0.510998902e-3;   // electron mass [GeV]

//_______________________________________________________________________
// =====================
//  class KKhhBases
// =====================
//-----------------------------------------------------------------------
class KKhhBases : public JSFBases {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  KKhhBases(const char *name  = "KKhhBases", 
            const char *title = "KKhh Bases");
  virtual ~KKhhBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
#if 0
  Int_t    GetID     ()           const { return fID;      }
  Double_t GetCharge ()           const { return fCharge;  }
#endif
  Double_t GetMass   ()           const { return fMass;    }
  Double_t GetEcmInit()           const { return fEcmInit; }

  Double_t GetCosth  ()           const { return fCosth;   }
  Double_t GetPhi    ()           const { return fPhi;     }
  Double_t GetZboost ()           const { return fZBoost;  }
  Double_t GetEcmIP  ()           const { return fEcmIP;   }

#if 0
  void     SetID     (Int_t id)         { fID      = id;   }
  void     SetCharge (Double_t c)       { fCharge  = c;    }
#endif
  void     SetMass   (Double_t m)       { fMass    = m;    }
  void     SetEcmInit(Double_t ecm)     { fEcmInit = ecm;  }
  void     SetISR    (Bool_t b = kTRUE) { fISR     = b;    }
  void     SetBeamStr(Bool_t b = kTRUE) { fBeamStr = b;    }

  // ----------------------
  //   Base class methods
  // ----------------------
  virtual void     Userin();   // Bases user initialization
  virtual void     Userout();  // Bases user output 

  Double_t Func();     // Bases integration function.

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fLambda;         // lambda
  Double_t fMScale;         // M_s Mass scale

#if 0
  Int_t    fID;             // Parton ID
  Double_t fCharge;         // Parton charge		 
#endif
  Double_t fMass;           // m_h mass of the higgs 

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?

  // ----------------
  //  Event info
  // ----------------
  Double_t fCosth;          // cos(theta_h) in cm  frame
  Double_t fPhi;            // phi_h        in cm  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t fR_BW_p;         //!                     for e+
  Double_t fR_BS_m;         //! beamstrahlung       for e-
  Double_t fR_BS_p;         //!                     for e+

  Double_t fR_ISR_var;      //! initial state radiation
  Double_t fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t fXU[2];          //! [0] = cos(theta)
  Double_t fXL[2];          //! [1] = phi

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(KKhhBases, 1) // Bases for e+e- -> XX process
};

class KKhhSpring;

//_______________________________________________________________________
// =====================
//  class KKhhSpringBuf
// =====================
//-----------------------------------------------------------------------
class KKhhSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  KKhhSpringBuf(const char *name   = "KKhhSpringBuf", 
                const char *title  = "KKhh Spring test event buffer",
	        KKhhSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~KKhhSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosth  ()           const { return fCosth;   }
  Double_t GetPhi    ()           const { return fPhi;     }
  Double_t GetZboost ()           const { return fZBoost;  }
  Double_t GetEcmIP  ()           const { return fEcmIP;   }

  // ----------------------
  //   Base class methods
  // ----------------------
  virtual Bool_t SetPartons();

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Event info
  // ----------------
  Double_t fCosth;          // cos(theta_h) in cm  frame
  Double_t fPhi;            // phi_h        in cm  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung


  ClassDef(KKhhSpringBuf, 1)  // KKhhSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class KKhhSpring
// =====================
//-----------------------------------------------------------------------
class KKhhSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   KKhhSpring(const char *name  = "KKhhSpring", 
              const char *title = "KKhh Spring test",
              KKhhBases  *bases = 0);
   virtual ~KKhhSpring();

   ClassDef(KKhhSpring, 1)  // KKhhSpring class
};
#endif
