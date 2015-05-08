#ifndef EEHSPRING_H
#define EEHSPRING_H
//*****************************************************************************
//* =====================
//*  EEHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> EEH generator
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

//_______________________________________________________________________
// =====================
//  class EEHBases
// =====================
//-----------------------------------------------------------------------
class EEHBases : public JSFBases {
friend class EEHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEHBases(const char *name  = "EEHBases", 
           const char *title = "EEH Bases");
  virtual ~EEHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetLambda   ()           const { return fLambda;    }
  Double_t GetA        ()           const { return fA;         }
  Double_t GetB        ()           const { return fB;         }
  Double_t GetBtilde   ()           const { return fBtilde;    }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass      = m;    }
  void     SetLambda   (Double_t l      ) { fLambda    = l;    }
  void     SetA        (Double_t a      ) { fA         = a;    }
  void     SetB        (Double_t b      ) { fB         = b;    }
  void     SetBtilde   (Double_t bt     ) { fBtilde    = bt;   }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPolem    (Double_t p      ) { fPolem     = p;    }
  void     SetPolep    (Double_t p      ) { fPolep     = p;    }

  // ----------------------
  //   Base class methods
  // ----------------------
  virtual void     Userin();   // Bases user initialization
  virtual void     Userout();  // Bases user output 

  Double_t Func();     // Bases integration function.

  // ----------------------
  //   Utility methods
  // ----------------------
private:
  void     SelectHelicities(Double_t &weight);

  Double_t  DSigmaDX     ();
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoEEH   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &e,
                          const HELFermion &eb,
                          const HELScalar  &hs);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fLambda;         // Lambda : scale of anomalous couplings
  Double_t fA;              // a
  Double_t fB;              // b
  Double_t fBtilde;         // b~

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPolem;          // electron polarization
  Double_t fPolep;          // positron polarization

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fXi;             // xi
  Double_t       fEta1;           // eta_1 for e- -> e-
  Double_t       fEta2;           // eta_2 for e+ -> e+
  Double_t       fPhi1;           // phi_1 for e- -> e-
  Double_t       fPhi2;           // phi_2 for e+ -> e+
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [3];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [h , e-, e+]
  Double_t       fM[3];           // [0,1,2] = [mh, mn, mn]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> e+
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> e+

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(EEHBases, 1) // Bases for e+e- -> XX process
};

class EEHSpring;

//_______________________________________________________________________
// =====================
//  class EEHSpringBuf
// =====================
//-----------------------------------------------------------------------
class EEHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEHSpringBuf(const char *name   = "EEHSpringBuf", 
               const char *title  = "EEH Spring test event buffer",
	        EEHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }

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
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung
  Double_t fXi;             // xi
  Double_t fEta1;           // eta_1 for e- -> e-
  Double_t fEta2;           // eta_2 for e+ -> e+
  Double_t fPhi1;           // phi_1 for e- -> e-
  Double_t fPhi2;           // phi_2 - phi_1 for e+ -> e+

  ClassDef(EEHSpringBuf, 1)  // EEHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class EEHSpring
// =====================
//-----------------------------------------------------------------------
class EEHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   EEHSpring(const char *name  = "EEHSpring", 
             const char *title = "EEH Spring test",
              EEHBases  *bases = 0);
   virtual ~EEHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(EEHSpring, 1)  // EEHSpring class
};
#endif
