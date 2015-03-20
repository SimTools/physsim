#ifndef NNHSPRING_H
#define NNHSPRING_H
//*****************************************************************************
//* =====================
//*  NNHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> NNH generator
//*
//* (Update Record)
//*    2008/12/29  K.Fujii	Original version.
//*    2015/02/21  K.Fujii	Implemented anomalous HWW couplings.
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
//  class NNHBases
// =====================
//-----------------------------------------------------------------------
class NNHBases : public JSFBases {
friend class NNHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNHBases(const char *name  = "NNHBases", 
           const char *title = "NNH Bases");
  virtual ~NNHBases();

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
  void     SetPole     (Double_t p      ) { fPole      = p;    }

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
  Complex_t AmpEEtoNNH   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &ne,
                          const HELFermion &neb,
                          const HELScalar  &hs);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fLambda;         // Lambda  : scale of anomalous couplings
  Double_t fA;              // a
  Double_t fB;              // b
  Double_t fBtilde;         // b~

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fXi;             // xi
  Double_t       fEta1;           // eta_1 for e- -> ne
  Double_t       fEta2;           // eta_2 for e+ -> neb
  Double_t       fPhi1;           // phi_1 for e- -> ne
  Double_t       fPhi2;           // phi_2 for e+ -> neb
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [3];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [h , ne, neb]
  Double_t       fM[3];           // [0,1,2] = [mh, mn, mn]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> ne
  Double_t       fXEta2;          // eta_2 for e+ -> neb
  Double_t       fXPhi1;          // phi_1 for e- -> ne
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> neb

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(NNHBases, 1) // Bases for e+e- -> XX process
};

class NNHSpring;

//_______________________________________________________________________
// =====================
//  class NNHSpringBuf
// =====================
//-----------------------------------------------------------------------
class NNHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNHSpringBuf(const char *name   = "NNHSpringBuf", 
               const char *title  = "NNH Spring test event buffer",
	        NNHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NNHSpringBuf() {}

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
  Double_t fEta1;           // eta_1 for e- -> ne
  Double_t fEta2;           // eta_2 for e+ -> neb
  Double_t fPhi1;           // phi_1 for e- -> ne
  Double_t fPhi2;           // phi_2 - phi_1 for e+ -> neb

  ClassDef(NNHSpringBuf, 1)  // NNHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class NNHSpring
// =====================
//-----------------------------------------------------------------------
class NNHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   NNHSpring(const char *name  = "NNHSpring", 
             const char *title = "NNH Spring test",
              NNHBases  *bases = 0);
   virtual ~NNHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(NNHSpring, 1)  // NNHSpring class
};
#endif
