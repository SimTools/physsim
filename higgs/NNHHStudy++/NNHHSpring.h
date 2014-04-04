#ifndef NNHHSPRING_H
#define NNHHSPRING_H
//*****************************************************************************
//* =====================
//*  NNHHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> NNHH generator
//*
//* (Update Record)
//*    2010/04/22  K.Fujii	Original version.
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
//  class NNHHBases
// =====================
//-----------------------------------------------------------------------
class NNHHBases : public JSFBases {
friend class NNHHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNHHBases(const char *name  = "NNHHBases", 
            const char *title = "NNHH Bases");
  virtual ~NNHHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetWidth    ()           const { return fWidth;     }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2HH     ()           const { return fQ2HH;      }
  Double_t GetCosH     ()           const { return fCosH;      }
  Double_t GetPhiH     ()           const { return fPhiH;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass      = m;    }
  void     SetWidth    (Double_t g      ) { fWidth     = g;    }
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
  Complex_t AmpEEtoNNHH   (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELScalar  &h1,
                           const HELScalar  &h2);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fWidth;          // Gam_h  : width of H

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
  Double_t       fQ2HH;           // q^2 of the HH system
  Double_t       fCosH;           // cos(theta_H) in the HH frame
  Double_t       fPhiH;           // phi_h        in the HH frame
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [ne, neb, h1, h2]
  Double_t       fM[4];           // [0,1,2,3] = [mn, mn , mh, mh]

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
  Double_t       fXQ2HH;          // q^2 of the HH system
  Double_t       fXCosH;          // cos(theta_H) in the HH frame
  Double_t       fXPhiH;          // phi_h        in the HH frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(NNHHBases, 1) // Bases for e+e- -> XX process
};

class NNHHSpring;

//_______________________________________________________________________
// =====================
//  class NNHHSpringBuf
// =====================
//-----------------------------------------------------------------------
class NNHHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNHHSpringBuf(const char *name   = "NNHHSpringBuf", 
                const char *title  = "NNHH Spring test event buffer",
	        NNHHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NNHHSpringBuf() {}

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
  Double_t GetQ2HH     ()           const { return fQ2HH;      }
  Double_t GetCosH     ()           const { return fCosH;      }
  Double_t GetPhiH     ()           const { return fPhiH;      }

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
  Double_t fQ2HH;           // q^2 of the HH system
  Double_t fCosH;           // cos(theta_H) in the HH frame
  Double_t fPhiH;           // phi_h        in the HH frame

  ClassDef(NNHHSpringBuf, 1)  // NNHHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class NNHHSpring
// =====================
//-----------------------------------------------------------------------
class NNHHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   NNHHSpring(const char *name  = "NNHHSpring", 
              const char *title = "NNHH Spring test",
              NNHHBases  *bases = 0);
   virtual ~NNHHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(NNHHSpring, 1)  // NNHHSpring class
};
#endif
