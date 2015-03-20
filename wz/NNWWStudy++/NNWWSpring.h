#ifndef NNWWSPRING_H
#define NNWWSPRING_H
//*****************************************************************************
//* =====================
//*  NNWWSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> NNWW generator
//*
//* (Update Record)
//*    2014/09/19  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class HBoson;

//_______________________________________________________________________
// =====================
//  class NNWWBases
// =====================
//-----------------------------------------------------------------------
class NNWWBases : public JSFBases {
friend class NNWWSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNWWBases(const char *name  = "NNWWBases", 
            const char *title = "NNWW Bases");
  virtual ~NNWWBases();

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
  Double_t GetQ2WW     ()           const { return fQ2WW;      }
  Double_t GetCosWm    ()           const { return fCosWm;     }
  Double_t GetPhiWm    ()           const { return fPhiWm;     }
  Double_t GetQ2Wm     ()           const { return fQ2Wm;      }
  Double_t GetCosWmF   ()           const { return fCosWmF;    }
  Double_t GetPhiWmF   ()           const { return fPhiWmF;    }
  Double_t GetQ2Wp     ()           const { return fQ2Wp;      }
  Double_t GetCosWpF   ()           const { return fCosWpF;    }
  Double_t GetPhiWpF   ()           const { return fPhiWpF;    }
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

  Double_t  DSigmaDX     (GENBranch &wwbranch);
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoNNWW   (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELVector  &wm,
                           const HELVector  &wp);

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
  Int_t    fWmModesLo;      // W- decay mode lo;
  Int_t    fWmModesHi;      // W- decay mode hi;
  Int_t    fWpModesLo;      // W+ decay mode lo;
  Int_t    fWpModesHi;      // W+ decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

  // ----------------
  //  Particle Data
  // ----------------
  HBoson       *fHBosonPtr;       //! PD table entry of "H0"
  GENPDTWBoson *fWmBosonPtr;      //! PD table entry of "W-"
  GENPDTWBoson *fWpBosonPtr;      //! PD table entry of "W+"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

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
  Double_t       fQ2WW;           // q^2 of the WW system
  Double_t       fCosWm;          // cos(theta_W-) in the WW frame
  Double_t       fPhiWm;          // phi_W-        in the WW frame
  Double_t       fQ2Wm;           // q^2 of the W- system
  Double_t       fCosWmF;         // cos(theta_fub) in the W- frame
  Double_t       fPhiWmF;         // phi_fub        in the W- frame
  Double_t       fQ2Wp;           // q^2 of the W+ system
  Double_t       fCosWpF;         // cos(theta_fu) in the W+ frame
  Double_t       fPhiWpF;         // phi_fu        in the W+ frame
  GENDecayMode  *fWmModePtr;      // pointer to W- decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st W- daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd W- daughter
  GENDecayMode  *fWpModePtr;      // pointer to W+ decay mode
  GENPDTEntry   *f5Ptr;           // point to 1st W+ daughter
  GENPDTEntry   *f6Ptr;           // point to 2nd W+ daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWmMode;         // W- decay mode
  Int_t          fWpMode;         // W+ decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [ne, neb, fub, fd, fu, fdb]
  Double_t       fM[6];           // [0,1,2,3,4,5] = [m1, m2 , m3 , m4, m5, m6 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fWmDecayMode;    // decay mode selector for W-
  Double_t       fWpDecayMode;    // decay mode selector for W+
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> ne
  Double_t       fXEta2;          // eta_2 for e+ -> neb
  Double_t       fXPhi1;          // phi_1 for e- -> ne
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> neb
  Double_t       fXQ2WW;          // q^2 of the WW system
  Double_t       fXCosWm;         // cos(theta_W-) in the WW frame
  Double_t       fXPhiWm;         // phi_W-        in the WW frame
  Double_t       fXQ2Wm;          // q^2 of the W- system
  Double_t       fXCosWmF;        // cos(theta_fub) in the WW frame
  Double_t       fXPhiWmF;        // phi_fub        in the WW frame
  Double_t       fXQ2Wp;          // q^2 of the W+ system
  Double_t       fXCosWpF;        // cos(theta_fu) in the WW frame
  Double_t       fXPhiWpF;        // phi_fu        in the WW frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(NNWWBases, 1) // Bases for e+e- -> XX process
};

class NNWWSpring;

//_______________________________________________________________________
// =====================
//  class NNWWSpringBuf
// =====================
//-----------------------------------------------------------------------
class NNWWSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNWWSpringBuf(const char *name   = "NNWWSpringBuf", 
                const char *title  = "NNWW Spring test event buffer",
	        NNWWSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NNWWSpringBuf() {}

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
  Double_t GetQ2WW     ()           const { return fQ2WW;      }
  Double_t GetCosWm    ()           const { return fCosWm;     }
  Double_t GetPhiWm    ()           const { return fPhiWm;     }
  Double_t GetQ2Wm     ()           const { return fQ2Wm;      }
  Double_t GetCosWmF   ()           const { return fCosWmF;    }
  Double_t GetPhiWmF   ()           const { return fPhiWmF;    }
  Double_t GetQ2Wp     ()           const { return fQ2Wp;      }
  Double_t GetCosWpF   ()           const { return fCosWpF;    }
  Double_t GetPhiWpF   ()           const { return fPhiWpF;    }

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
  Double_t fQ2WW;           // q^2 of the WW system
  Double_t fCosWm;          // cos(theta_W-) in the WW frame
  Double_t fPhiWm;          // phi_W-        in the WW frame
  Double_t fQ2Wm;           // q^2 of the W- system
  Double_t fCosWmF;         // cos(theta_fub) in the W- frame
  Double_t fPhiWmF;         // phi_fub        in the W- frame
  Double_t fQ2Wp;           // q^2 of the W+ system
  Double_t fCosWpF;         // cos(theta_fu) in the W+ frame
  Double_t fPhiWpF;         // phi_fu        in the W+ frame

  ClassDef(NNWWSpringBuf, 1)  // NNWWSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class NNWWSpring
// =====================
//-----------------------------------------------------------------------
class NNWWSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   NNWWSpring(const char *name  = "NNWWSpring", 
              const char *title = "NNWW Spring test",
              NNWWBases  *bases = 0);
   virtual ~NNWWSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(NNWWSpring, 1)  // NNWWSpring class
};
#endif
