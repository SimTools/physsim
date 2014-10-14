#ifndef EEWWSPRING_H
#define EEWWSPRING_H
//*****************************************************************************
//* =====================
//*  EEWWSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> EEWW generator
//*
//* (Update Record)
//*    2014/10/13  K.Fujii	Original version.
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
//  class EEWWBases
// =====================
//-----------------------------------------------------------------------
class EEWWBases : public JSFBases {
friend class EEWWSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEWWBases(const char *name  = "EEWWBases", 
            const char *title = "EEWW Bases");
  virtual ~EEWWBases();

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

  Double_t  DSigmaDX     (GENBranch &wwbranch);
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoEEWW  (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &emf,
                          const HELFermion &epf,
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
  Double_t fWidth;          // Gam_h  : width of H

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPolem;          // electron polarization
  Double_t fPolep;          // positron polarization
  Int_t    fWmModesLo;      // Wm decay mode lo;
  Int_t    fWmModesHi;      // Wm decay mode hi;
  Int_t    fWpModesLo;      // Wp decay mode lo;
  Int_t    fWpModesHi;      // Wp decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTWBoson *fWmBosonPtr;      //! PD table entry of "Wm"
  GENPDTWBoson *fWpBosonPtr;      //! PD table entry of "Wp"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z0"

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
  Double_t       fQ2WW;           // q^2 of the WW system
  Double_t       fCosWm;          // cos(theta_Wm) in the WW frame
  Double_t       fPhiWm;          // phi_Wm        in the WW frame
  Double_t       fQ2Wm;           // q^2 of the Wm system
  Double_t       fCosWmF;         // cos(theta_fub) in the Wm frame
  Double_t       fPhiWmF;         // phi_fub        in the Wm frame
  Double_t       fQ2Wp;           // q^2 of the Wp system
  Double_t       fCosWpF;         // cos(theta_fu) in the Wp frame
  Double_t       fPhiWpF;         // phi_fu        in the Wp frame
  GENDecayMode  *fWmModePtr;      // pointer to Wm decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Wm daughter (fub)
  GENPDTEntry   *f4Ptr;           // point to 2nd Wm daughter (fd)
  GENDecayMode  *fWpModePtr;      // pointer to Wp decay mode
  GENPDTEntry   *f5Ptr;           // point to 1st Wp daughter (fu)
  GENPDTEntry   *f6Ptr;           // point to 2nd Wp daughter (fdb)
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWmMode;         // Wm decay mode
  Int_t          fWpMode;         // Wp decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [e-, e+, fub, fd, fu, fdb]
  Double_t       fM[6];           // [0,1,2,3,4,5] = [m1, m2, m3, m4 , m5, m6 ]
  Double_t       fSh1;            // sh1
  Double_t       fCh1;            // ch1
  Double_t       fSh2;            // sh2
  Double_t       fCh2;            // ch2

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fWmDecayMode;    // decay mode selector for Wm
  Double_t       fWpDecayMode;    // decay mode selector for Wp
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> e+
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> e+
  Double_t       fXQ2WW;          // q^2 of the WW system
  Double_t       fXQ2Wm;          // q^2 of the Wm system
  Double_t       fXQ2Wp;          // q^2 of the Wp system

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1] = (coswm, phiwm, coswmf, 
  Double_t       fXL[6];          //!          phiwmf, coswpf, phiwpf)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(EEWWBases, 1) // Bases for e+e- -> XX process
};

class EEWWSpring;

//_______________________________________________________________________
// =====================
//  class EEWWSpringBuf
// =====================
//-----------------------------------------------------------------------
class EEWWSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEWWSpringBuf(const char *name   = "EEWWSpringBuf", 
                const char *title  = "EEWW Spring test event buffer",
	        EEWWSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEWWSpringBuf() {}

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
  Double_t fEta1;           // eta_1 for e- -> e-
  Double_t fEta2;           // eta_2 for e+ -> e+
  Double_t fPhi1;           // phi_1 for e- -> e-
  Double_t fPhi2;           // phi_2 - phi_1 for e+ -> e+
  Double_t fQ2WW;           // q^2 of the WW system
  Double_t fCosWm;          // cos(theta_Wm) in the WW frame
  Double_t fPhiWm;          // phi_Wm        in the WW frame
  Double_t fQ2Wm;           // q^2 of the Wm system
  Double_t fCosWmF;         // cos(theta_fub) in the Wm frame
  Double_t fPhiWmF;         // phi_fub        in the Wm frame
  Double_t fQ2Wp;           // q^2 of the Wp system
  Double_t fCosWpF;         // cos(theta_fu) in the Wp frame
  Double_t fPhiWpF;         // phi_fu        in the Wp frame

  ClassDef(EEWWSpringBuf, 1)  // EEWWSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class EEWWSpring
// =====================
//-----------------------------------------------------------------------
class EEWWSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   EEWWSpring(const char *name  = "EEWWSpring", 
              const char *title = "EEWW Spring test",
              EEWWBases  *bases = 0);
   virtual ~EEWWSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(EEWWSpring, 1)  // EEWWSpring class
};
#endif
