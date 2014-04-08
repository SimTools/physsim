#ifndef EEZSPRING_H
#define EEZSPRING_H
//*****************************************************************************
//* =====================
//*  EEZSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> EEZ generator
//*
//* (Update Record)
//*    2014/04/03  K.Fujii	Original version.
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
//  class EEZBases
// =====================
//-----------------------------------------------------------------------
class EEZBases : public JSFBases {
friend class EEZSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEZBases(const char *name  = "EEZBases", 
           const char *title = "EEZ Bases");
  virtual ~EEZBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosF;      }
  Double_t GetPhiF     ()           const { return fPhiF;      }
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

  Double_t  DSigmaDX     ();
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoEEZ   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &e,
                          const HELFermion &eb,
                          const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPolem;          // electron polarization
  Double_t fPolep;          // positron polarization
  Int_t    fZModesLo;       // Z decay mode lo;
  Int_t    fZModesHi;       // Z decay mode hi;

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
  Double_t       fQ2Z;            // q^2 of the final state Z
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Double_t       fCosF;           // cos(theta_f) in the Z frame
  Double_t       fPhiF;           // phi_f        in the Z frame
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [e-, e+,  f, fb]
  Double_t       fM[4];           // [0,1,2,3] = [me, me, mf, mf]
  Double_t       fSh1;            // sh1
  Double_t       fCh1;            // ch1
  Double_t       fSh2;            // sh2
  Double_t       fCh2;            // ch2

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> e+
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> e+
  Double_t       fXQ2Z;           // q^2 of the final state Z
  Double_t       fXCosF;          // cos(theta_f) in the Z frame
  Double_t       fXPhiF;          // phi_f        in the Z frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[2];          //! [0,1] = (cosf, phif)
  Double_t       fXL[2];          //!

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(EEZBases, 1) // Bases for e+e- -> XX process
};

class EEZSpring;

//_______________________________________________________________________
// =====================
//  class EEZSpringBuf
// =====================
//-----------------------------------------------------------------------
class EEZSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEZSpringBuf(const char *name   = "EEZSpringBuf", 
               const char *title  = "EEZ Spring test event buffer",
	        EEZSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEZSpringBuf() {}

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
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosF;      }
  Double_t GetPhiF     ()           const { return fPhiF;      }

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
  Double_t fQ2Z;            // q^2 of the final state Z
  Double_t fCosF;           // cos(theta_f) in the Z frame
  Double_t fPhiF;           // phi_f        in the Z frame

  ClassDef(EEZSpringBuf, 1)  // EEZSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class EEZSpring
// =====================
//-----------------------------------------------------------------------
class EEZSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   EEZSpring(const char *name  = "EEZSpring", 
             const char *title = "EEZ Spring test",
              EEZBases  *bases = 0);
   virtual ~EEZSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(EEZSpring, 1)  // EEZSpring class
};
#endif
