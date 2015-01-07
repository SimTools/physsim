#ifndef NNZHSPRING_H
#define NNZHSPRING_H
//*****************************************************************************
//* =====================
//*  NNZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> NNZH generator
//*
//* (Update Record)
//*    2015/01/07  K.Fujii	Original version.
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
//  class NNZHBases
// =====================
//-----------------------------------------------------------------------
class NNZHBases : public JSFBases {
friend class NNZHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNZHBases(const char *name  = "NNZHBases", 
            const char *title = "NNZH Bases");
  virtual ~NNZHBases();

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
  Double_t GetQ2ZH     ()           const { return fQ2ZH;      }
  Double_t GetCosZ     ()           const { return fCosZ;      }
  Double_t GetPhiZ     ()           const { return fPhiZ;      }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosZF    ()           const { return fCosZF;     }
  Double_t GetPhiZF    ()           const { return fPhiZF;     }
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

  Double_t  DSigmaDX     (GENBranch &wwbranch);
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoNNZH   (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELVector  &z,
                           const HELScalar  &h);
  Complex_t AmpEEtoNNZ    (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELVector  &z);
  Complex_t AmpEEtoNNH    (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELScalar  &h);

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
  Int_t    fZModesLo;       // Z decay mode lo;
  Int_t    fZModesHi;       // Z decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTZBoson *fZBosonPtr;      //! PD table entry of "Z1"
  GENPDTWBoson *fWBosonPtr;      //! PD table entry of "W"

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
  Double_t       fQ2ZH;           // q^2 of the ZH system
  Double_t       fCosZ;           // cos(theta_Z) in the ZH frame
  Double_t       fPhiZ;           // phi_Z        in the ZH frame
  Double_t       fQ2Z;            // q^2 of the Z system
  Double_t       fCosZF;          // cos(theta_F) in the Z frame
  Double_t       fPhiZF;          // phi_F        in the Z frame
  GENDecayMode  *fZModePtr ;      // pointer to Z decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [5];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[5];           // [0,1,2,3,4,5] = [ne, neb, f3, f4b, h]
  Double_t       fM[5];           // [0,1,2,3,4,5] = [m1, m2 , m3, m4 , m5]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z1
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> ne
  Double_t       fXEta2;          // eta_2 for e+ -> neb
  Double_t       fXPhi1;          // phi_1 for e- -> ne
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> neb
  Double_t       fXQ2ZH;          // q^2 of the Zh system
  Double_t       fXCosZ;          // cos(theta_Z) in the Zh frame
  Double_t       fXPhiZ;          // phi_Z        in the ZH frame
  Double_t       fXQ2Z;           // q^2 of the Z system
  Double_t       fXCosZF;         // cos(theta_Z1) in the Z frame
  Double_t       fXPhiZF;         // phi_Z1        in the Z frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(NNZHBases, 1) // Bases for e+e- -> XX process
};

class NNZHSpring;

//_______________________________________________________________________
// =====================
//  class NNZHSpringBuf
// =====================
//-----------------------------------------------------------------------
class NNZHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NNZHSpringBuf(const char *name   = "NNZHSpringBuf", 
                const char *title  = "NNZH Spring test event buffer",
                NNZHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NNZHSpringBuf() {}

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
  Double_t GetQ2ZH     ()           const { return fQ2ZH;      }
  Double_t GetCosZ     ()           const { return fCosZ;      }
  Double_t GetPhiZ     ()           const { return fPhiZ;      }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosZF    ()           const { return fCosZF;     }
  Double_t GetPhiZF    ()           const { return fPhiZF;     }

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
  Double_t fQ2ZH;           // q^2 of the ZZ system
  Double_t fCosZ;           // cos(theta_Z) in the ZH frame
  Double_t fPhiZ;           // phi_Z        in the ZH frame
  Double_t fQ2Z;            // q^2 of the Z system
  Double_t fCosZF;          // cos(theta_F) in the Z frame
  Double_t fPhiZF;          // phi_F        in the Z frame

  ClassDef(NNZHSpringBuf, 1)  // NNZHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class NNZHSpring
// =====================
//-----------------------------------------------------------------------
class NNZHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   NNZHSpring(const char *name  = "NNZHSpring", 
             const char *title = "NNZH Spring test",
             NNZHBases   *bases = 0);
   virtual ~NNZHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(NNZHSpring, 1)  // NNZHSpring class
};
#endif
