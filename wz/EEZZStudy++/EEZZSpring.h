#ifndef EEZZSPRING_H
#define EEZZSPRING_H
//*****************************************************************************
//* =====================
//*  EEZZSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> EEZZ generator
//*
//* (Update Record)
//*    2014/09/22  K.Fujii	Original version.
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
//  class EEZZBases
// =====================
//-----------------------------------------------------------------------
class EEZZBases : public JSFBases {
friend class EEZZSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEZZBases(const char *name  = "EEZZBases", 
            const char *title = "EEZZ Bases");
  virtual ~EEZZBases();

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
  Double_t GetQ2ZZ     ()           const { return fQ2ZZ;      }
  Double_t GetCosZ1    ()           const { return fCosZ1;     }
  Double_t GetPhiZ1    ()           const { return fPhiZ1;     }
  Double_t GetQ2Z1     ()           const { return fQ2Z1;      }
  Double_t GetCosZ1F   ()           const { return fCosZ1F;    }
  Double_t GetPhiZ1F   ()           const { return fPhiZ1F;    }
  Double_t GetQ2Z2     ()           const { return fQ2Z2;      }
  Double_t GetCosZ2F   ()           const { return fCosZ2F;    }
  Double_t GetPhiZ2F   ()           const { return fPhiZ2F;    }
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

  Double_t  DSigmaDX     (GENBranch &zzbranch);
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoEEZZ  (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &emf,
                          const HELFermion &epf,
                          const HELVector  &z1,
                          const HELVector  &z2);

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
  Int_t    fZ1ModesLo;      // Z1 decay mode lo;
  Int_t    fZ1ModesHi;      // Z1 decay mode hi;
  Int_t    fZ2ModesLo;      // Z2 decay mode lo;
  Int_t    fZ2ModesHi;      // Z2 decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTZBoson *fZ1BosonPtr;      //! PD table entry of "Z1"
  GENPDTZBoson *fZ2BosonPtr;      //! PD table entry of "Z2"

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
  Double_t       fQ2ZZ;           // q^2 of the ZZ system
  Double_t       fCosZ1;          // cos(theta_Z1) in the ZZ frame
  Double_t       fPhiZ1;          // phi_Z1        in the ZZ frame
  Double_t       fQ2Z1;           // q^2 of the Z1 system
  Double_t       fCosZ1F;         // cos(theta_fub) in the Z1 frame
  Double_t       fPhiZ1F;         // phi_fub        in the Z1 frame
  Double_t       fQ2Z2;           // q^2 of the Z2 system
  Double_t       fCosZ2F;         // cos(theta_fu) in the Z2 frame
  Double_t       fPhiZ2F;         // phi_fu        in the Z2 frame
  GENDecayMode  *fZ1ModePtr;      // pointer to Z1 decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z1 daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z1 daughter
  GENDecayMode  *fZ2ModePtr;      // pointer to Z2 decay mode
  GENPDTEntry   *f5Ptr;           // point to 1st Z2 daughter
  GENPDTEntry   *f6Ptr;           // point to 2nd Z2 daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZ1Mode;         // Z1 decay mode
  Int_t          fZ2Mode;         // Z2 decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [e-, e+, f3, f4b, f5, f6b]
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
  Double_t       fZ1DecayMode;    // decay mode selector for Z1
  Double_t       fZ2DecayMode;    // decay mode selector for Z2
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> e+
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> e+
  Double_t       fXQ2ZZ;          // q^2 of the ZZ system
  Double_t       fXQ2Z1;          // q^2 of the Z1 system
  Double_t       fXQ2Z2;          // q^2 of the Z2 system

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1] = (cosz1, phiz1, cosz1f, 
  Double_t       fXL[6];          //!          phiz1f, cosz2f, phiz2f)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(EEZZBases, 1) // Bases for e+e- -> XX process
};

class EEZZSpring;

//_______________________________________________________________________
// =====================
//  class EEZZSpringBuf
// =====================
//-----------------------------------------------------------------------
class EEZZSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  EEZZSpringBuf(const char *name   = "EEZZSpringBuf", 
                const char *title  = "EEZZ Spring test event buffer",
	        EEZZSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEZZSpringBuf() {}

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
  Double_t GetQ2ZZ     ()           const { return fQ2ZZ;      }
  Double_t GetCosZ1    ()           const { return fCosZ1;     }
  Double_t GetPhiZ1    ()           const { return fPhiZ1;     }
  Double_t GetQ2Z1     ()           const { return fQ2Z1;      }
  Double_t GetCosZ1F   ()           const { return fCosZ1F;    }
  Double_t GetPhiZ1F   ()           const { return fPhiZ1F;    }
  Double_t GetQ2Z2     ()           const { return fQ2Z2;      }
  Double_t GetCosZ2F   ()           const { return fCosZ2F;    }
  Double_t GetPhiZ2F   ()           const { return fPhiZ2F;    }

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
  Double_t fQ2ZZ;           // q^2 of the ZZ system
  Double_t fCosZ1;          // cos(theta_Z1) in the ZZ frame
  Double_t fPhiZ1;          // phi_Z1        in the ZZ frame
  Double_t fQ2Z1;           // q^2 of the Z1 system
  Double_t fCosZ1F;         // cos(theta_fub) in the Z1 frame
  Double_t fPhiZ1F;         // phi_fub        in the Z1 frame
  Double_t fQ2Z2;           // q^2 of the Z2 system
  Double_t fCosZ2F;         // cos(theta_fu) in the Z2 frame
  Double_t fPhiZ2F;         // phi_fu        in the Z2 frame

  ClassDef(EEZZSpringBuf, 1)  // EEZZSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class EEZZSpring
// =====================
//-----------------------------------------------------------------------
class EEZZSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   EEZZSpring(const char *name  = "EEZZSpring", 
              const char *title = "EEZZ Spring test",
              EEZZBases  *bases = 0);
   virtual ~EEZZSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(EEZZSpring, 1)  // EEZZSpring class
};
#endif
