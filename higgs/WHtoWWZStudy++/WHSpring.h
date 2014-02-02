#ifndef WHSPRING_H
#define WHSPRING_H
//*****************************************************************************
//* =====================
//*  WHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> WH generator
//*
//* (Update Record)
//*    2014/01/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class XBoson;

//_______________________________________________________________________
// =====================
//  class WHBases
// =====================
//-----------------------------------------------------------------------
class WHBases : public JSFBases {
friend class WHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WHBases(const char *name  = "WHBases", 
          const char *title = "WH Bases");
  virtual ~WHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass       ()         const { return fMass;        }
  Double_t GetFhwz       ()         const { return fFhwz;        }
  Double_t GetFhwa       ()         const { return fFhwa;        }
  Double_t GetEcmInit    ()         const { return fEcmInit;     }
  Double_t GetPole       ()         const { return fPole;        }
  Double_t GetPolp       ()         const { return fPolp;        }

  Double_t GetQ2WH       ()         const { return fQ2WH;        }
  Double_t GetCosTheta   ()         const { return fCosTheta;    }
  Double_t GetPhi        ()         const { return fPhi;         }
  Double_t GetQ2W1       ()         const { return fQ2W1;        }
  Double_t GetCosThetaW1F()         const { return fCosThetaW1F; }
  Double_t GetPhiW1F     ()         const { return fPhiW1F;      }
  Double_t GetQ2H        ()         const { return fQ2H;         }
  Double_t GetCosThetaW2 ()         const { return fCosThetaW2;  }
  Double_t GetPhiW2      ()         const { return fPhiW2;       }
  Double_t GetQ2W2       ()         const { return fQ2W2;        }
  Double_t GetCosThetaW2F()         const { return fCosThetaW2F; }
  Double_t GetPhiW2F     ()         const { return fPhiW2F;      }
  Double_t GetQ2Z        ()         const { return fQ2Z;         }
  Double_t GetCosThetaZF ()         const { return fCosThetaZF;  }
  Double_t GetPhiZF      ()         const { return fPhiZF;       }

  Double_t GetZBoost     ()         const { return fZBoost;      }
  Double_t GetEcmIP      ()         const { return fEcmIP;       }
  Int_t    GetCP         ()         const { return fCP;          }

  void     SetMass     (Double_t m      ) { fMass      = m;    }
  void     SetFhwz     (Double_t f      ) { fFhwz      = f;    }
  void     SetFhwa     (Double_t f      ) { fFhwa      = f;    }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPole     (Double_t p      ) { fPole      = p;    }
  void     SetPolp     (Double_t p      ) { fPolp      = p;    }

  // ----------------------
  //   Base class methods
  // ----------------------
  virtual void     Userin ();  // Bases user initialization
  virtual void     Userout();  // Bases user output 

  Double_t Func();     // Bases integration function.

  // ----------------------
  //   Utility methods
  // ----------------------
private:
  void     SelectHelicities(Double_t &weight);

  Double_t  DSigmaDX     (GENBranch &cmbranch);
  Double_t  AmpSquared   (GENBranch &cmbranch);
  Complex_t FullAmplitude();
  Complex_t AmpEEtoWH    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &hs,
                          const HELVector  &wf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fFhwz;           // g_w M_w fFhwz : HWZ coupling
  Double_t fFhwa;           // g_w M_w fFhwa : HWA coupling

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization
  Double_t fPolp;           // positron polarization
  Double_t fFixCP;          // CP combination (W+H-,W-H+)=(1,-1)
  Int_t    fW1ModesLo;      // W1 decay mode lo;
  Int_t    fW1ModesHi;      // W1 decay mode hi;
  Int_t    fW2ModesLo;      // W2 decay mode lo;
  Int_t    fW2ModesHi;      // W2 decay mode hi;
  Int_t    fZModesLo;       // Z  decay mode lo;
  Int_t    fZModesHi;       // Z  decay mode hi;

  // ----------------
  //  Particle Data
  // ----------------
  XBoson       *fXBosonPtr;       //! PD table entry of "X"
  GENPDTWBoson *fW1BosonPtr;      //! PD table entry of "W"
  GENPDTWBoson *fW2BosonPtr;      //! PD table entry of "W"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2WH;           // q^2 of WH system
  Double_t       fQ2H;            // q^2 of final state H
  Double_t       fQ2W1;           // q^2 of final state W
  Double_t       fQ2W2;           // q^2 of final state W
  Double_t       fQ2Z;            // q^2 of final state W
  Int_t          fCP;             // CP combination of final state (W+H-,W-H+)=(+1,-1)
  GENDecayMode  *fW1ModePtr;      // pointer to W1 decay mode
  GENDecayMode  *fW2ModePtr;      // pointer to W2 decay mode
  GENDecayMode  *fZModePtr;       // pointer to Z  decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st W1 daughter (up type fermion)
  GENPDTEntry   *f2Ptr;           // point to 2nd W1 daughter (down type anti-fermion)
  GENPDTEntry   *f3Ptr;           // point to 1st W2 daughter (up type anti-fermion)
  GENPDTEntry   *f4Ptr;           // point to 2nd W2 daughter (down type fermion)
  GENPDTEntry   *f5Ptr;           // point to 1st Z  daughter (fermion)
  GENPDTEntry   *f6Ptr;           // point to 2nd Z  daughter (anti-fermion)
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fW1Mode;         // W1 decay mode
  Int_t          fW2Mode;         // W2 decay mode
  Int_t          fZMode;          // Z  decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [fw1, fbw1, fbw2, fw2, fz, fbz]
  Double_t       fM[6];           // [0,1,2,3,4,5] = [ m1,   m2,   m3,  m4, m5,  m6]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fCPFinal;        // final   state CP combination
  Double_t       fW1DecayMode;    // decay mode selector for W1
  Double_t       fW2DecayMode;    // decay mode selector for W2
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2W1;          // q^2 of final state W1
  Double_t       fCosThetaW1F;    // cos(theta_f) in W1  frame
  Double_t       fPhiW1F;         // phi_f        in W1  frame
  Double_t       fXQ2H;           // q^2 of final state H
  Double_t       fCosThetaW2;     // cos(theta_W) in H   frame
  Double_t       fPhiW2;          // phi_W        in H   frame
  Double_t       fXQ2W2;          // q^2 of final state W2
  Double_t       fCosThetaW2F;    // cos(theta_f) in W2  frame
  Double_t       fPhiW2F;         // phi_f        in W2  frame
  Double_t       fXQ2Z;           // q^2 of final state Z
  Double_t       fCosThetaZF;     // cos(theta_f) in Z   frame
  Double_t       fPhiZF;          // phi_f        in Z   frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[10];         //! (cosx, phix, cosf1, phif1, cosw, phiw, cosf3, phif3, cosf5, phif5)
  Double_t       fXL[10];         //!

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(WHBases, 1) // Bases for e+e- -> XX process
};

class WHSpring;

//_______________________________________________________________________
// =====================
//  class WHSpringBuf
// =====================
//-----------------------------------------------------------------------
class WHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WHSpringBuf(const char *name   = "WHSpringBuf", 
              const char *title  = "WH Spring test event buffer",
	        WHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~WHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2WH       ()         const { return fQ2WH;        }
  Double_t GetCosTheta   ()         const { return fCosTheta;    }
  Double_t GetPhi        ()         const { return fPhi;         }
  Double_t GetQ2W1       ()         const { return fQ2W1;        }
  Double_t GetCosThetaW1F()         const { return fCosThetaW1F; }
  Double_t GetPhiW1F     ()         const { return fPhiW1F;      }
  Double_t GetQ2H        ()         const { return fQ2H;         }
  Double_t GetCosThetaW2 ()         const { return fCosThetaW2;  }
  Double_t GetPhiW2      ()         const { return fPhiW2;       }
  Double_t GetQ2W2       ()         const { return fQ2W2;        }
  Double_t GetCosThetaW2F()         const { return fCosThetaW2F; }
  Double_t GetPhiW2F     ()         const { return fPhiW2F;      }
  Double_t GetQ2Z        ()         const { return fQ2Z;         }
  Double_t GetCosThetaZF ()         const { return fCosThetaZF;  }
  Double_t GetPhiZF      ()         const { return fPhiZF;       }
  Double_t GetZBoost     ()         const { return fZBoost;      }
  Double_t GetEcmIP      ()         const { return fEcmIP;       }
  Int_t    GetCP         ()         const { return fCP;          }

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
  Double_t fQ2WH;           // q^2 of WH system
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame

  Double_t fQ2W1;           // q^2 of final state W
  Double_t fCosThetaW1F;    // cos(theta_f) in W   frame
  Double_t fPhiW1F;         // phi_f        in W   frame

  Double_t fQ2H;            // q^2 of final state H
  Double_t fCosThetaW2;     // cos(theta_W) in H   frame
  Double_t fPhiW2;          // phi_W        in H   frame

  Double_t fQ2W2;           // q^2 of final state W
  Double_t fCosThetaW2F;    // cos(theta_f) in W2  frame
  Double_t fPhiW2F;         // phi_f        in W2  frame

  Double_t fQ2Z;            // q^2 of final state W
  Double_t fCosThetaZF;     // cos(theta_f) in Z   frame
  Double_t fPhiZF;          // phi_f        in Z   frame

  Int_t    fCP;             // CP combination of final state (W+H-,W-H+)=(+1,-1)
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(WHSpringBuf, 1)  // WHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class WHSpring
// =====================
//-----------------------------------------------------------------------
class WHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   WHSpring(const char *name  = "WHSpring", 
            const char *title = "WH Spring test",
              WHBases  *bases = 0);
   virtual ~WHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(WHSpring, 1)  // WHSpring class
};
#endif
