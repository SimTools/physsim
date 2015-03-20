#ifndef WWZSPRING_H
#define WWZSPRING_H
//*****************************************************************************
//* =====================
//*  WWZSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> WWZ generator
//*
//* (Update Record)
//*    2014/09/16  K.Fujii	Original version.
//*    2015/03/20  K.Fujii	Implemented anomalous HVV couplings.
//*
//*****************************************************************************

#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class HBoson;

//_______________________________________________________________________
// =====================
//  class WWZBases
// =====================
//-----------------------------------------------------------------------
class WWZBases : public JSFBases {
friend class WWZSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WWZBases(const char *name  = "WWZBases", 
          const char *title = "WWZ Bases");
  virtual ~WWZBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass       ()           const { return fMass;        }
  Double_t GetLambda     ()           const { return fLambda;      }
  Double_t GetA          ()           const { return fA;           }
  Double_t GetB          ()           const { return fB;           }
  Double_t GetBtilde     ()           const { return fBtilde;      }
  Double_t GetEcmInit    ()           const { return fEcmInit;     }

  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2WW       ()           const { return fQ2WW;        }
  Double_t GetCosThetaWm ()           const { return fCosThetaWm;  }
  Double_t GetPhiWm      ()           const { return fPhiWm;       }
  Double_t GetQ2Wm       ()           const { return fQ2Wm;        }
  Double_t GetCosThetaWmF()           const { return fCosThetaWmF; }
  Double_t GetPhiWmF     ()           const { return fPhiWmF;      }
  Double_t GetQ2Wp       ()           const { return fQ2Wp;        }
  Double_t GetCosThetaWpF()           const { return fCosThetaWpF; }
  Double_t GetPhiWpF     ()           const { return fPhiWpF;      }
  Double_t GetQ2Z        ()           const { return fQ2Z;         }
  Double_t GetCosThetaZF ()           const { return fCosThetaZF;  }
  Double_t GetPhiZF      ()           const { return fPhiZF;       }
  Double_t GetZBoost     ()           const { return fZBoost;      }
  Double_t GetEcmIP      ()           const { return fEcmIP;       }

  void     SetMass       (Double_t m      ) { fMass      = m;      }
  void     SetLambda     (Double_t l      ) { fLambda    = l;      }
  void     SetA          (Double_t a      ) { fA         = a;      }
  void     SetB          (Double_t b      ) { fB         = b;      }
  void     SetBtilde     (Double_t bt     ) { fBtilde    = bt;     }
  void     SetEcmInit    (Double_t ecm    ) { fEcmInit   = ecm;    }
  void     SetISR        (Bool_t b = kTRUE) { fISR       = b;      }
  void     SetBeamStr    (Bool_t b = kTRUE) { fBeamStr   = b;      }
  void     SetBeamWidth  (Double_t w      ) { fBeamWidth = w;      }
  void     SetPole       (Double_t p      ) { fPole      = p;      }

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

  Double_t  DSigmaDX     (GENBranch &cmbranch);
  Double_t  AmpSquared   (GENBranch &cmbranch);
  Complex_t FullAmplitude();
  Complex_t AmpEEtoWWZ   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELVector  &wm,
                          const HELVector  &wp,
			  const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h     : mass  of H
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
  Int_t    fZModesLo;       // Z  decay mode lo;
  Int_t    fZModesHi;       // Z  decay mode hi;
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
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z0"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2WW;           // q^2 of final state W-W+ system
  Double_t       fQ2Wm;           // q^2 of final state W-
  Double_t       fQ2Wp;           // q^2 of final state W+
  Double_t       fQ2Z;            // q^2 of final state Z 
  GENDecayMode  *fWmModePtr;      // pointer to W- decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st W- daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd W- daughter
  GENDecayMode  *fWpModePtr;      // pointer to W+ decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st W+ daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd W+ daughter
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f5Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f6Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWmMode;         // W- decay mode
  Int_t          fWpMode;         // W+ decay mode
  Int_t          fZMode;          // Z  decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [fwm1,fwm2, fwp1,fwp2, fz1, fz2]
  Double_t       fM[6];           // [0,1,2,3,4,5] = [  m1,  m2,   m3,  m4,  m5,  m6]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fWmDecayMode;    // decay mode selector for W-
  Double_t       fWpDecayMode;    // decay mode selector for W+
  Double_t       fZDecayMode;     // decay mode selector for Z 
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2WW;          // q^2 of final state W-W+ system
  Double_t       fCosThetaWm;     // cos(theta_W-) in W-W+  frame
  Double_t       fPhiWm;          // phi_W-        in W-W+  frame
  Double_t       fXQ2Wm;          // q^2 of final state W-
  Double_t       fCosThetaWmF;    // cos(theta_f) in W-  frame
  Double_t       fPhiWmF;         // phi_f        in W-  frame
  Double_t       fXQ2Wp;          // q^2 of final state W+
  Double_t       fCosThetaWpF;    // cos(theta_f) in W+  frame
  Double_t       fPhiWpF;         // phi_f        in W+  frame
  Double_t       fXQ2Z;           // q^2 of final state Z 
  Double_t       fCosThetaZF;     // cos(theta_f) in Z  frame
  Double_t       fPhiZF;          // phi_f        in Z  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[10];         //! [0,1,2,3,4,5,6,7,8,9] = (cosx, phix, cosw, phiw,
  Double_t       fXL[10];         //! coswmf, phiwmf, coswpf, phiwpf, coszf, phizf)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(WWZBases, 1) // Bases for e+e- -> XX process
};

class WWZSpring;

//_______________________________________________________________________
// =====================
//  class WWZSpringBuf
// =====================
//-----------------------------------------------------------------------
class WWZSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WWZSpringBuf(const char *name   = "WWZSpringBuf", 
               const char *title  = "WWZ Spring test event buffer",
	        WWZSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~WWZSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2WW       ()           const { return fQ2WW;        }
  Double_t GetCosThetaWm ()           const { return fCosThetaWm;  }
  Double_t GetPhiWm      ()           const { return fPhiWm;       }
  Double_t GetQ2Wm       ()           const { return fQ2Wm;        }
  Double_t GetCosThetaWmF()           const { return fCosThetaWmF; }
  Double_t GetPhiWmF     ()           const { return fPhiWmF;      }
  Double_t GetQ2Wp       ()           const { return fQ2Wp;        }
  Double_t GetCosThetaWpF()           const { return fCosThetaWpF; }
  Double_t GetPhiWpF     ()           const { return fPhiWpF;      }
  Double_t GetQ2Z        ()           const { return fQ2Z;         }
  Double_t GetCosThetaZF ()           const { return fCosThetaZF;  }
  Double_t GetPhiZF      ()           const { return fPhiZF;       }
  Double_t GetZBoost     ()           const { return fZBoost;      }
  Double_t GetEcmIP      ()           const { return fEcmIP;       }

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
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fQ2WW;           // q^2 of final state W-W+ system
  Double_t fCosThetaWm;     // cos(theta_W-) in W-W+  frame
  Double_t fPhiWm;          // phi_W-        in W-W+  frame
  Double_t fQ2Wm;           // q^2 of final state W-
  Double_t fCosThetaWmF;    // cos(theta_f) in W-  frame
  Double_t fPhiWmF;         // phi_f        in W-  frame
  Double_t fQ2Wp;           // q^2 of final state W-
  Double_t fCosThetaWpF;    // cos(theta_f) in W+  frame
  Double_t fPhiWpF;         // phi_f        in W+  frame
  Double_t fQ2Z;            // q^2 of final state Z
  Double_t fCosThetaZF;     // cos(theta_f) in Z  frame
  Double_t fPhiZF;          // phi_f        in Z  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(WWZSpringBuf, 1)  // WWZSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class WWZSpring
// =====================
//-----------------------------------------------------------------------
class WWZSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   WWZSpring(const char *name  = "WWZSpring", 
             const char *title = "WWZ Spring test",
               WWZBases *bases = 0);
   virtual ~WWZSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(WWZSpring, 1)  // WWZSpring class
};
#endif
