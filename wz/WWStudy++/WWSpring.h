#ifndef WWSPRING_H
#define WWSPRING_H
//*****************************************************************************
//* =====================
//*  WWSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> WW generator
//*
//* (Update Record)
//*    2010/04/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

//_______________________________________________________________________
// =====================
//  class WWBases
// =====================
//-----------------------------------------------------------------------
class WWBases : public JSFBases {
friend class WWSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WWBases(const char *name  = "WWBases", 
          const char *title = "WW Bases");
  virtual ~WWBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetEcmInit    ()           const { return fEcmInit;     }

  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2Wm       ()           const { return fQ2Wm;        }
  Double_t GetCosThetaWmF()           const { return fCosThetaWmF; }
  Double_t GetPhiWmF     ()           const { return fPhiWmF;      }
  Double_t GetQ2Wp       ()           const { return fQ2Wp;        }
  Double_t GetCosThetaWpF()           const { return fCosThetaWpF; }
  Double_t GetPhiWpF     ()           const { return fPhiWpF;      }
  Double_t GetZBoost     ()           const { return fZBoost;      }
  Double_t GetEcmIP      ()           const { return fEcmIP;       }

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
  Complex_t AmpEEtoWW    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELVector  &wm,
                          const HELVector  &wp);

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
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z0"
  GENPDTWBoson *fWmBosonPtr;      //! PD table entry of "W-"
  GENPDTWBoson *fWpBosonPtr;      //! PD table entry of "W+"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2Wm;           // q^2 of final state W-
  Double_t       fQ2Wp;           // q^2 of final state W+
  GENDecayMode  *fWmModePtr;      // pointer to W- decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st W- daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd W- daughter
  GENDecayMode  *fWpModePtr;      // pointer to W+ decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st W+ daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd W+ daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWmMode;         // W- decay mode
  Int_t          fWpMode;         // W+ decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [fwm1,fwm2, fwp1,fwp2]
  Double_t       fM[4];           // [0,1,2,3] = [  m1,  m2,   m3,  m4]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fWmDecayMode;    // decay mode selector for W-
  Double_t       fWpDecayMode;    // decay mode selector for W+
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2Wm;          // q^2 of final state W-
  Double_t       fCosThetaWmF;    // cos(theta_f) in W-  frame
  Double_t       fPhiWmF;         // phi_f        in W-  frame
  Double_t       fXQ2Wp;          // q^2 of final state W+
  Double_t       fCosThetaWpF;    // cos(theta_f) in W+  frame
  Double_t       fPhiWpF;         // phi_f        in W+  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1,2,3,4,5] = (cosw, phiw,
  Double_t       fXL[6];          //! coswmf, phiwmf, coswpf, phiwpf)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(WWBases, 1) // Bases for e+e- -> XX process
};

class WWSpring;

//_______________________________________________________________________
// =====================
//  class WWSpringBuf
// =====================
//-----------------------------------------------------------------------
class WWSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WWSpringBuf(const char *name   = "WWSpringBuf", 
              const char *title  = "WW Spring test event buffer",
	        WWSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~WWSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2Wm       ()           const { return fQ2Wm;        }
  Double_t GetCosThetaWmF()           const { return fCosThetaWmF; }
  Double_t GetPhiWmF     ()           const { return fPhiWmF;      }
  Double_t GetQ2Wp       ()           const { return fQ2Wp;        }
  Double_t GetCosThetaWpF()           const { return fCosThetaWpF; }
  Double_t GetPhiWpF     ()           const { return fPhiWpF;      }
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
  Double_t fQ2Wm;           // q^2 of final state W-
  Double_t fCosThetaWmF;    // cos(theta_f) in W-  frame
  Double_t fPhiWmF;         // phi_f        in W-  frame
  Double_t fQ2Wp;           // q^2 of final state W-
  Double_t fCosThetaWpF;    // cos(theta_f) in W+  frame
  Double_t fPhiWpF;         // phi_f        in W+  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(WWSpringBuf, 1)  // WWSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class WWSpring
// =====================
//-----------------------------------------------------------------------
class WWSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   WWSpring(const char *name  = "WWSpring", 
            const char *title = "WW Spring test",
               WWBases *bases = 0);
   virtual ~WWSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(WWSpring, 1)  // WWSpring class
};
#endif
