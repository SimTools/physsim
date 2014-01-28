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
//*    2014/01/28  K.Fujii	Original version.
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
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetFhwz     ()           const { return fFhwz;      }
  Double_t GetFhwa     ()           const { return fFhwa;      }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetPole     ()           const { return fPole;      }
  Double_t GetPolp     ()           const { return fPolp;      }

  Double_t GetQ2WH     ()           const { return fQ2WH;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2W      ()           const { return fQ2W;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }
  Int_t    GetCP       ()           const { return fCP;        }

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
  Int_t    fWModesLo;       // W decay mode lo;
  Int_t    fWModesHi;       // W decay mode hi;

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2WH;           // q^2 of WH system
  Double_t       fQ2W;            // q^2 of final state W
  Int_t          fCP;             // CP combination of final state (W+H-,W-H+)=(+1,-1)
  GENDecayMode  *fWModePtr;       // pointer to W decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st W daughter (up type fermion)
  GENPDTEntry   *f4Ptr;           // point to 2nd W daughter (down type anti-fermion)
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [2];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWMode;          // W decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [mh, fw1, fw2]
  Double_t       fM[3];           // [0,1,2] = [mh, m3 , m4 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fCPFinal;        // final   state CP combination
  Double_t       fWDecayMode;     // decay mode selector for W
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2W;           // q^2 of final state W
  Double_t       fCosThetaF;      // cos(theta_f) in W   frame
  Double_t       fPhiF;           // phi_f        in W   frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[4];          //! [0,1,2,3] = (cosx, phix, cosf, phif)
  Double_t       fXL[4];          //!

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
  Double_t GetQ2WH     ()           const { return fQ2WH;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2W      ()           const { return fQ2W;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

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
  Double_t fQ2W;            // q^2 of final state W
  Double_t fCosThetaF;      // cos(theta_f) in W   frame
  Double_t fPhiF;           // phi_f        in W   frame
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
