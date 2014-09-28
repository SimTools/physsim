#ifndef ENWSPRING_H
#define ENWSPRING_H
//*****************************************************************************
//* =====================
//*  ENWSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> enW generator
//*
//* (Update Record)
//*    2014/09/25  K.Fujii	Original version.
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
//  class ENWBases
// =====================
//-----------------------------------------------------------------------
class ENWBases : public JSFBases {
friend class ENWSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ENWBases(const char *name  = "ENWBases", 
           const char *title = "ENW Bases");
  virtual ~ENWBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Int_t    GetCP       ()           const { return fCP;        }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2W      ()           const { return fQ2W;       }
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
  Complex_t AmpEEtoENW   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &ef,
                          const HELFermion &nu,
                          const HELVector  &wf);

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
  Int_t    fWModesLo;       // W decay mode lo;
  Int_t    fWModesHi;       // W decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

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
  Double_t       fXi;             // xi
  Double_t       fEta1;           // eta_1 for e- -> e-
  Double_t       fEta2;           // eta_2 for e+ -> nb
  Double_t       fPhi1;           // phi_1 for e- -> e-
  Double_t       fPhi2;           // phi_2 for e+ -> nb
  Double_t       fQ2W;            // q^2 of the final state W
  GENDecayMode  *fWModePtr;       // pointer to W decay mode
  GENPDTEntry   *f3Ptr;           // point to up   type W daughter
  GENPDTEntry   *f4Ptr;           // point to down type W daughter
  Double_t       fCosF;           // cos(theta_fu) in the W frame
  Double_t       fPhiF;           // phi_fu        in the W frame
  Int_t          fCP;             // CP flip: (-1,+1) = (nu e+ W-, e- nb W+)
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWMode;          // W decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [e-, nb,  fu, fdb]
  Double_t       fM[4];           // [0,1,2,3] = [me, mn, mfu, mfd]
  Double_t       fSh1;            // sh1
  Double_t       fCh1;            // ch1

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fWDecayMode;     // decay mode selector for W
  Double_t       fXCP;            // CP flip: (nu e+ W-, e- nb W+)
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> nb
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> nb
  Double_t       fXQ2W;           // q^2 of the final state W
  Double_t       fXCosF;          // cos(theta_f) in the W frame
  Double_t       fXPhiF;          // phi_f        in the W frame

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

  ClassDef(ENWBases, 1) // Bases for e+e- -> XX process
};

class ENWSpring;

//_______________________________________________________________________
// =====================
//  class ENWSpringBuf
// =====================
//-----------------------------------------------------------------------
class ENWSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ENWSpringBuf(const char *name   = "ENWSpringBuf", 
               const char *title  = "ENW Spring test event buffer",
	        ENWSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ENWSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Int_t    GetCP       ()           const { return fCP;        }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2W      ()           const { return fQ2W;       }
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
  Int_t    fCP;             // CP flip: (-1,+1) = (nu e+ W-, e- nb W+)
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung
  Double_t fXi;             // xi
  Double_t fEta1;           // eta_1 for e- -> e-
  Double_t fEta2;           // eta_2 for e+ -> e+
  Double_t fPhi1;           // phi_1 for e- -> e-
  Double_t fPhi2;           // phi_2 - phi_1 for e+ -> e+
  Double_t fQ2W;            // q^2 of the final state W
  Double_t fCosF;           // cos(theta_f) in the W frame
  Double_t fPhiF;           // phi_f        in the W frame

  ClassDef(ENWSpringBuf, 1)  // ENWSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ENWSpring
// =====================
//-----------------------------------------------------------------------
class ENWSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ENWSpring(const char *name  = "ENWSpring", 
             const char *title = "ENW Spring test",
              ENWBases  *bases = 0);
   virtual ~ENWSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ENWSpring, 1)  // ENWSpring class
};
#endif
