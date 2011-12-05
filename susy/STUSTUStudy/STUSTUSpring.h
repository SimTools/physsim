#ifndef STUSTUSPRING_H
#define STUSTUSPRING_H
//*****************************************************************************
//* =====================
//*  STUSTUSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> STUSTU generator
//*
//* (Update Record)
//*    2011/02/13  K.Fujii      Original version.
//*    2011/12/05  T.Suehara    Option to set random seed
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
//  class STUSTUBases
// =====================
//-----------------------------------------------------------------------
class STUSTUBases : public JSFBases {
friend class STUSTUSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  STUSTUBases(const char *name  = "STUSTUBases", 
              const char *title = "STUSTU Bases");
  virtual ~STUSTUBases();

  // ----------------------
  //  Set random seed
  // ----------------------
  void SetSeed(int seed){GetRan1()->SetSeed(seed);}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMassStau ()           const { return fMassStau;  }
  Double_t GetMassDM   ()           const { return fMassDM;    }
  Double_t GetThetaMix ()           const { return fThetaMix;  }
  Double_t GetHelTauM  ()           const { return fHelTauM;   }
  Double_t GetHelTauP  ()           const { return fHelTauP;   }
  Double_t GetCTau     ()           const { return fCTau;      }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMassStau (Double_t m      ) { fMassStau  = m;    }
  void     SetMassDM   (Double_t m      ) { fMassDM    = m;    }
  void     SetThetaMix (Double_t t      ) { fThetaMix  = t;    }
  void     SetHelTauM  (Double_t pol    ) { fHelTauM   = pol;  }
  void     SetHelTauP  (Double_t pol    ) { fHelTauP   = pol;  }
  void     SetCTau     (Double_t ctau   ) { fCTau      = ctau; }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPole     (Double_t p      ) { fPole      = p;    }

  static void EnableStauDecay(Bool_t b = kTRUE) { fgEnableStauDecay = b; }

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
  Complex_t AmpEEtoSTUSTU(const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &stm,
                          const HELScalar  &stp);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMassStau;       // stau mass
  Double_t fMassDM;         // dark matter mass
  Double_t fThetaMix;       // mixing angle: st_1 = stR*cos(theta)+stL*sin(theta)
  Double_t fHelTauM;        // tau- polarization
  Double_t fHelTauP;        // tau- polarization
  Double_t fCTau;           // (c tau) [cm] of stau
  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization

  static Bool_t   fgEnableStauDecay;// (0,1)=(no decay, enable stau -> tau DM)

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTEntry  *fFPtr;            //! PD table entry of "f"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [2];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[2];           // [0,1] = [stau-, stau+]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[2];          //! [0,1] = (cosx, phix)
  Double_t       fXL[2];          //!

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(STUSTUBases, 1) // Bases for e+e- -> XX process
};

class STUSTUSpring;

//_______________________________________________________________________
// =====================
//  class STUSTUSpringBuf
// =====================
//-----------------------------------------------------------------------
class STUSTUSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  STUSTUSpringBuf(const char *name   = "STUSTUSpringBuf", 
                const char *title  = "STUSTU Spring test event bustustuer",
	        STUSTUSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~STUSTUSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
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
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(STUSTUSpringBuf, 1)  // STUSTUSpring event bustustuer
};

//_______________________________________________________________________
// =====================
//  class STUSTUSpring
// =====================
//-----------------------------------------------------------------------
class STUSTUSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   STUSTUSpring(const char *name  = "STUSTUSpring", 
                const char *title = "STUSTU Spring test",
                STUSTUBases  *bases = 0);
   virtual ~STUSTUSpring();

   ClassDef(STUSTUSpring, 1)  // STUSTUSpring class
};
#endif
