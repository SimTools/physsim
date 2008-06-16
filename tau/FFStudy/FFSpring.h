#ifndef __FFSPRING__
#define __FFSPRING__
//*****************************************************************************
//* =====================
//*  FFSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> FF generator
//*
//* (Update Record)
//*    2008/06/05  K.Fujii	Original version.
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
//  class FFBases
// =====================
//-----------------------------------------------------------------------
class FFBases : public JSFBases {
friend class FFSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  FFBases(const char *name  = "FFBases", 
            const char *title = "FF Bases");
  virtual ~FFBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Int_t    GePID       ()           const { return fPID;       }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetPID      (Int_t    id     ) { fPID     = id;     }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit = ecm;    }
  void     SetISR      (Bool_t b = kTRUE) { fISR     = b;      }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr = b;      }
  void     SetPole     (Double_t p      ) { fPole    = p;      }

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
  Complex_t AmpEEtoFF    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &fm,
                          const HELFermion &fp);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Int_t    fPID;            // PID
  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fPole;           // electron polarization

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
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[2];           // [0,1,2] = [f, fbar]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
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

  ClassDef(FFBases, 1) // Bases for e+e- -> XX process
};

class FFSpring;

//_______________________________________________________________________
// =====================
//  class FFSpringBuf
// =====================
//-----------------------------------------------------------------------
class FFSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  FFSpringBuf(const char *name   = "FFSpringBuf", 
                const char *title  = "FF Spring test event buffer",
	        FFSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~FFSpringBuf() {}

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

  ClassDef(FFSpringBuf, 1)  // FFSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class FFSpring
// =====================
//-----------------------------------------------------------------------
class FFSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   FFSpring(const char *name  = "FFSpring", 
              const char *title = "FF Spring test",
              FFBases  *bases = 0);
   virtual ~FFSpring();

   ClassDef(FFSpring, 1)  // FFSpring class
};
#endif
