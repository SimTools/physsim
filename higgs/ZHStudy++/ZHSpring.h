#ifndef __ZHSPRING__
#define __ZHSPRING__
//*****************************************************************************
//* =====================
//*  ZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ZH generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
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
//  class ZHBases
// =====================
//-----------------------------------------------------------------------
class ZHBases : public JSFBases {
friend class ZHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHBases(const char *name  = "ZHBases", 
            const char *title = "ZH Bases");
  virtual ~ZHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetMassA    ()           const { return fMassA;     }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetQ2ZH     ()           const { return fQ2ZH;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass    = m;      }
  void     SetMassA    (Double_t m      ) { fMassA   = m;      }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit = ecm;    }
  void     SetISR      (Bool_t b = kTRUE) { fISR     = b;      }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr = b;      }
  void     SetPole     (Double_t p      ) { fPole    = p;      }

  static void EnableHtoAA(Bool_t b = kTRUE) { fgEnableHtoAA = b; }

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
  Complex_t AmpEEtoZH    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &hs,
                          const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fMassA;          // m_A    : mass  of A

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fPole;           // electron polarization

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2ZH;           // q^2 of ZH system
  Double_t       fQ2Z;            // q^2 of final state Z
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [2];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [mh, fz1, fz2]
  Double_t       fM[3];           // [0,1,2] = [mh, m3 , m4 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2Z;           // q^2 of final state Z
  Double_t       fCosThetaF;      // cos(theta_f) in Z   frame
  Double_t       fPhiF;           // phi_f        in Z   frame

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

  static Bool_t  fgEnableHtoAA;    //! H --> AA decay

  ClassDef(ZHBases, 1) // Bases for e+e- -> XX process
};

class ZHSpring;

//_______________________________________________________________________
// =====================
//  class ZHSpringBuf
// =====================
//-----------------------------------------------------------------------
class ZHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHSpringBuf(const char *name   = "ZHSpringBuf", 
                const char *title  = "ZH Spring test event buffer",
	        ZHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2ZH     ()           const { return fQ2ZH;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
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
  Double_t fQ2ZH;           // q^2 of ZH system
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fQ2Z;            // q^2 of final state Z
  Double_t fCosThetaF;      // cos(theta_f) in Z   frame
  Double_t fPhiF;           // phi_f        in Z   frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(ZHSpringBuf, 1)  // ZHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ZHSpring
// =====================
//-----------------------------------------------------------------------
class ZHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ZHSpring(const char *name  = "ZHSpring", 
              const char *title = "ZH Spring test",
              ZHBases  *bases = 0);
   virtual ~ZHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ZHSpring, 1)  // ZHSpring class
};
#endif
