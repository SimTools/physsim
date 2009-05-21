#ifndef ZDDSPRING_H
#define ZDDSPRING_H
//*****************************************************************************
//* =====================
//*  ZDDSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ZDD generator
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

class DMBoson;

//_______________________________________________________________________
// =====================
//  class ZDDBases
// =====================
//-----------------------------------------------------------------------
class ZDDBases : public JSFBases {
friend class ZDDSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZDDBases(const char *name  = "ZDDBases", 
           const char *title = "ZDD Bases");
  virtual ~ZDDBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetMassH    ()           const { return fMassH;     }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetQ2ZDD    ()           const { return fQ2ZDD;     }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetQ2HH     ()           const { return fQ2HH;      }
  Double_t GetCosThetaH()           const { return fCosThetaH; }
  Double_t GetPhiH     ()           const { return fPhiH;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass      = m;    }
  void     SetMassH    (Double_t m      ) { fMassH     = m;    }
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

  Double_t  DSigmaDX     (GENBranch &cmbranch);
  Double_t  AmpSquared   (GENBranch &cmbranch);
  Complex_t FullAmplitude();
  Complex_t AmpEEtoZDD   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELVector  &h1,
                          const HELVector  &h2,
                          const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_DM   : mass  of DM
  Double_t fMassH;          // m_h    : mass  of H

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization
  Int_t    fZModesLo;       // Z decay mode lo;
  Int_t    fZModesHi;       // Z decay mode hi;
  Double_t fCv;             // Cv

  // ----------------
  //  Particle Data
  // ----------------
  DMBoson      *fDMBosonPtr;      //! PD table entry of "DM"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2ZDD;          // q^2 of ZDD system
  Double_t       fQ2Z;            // q^2 of final state Z
  Double_t       fQ2HH;           // q^2 of final state HH
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [18]; // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [h1, h2, fz1, fz2]
  Double_t       fM[4];           // [0,1,2,4] = [mh, mh, m3 , m4 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2Z;           // q^2 of final state Z
  Double_t       fCosThetaF;      // cos(theta_H) in HH  frame
  Double_t       fPhiF;           // phi_H        in HH  frame
  Double_t       fXQ2HH;          // q^2 of final state Z
  Double_t       fCosThetaH;      // cos(theta_H) in HH  frame
  Double_t       fPhiH;           // phi_H        in HH  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1,2,3,4,5] = (cosx, phix, cosf, phif, cosh, phih)
  Double_t       fXL[6];          //!

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ZDDBases, 1) // Bases for e+e- -> XX process
};

class ZDDSpring;

//_______________________________________________________________________
// =====================
//  class ZDDSpringBuf
// =====================
//-----------------------------------------------------------------------
class ZDDSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZDDSpringBuf(const char *name   = "ZDDSpringBuf", 
               const char *title  = "ZDD Spring test event buffer",
                ZDDSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZDDSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2ZDD    ()           const { return fQ2ZDD;     }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetQ2HH     ()           const { return fQ2HH;      }
  Double_t GetCosThetaH()           const { return fCosThetaH; }
  Double_t GetPhiH     ()           const { return fPhiH;      }
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
  Double_t fQ2ZDD;          // q^2 of ZDD system
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fQ2Z;            // q^2 of final state Z
  Double_t fCosThetaF;      // cos(theta_f) in Z   frame
  Double_t fPhiF;           // phi_f        in Z   frame
  Double_t fQ2HH;           // q^2 of final state HH
  Double_t fCosThetaH;      // cos(theta_f) in HH  frame
  Double_t fPhiH;           // phi_f        in HH  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(ZDDSpringBuf, 1)  // ZDDSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ZDDSpring
// =====================
//-----------------------------------------------------------------------
class ZDDSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ZDDSpring(const char *name  = "ZDDSpring", 
             const char *title = "ZDD Spring test",
              ZDDBases  *bases = 0);
   virtual ~ZDDSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ZDDSpring, 1)  // ZDDSpring class
};

//_______________________________________________________________________
// =====================
//  class DMBoson
// =====================
//-----------------------------------------------------------------------

class DMBoson: public GENPDTEntry {
public:
   DMBoson(Double_t m      = 70.);
   virtual ~DMBoson() {}

   ClassDef(DMBoson, 1)  // eta_I boson class
};
#endif
