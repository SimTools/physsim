#ifndef ZZHSPRING_H
#define ZZHSPRING_H
//*****************************************************************************
//* =====================
//*  ZZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ZZH generator
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
//  class ZZHBases
// =====================
//-----------------------------------------------------------------------
class ZZHBases : public JSFBases {
friend class ZZHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZZHBases(const char *name  = "ZZHBases", 
           const char *title = "ZZH Bases");
  virtual ~ZZHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass       ()           const { return fMass;        }
  Double_t GetEcmInit    ()           const { return fEcmInit;     }

  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2ZZ       ()           const { return fQ2ZZ;        }
  Double_t GetCosThetaZ  ()           const { return fCosThetaZ;   }
  Double_t GetPhiZ       ()           const { return fPhiZ;        }
  Double_t GetQ2Z1       ()           const { return fQ2Z1;        }
  Double_t GetCosThetaZ1F()           const { return fCosThetaZ1F; }
  Double_t GetPhiZ1F     ()           const { return fPhiZ1F;      }
  Double_t GetQ2Z2       ()           const { return fQ2Z2;        }
  Double_t GetCosThetaZ2F()           const { return fCosThetaZ2F; }
  Double_t GetPhiZ2F     ()           const { return fPhiZ2F;      }
  Double_t GetZBoost     ()           const { return fZBoost;      }
  Double_t GetEcmIP      ()           const { return fEcmIP;       }

  void     SetMass       (Double_t m      ) { fMass      = m;      }
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
  Complex_t AmpEEtoZZH   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELVector  &z1,
                          const HELVector  &z2,
                          const HELScalar  &hs);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization
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
  Double_t       fQ2ZZ;           // q^2 of final state Z2Z1
  Double_t       fQ2Z1;           // q^2 of final state Z1
  Double_t       fQ2Z2;           // q^2 of final state Z2
  GENDecayMode  *fZ1ModePtr;      // pointer to Z1 decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st Z1 daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd Z1 daughter
  GENDecayMode  *fZ2ModePtr;      // pointer to Z2 decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z2 daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z2 daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [5];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZ1Mode;         // Z1 decay mode
  Int_t          fZ2Mode;         // Z2 decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[5];           // [0,1,2,3,4] = [h , fz11,fz12, fz21,fz22]
  Double_t       fM[5];           // [0,1,2,3,4] = [mh,   m1,  m2,   m3,  m4]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZ1DecayMode;    // decay mode selector for Z1
  Double_t       fZ2DecayMode;    // decay mode selector for Z2
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2ZZ;          // q^2 of final state ZZ
  Double_t       fCosThetaZ;      // cos(theta_Z1) in ZZ frame
  Double_t       fPhiZ;           // phi_Z1        in ZZ frame
  Double_t       fXQ2Z1;          // q^2 of final state Z1
  Double_t       fCosThetaZ1F;    // cos(theta_f) in Z1  frame
  Double_t       fPhiZ1F;         // phi_f        in Z1  frame
  Double_t       fXQ2Z2;          // q^2 of final state Z2
  Double_t       fCosThetaZ2F;    // cos(theta_f) in Z2  frame
  Double_t       fPhiZ2F;         // phi_f        in Z2  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[8];          //! [0,1,2,3,4,5,6,7] = (cosx, phix,
  Double_t       fXL[8];          //! cosw, phiw, cosz1f, phiz1f, cosz2f, phiz2f)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ZZHBases, 1) // Bases for e+e- -> XX process
};

class ZZHSpring;

//_______________________________________________________________________
// =====================
//  class ZZHSpringBuf
// =====================
//-----------------------------------------------------------------------
class ZZHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZZHSpringBuf(const char *name   = "ZZHSpringBuf", 
               const char *title  = "ZZH Spring test event buffer",
	        ZZHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZZHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosTheta   ()           const { return fCosTheta;    }
  Double_t GetPhi        ()           const { return fPhi;         }
  Double_t GetQ2ZZ       ()           const { return fQ2ZZ;        }
  Double_t GetCosThetaZ  ()           const { return fCosThetaZ;   }
  Double_t GetPhiZ       ()           const { return fPhiZ;        }
  Double_t GetQ2Z1       ()           const { return fQ2Z1;        }
  Double_t GetCosThetaZ1F()           const { return fCosThetaZ1F; }
  Double_t GetPhiZ1F     ()           const { return fPhiZ1F;      }
  Double_t GetQ2Z2       ()           const { return fQ2Z2;        }
  Double_t GetCosThetaZ2F()           const { return fCosThetaZ2F; }
  Double_t GetPhiZ2F     ()           const { return fPhiZ2F;      }
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
  Double_t fQ2ZZ;           // q^2 of ZZ system
  Double_t fCosThetaZ;      // cos(theta_Z1) in H  frame
  Double_t fPhiZ;           // phi_Z1        in H  frame
  Double_t fQ2Z1;           // q^2 of final state Z1
  Double_t fCosThetaZ1F;    // cos(theta_f) in Z1  frame
  Double_t fPhiZ1F;         // phi_f        in Z1  frame
  Double_t fQ2Z2;           // q^2 of final state Z1
  Double_t fCosThetaZ2F;    // cos(theta_f) in Z2  frame
  Double_t fPhiZ2F;         // phi_f        in Z2  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(ZZHSpringBuf, 1)  // ZZHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ZZHSpring
// =====================
//-----------------------------------------------------------------------
class ZZHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ZZHSpring(const char *name  = "ZZHSpring", 
             const char *title = "ZZH Spring test",
               ZZHBases *bases = 0);
   virtual ~ZZHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ZZHSpring, 1)  // ZZHSpring class
};
#endif
