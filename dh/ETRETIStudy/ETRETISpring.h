#ifndef ETRETISPRING_H
#define ETRETISPRING_H
//*****************************************************************************
//* =====================
//*  ETRETISpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- -> eta_r eta_i generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/18  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class ETRBoson;
class ETIBoson;

//_______________________________________________________________________
// =====================
//  class ETRETIBases
// =====================
//-----------------------------------------------------------------------
class ETRETIBases : public JSFBases {
friend class ETRETISpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ETRETIBases(const char *name  = "ETRETIBases", 
              const char *title = "ETRETI Bases");
  virtual ~ETRETIBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass      ()           const { return fMass;       }
  Double_t GetMassDM    ()           const { return fMassDM;     }
  Double_t GetEcmInit   ()           const { return fEcmInit;    }

  Double_t GetQ2XD      ()           const { return fQ2XD;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X       ()           const { return fQ2X;        }
  Double_t GetCosThetaZ ()           const { return fCosThetaZ;  }
  Double_t GetPhiZ      ()           const { return fPhiZ;       }
  Double_t GetQ2Z       ()           const { return fQ2Z;        }
  Double_t GetCosThetaF ()           const { return fCosThetaF;  }
  Double_t GetPhiF      ()           const { return fPhiF;       }
  Double_t GetZBoost    ()           const { return fZBoost;     }
  Double_t GetEcmIP     ()           const { return fEcmIP;      }

  void     SetMass     (Double_t m      )  { fMass    = m;       }
  void     SetMassDM   (Double_t m      )  { fMassDM  = m;       }
  void     SetEcmInit  (Double_t ecm    )  { fEcmInit = ecm;     }
  void     SetISR      (Bool_t b = kTRUE)  { fISR     = b;       }
  void     SetBeamStr  (Bool_t b = kTRUE)  { fBeamStr = b;       }
  void     SetPole     (Double_t p      )  { fPole    = p;       }

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
  Complex_t AmpEEtoXD    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &xr,
                          const HELScalar  &xi);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_x    : mass  of X
  Double_t fMassDM;         // m_dm   : mass  of DM

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fPole;           // electron polarization

  // ----------------
  //  Particle Data
  // ----------------
  ETRBoson     *fXBosonPtr;       //! PD table entry of "X"
  ETIBoson     *fDMBosonPtr;      //! PD table entry of "DM"
  GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"
  GENPDTPhoton *fPhotonPtr;       //! PD table entry of "Gamma"
  GENPDTGluon  *fGluonPtr;        //! PD table entry of "Gluon"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2XD;           // q^2 of XD system
  Double_t       fQ2X;            // q^2 of final state X
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd Z daughter
  Double_t       fQ2Z;            // q^2 of final state Z
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2] = [et01, fz1, fz2, eta02]
  Double_t       fM[4];           // [0,1,2] = [mdm , m1  , m2, mdm  ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaZ;      // cos(theta_Z) in X  frame
  Double_t       fPhiZ;           // phi_Z        in X  frame
  Double_t       fXQ2X;           // q^2 of final state X
  Double_t       fCosThetaF;      // cos(theta_f1) in Z  frame
  Double_t       fPhiF;           // phi_f1        in Z  frame
  Double_t       fXQ2Z;           // q^2 of final state Z

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1,2,3,4,5] = (cosx, phix,
  Double_t       fXL[6];          //!      cosz, phiz, cosf1, phif2)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ETRETIBases, 1) // Bases for e+e- -> XD process
};

class ETRETISpring;

//_______________________________________________________________________
// =====================
//  class ETRETISpringBuf
// =====================
//-----------------------------------------------------------------------
class ETRETISpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ETRETISpringBuf(const char   *name   = "ETRETISpringBuf", 
                  const char   *title  = "ETRETI Spring test event buffer",
	          ETRETISpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ETRETISpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2XD      ()           const { return fQ2XD;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X       ()           const { return fQ2X;        }
  Double_t GetCosThetaZ ()           const { return fCosThetaZ;  }
  Double_t GetPhiZ      ()           const { return fPhiZ;       }
  Double_t GetQ2Z       ()           const { return fQ2Z;        }
  Double_t GetCosThetaF ()           const { return fCosThetaF;  }
  Double_t GetPhiF      ()           const { return fPhiF;       }
  Double_t GetZBoost    ()           const { return fZBoost;     }
  Double_t GetEcmIP     ()           const { return fEcmIP;      }

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
  Double_t fQ2XD;            // q^2 of XD system
  Double_t fCosTheta;        // cos(theta_x1) in cm  frame
  Double_t fPhi;             // phi_x1        in cm  frame
  Double_t fQ2X;             // q^2 of final state X
  Double_t fCosThetaZ;       // cos(theta_Z) in X  frame
  Double_t fPhiZ;            // phi_Z        in X  frame
  Double_t fQ2Z;             // q^2 of final state Z
  Double_t fCosThetaF;       // cos(theta_f1) in Z  frame
  Double_t fPhiF;            // phi_f1        in Z  frame
  Double_t fZBoost;          // p_z(cm)      in lab frame
  Double_t fEcmIP;           // Ecm after B-strahlung

  ClassDef(ETRETISpringBuf, 1)  // ETRETISpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ETRETISpring
// =====================
//-----------------------------------------------------------------------
class ETRETISpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ETRETISpring(const char   *name  = "ETRETISpring", 
                const char   *title = "ETRETI Spring test",
                ETRETIBases  *bases = 0);
   virtual ~ETRETISpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ETRETISpring, 1)  // ETRETISpring class
};

//_______________________________________________________________________
// =====================
//  class ETRBoson
// =====================
//-----------------------------------------------------------------------

class ETRBoson: public GENPDTEntry {
public:
   ETRBoson(Double_t m    = 180.,
            Double_t mdm  =  60.);
  virtual ~ETRBoson() {}

private:
   void     Initialize();
   Double_t GamToSV(Double_t m1, Double_t m2, Double_t a, Double_t color);

private:
   Double_t fMassDM; // daughter eta_I mass

   ClassDef(ETRBoson, 1)  // eta_R boson class
};

//_______________________________________________________________________
// =====================
//  class ETIBoson
// =====================
//-----------------------------------------------------------------------

class ETIBoson: public GENPDTEntry {
public:
   ETIBoson(Double_t m      = 60.);
  virtual ~ETIBoson() {}

   ClassDef(ETIBoson, 1)  // eta_I boson class
};
#endif
