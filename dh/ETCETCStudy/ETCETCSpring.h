#ifndef ETCETCSPRING_H
#define ETCETCSPRING_H
//*****************************************************************************
//* =====================
//*  ETCETCSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- -> XX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/17  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class ETCBoson;
class ETIBoson;

//_______________________________________________________________________
// =====================
//  class ETCETCBases
// =====================
//-----------------------------------------------------------------------
class ETCETCBases : public JSFBases {
friend class ETCETCSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ETCETCBases(const char *name  = "ETCETCBases", 
              const char *title = "ETCETC Bases");
  virtual ~ETCETCBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass      ()           const { return fMass;       }
  Double_t GetMassDM    ()           const { return fMassDM;     }
  Double_t GetEcmInit   ()           const { return fEcmInit;    }

  Double_t GetQ2XX      ()           const { return fQ2XX;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X1      ()           const { return fQ2X1;       }
  Double_t GetCosThetaW1()           const { return fCosThetaW1; }
  Double_t GetPhiW1     ()           const { return fPhiW1;      }
  Double_t GetQ2W1      ()           const { return fQ2W1;       }
  Double_t GetCosThetaF1()           const { return fCosThetaF1; }
  Double_t GetPhiF1     ()           const { return fPhiF1;      }
  Double_t GetQ2X2      ()           const { return fQ2X2;       }
  Double_t GetCosThetaW2()           const { return fCosThetaW2; }
  Double_t GetPhiW2     ()           const { return fPhiW2;      }
  Double_t GetQ2W2      ()           const { return fQ2W2;       }
  Double_t GetCosThetaF2()           const { return fCosThetaF2; }
  Double_t GetPhiF2     ()           const { return fPhiF2;      }
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
  Complex_t AmpEEtoXX    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &xm,
                          const HELScalar  &xp);

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
  ETCBoson     *fXBosonPtr;       //! PD table entry of "X"
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
  Double_t       fQ2XX;           // q^2 of XX system
  Double_t       fQ2X1;           // q^2 of final state X1
  Double_t       fQ2X2;           // q^2 of final state X2
  GENDecayMode  *fW1ModePtr;      // pointer to W- decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st W- daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd W- daughter
  Double_t       fQ2W1;           // q^2 of final state W-
  GENDecayMode  *fW2ModePtr;      // pointer to W+ decay mode
  GENPDTEntry   *f4Ptr;           // point to 1st W+ daughter
  GENPDTEntry   *f5Ptr;           // point to 2nd W+ daughter
  Double_t       fQ2W2;           // q^2 of final state W+
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fW1Mode;         // W- decay mode
  Int_t          fW2Mode;         // W+ decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3] = [et01, fwm1, fwm2, et02, fwp1, fwp2]
  Double_t       fM[6];           // [0,1,2,3] = [mdm , m1  , m2  , mdm , m3  , m4  ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fW1DecayMode;    // decay mode selector for W-
  Double_t       fW2DecayMode;    // decay mode selector for W+
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaW1;     // cos(theta_W-) in X1  frame
  Double_t       fPhiW1;          // phi_W-        in X1  frame
  Double_t       fXQ2X1;          // q^2 of final state X1
  Double_t       fCosThetaF1;     // cos(theta_f1) in W-  frame
  Double_t       fPhiF1;          // phi_f1        in W-  frame
  Double_t       fXQ2W1;          // q^2 of final state W-
  Double_t       fCosThetaW2;     // cos(theta_W+) in X2  frame
  Double_t       fPhiW2;          // phi_W+        in X2  frame
  Double_t       fXQ2X2;          // q^2 of final state X2
  Double_t       fCosThetaF2;     // cos(theta_f1) in W+  frame
  Double_t       fPhiF2;          // phi_f1        in W+  frame
  Double_t       fXQ2W2;          // q^2 of final state W+

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[10];         //! [0,1,2,3,4,5] = (cosx, phix,
  Double_t       fXL[10];         //!      cosf1, phif1, cosf2, phif2)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ETCETCBases, 1) // Bases for e+e- -> XX process
};

class ETCETCSpring;

//_______________________________________________________________________
// =====================
//  class ETCETCSpringBuf
// =====================
//-----------------------------------------------------------------------
class ETCETCSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ETCETCSpringBuf(const char   *name   = "ETCETCSpringBuf", 
                  const char   *title  = "ETCETC Spring test event buffer",
	          ETCETCSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ETCETCSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2XX      ()           const { return fQ2XX;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X1      ()           const { return fQ2X1;       }
  Double_t GetCosThetaW1()           const { return fCosThetaW1; }
  Double_t GetPhiW1     ()           const { return fPhiW1;      }
  Double_t GetQ2W1      ()           const { return fQ2W1;       }
  Double_t GetCosThetaF1()           const { return fCosThetaF1; }
  Double_t GetPhiF1     ()           const { return fPhiF1;      }
  Double_t GetQ2X2      ()           const { return fQ2X2;       }
  Double_t GetCosThetaW2()           const { return fCosThetaW2; }
  Double_t GetPhiW2     ()           const { return fPhiW2;      }
  Double_t GetQ2W2      ()           const { return fQ2W2;       }
  Double_t GetCosThetaF2()           const { return fCosThetaF2; }
  Double_t GetPhiF2     ()           const { return fPhiF2;      }
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
  Double_t fQ2XX;            // q^2 of XX system
  Double_t fCosTheta;        // cos(theta_x1) in cm  frame
  Double_t fPhi;             // phi_x1        in cm  frame
  Double_t fQ2X1;            // q^2 of final state X1
  Double_t fCosThetaW1;      // cos(theta_W-) in X1  frame
  Double_t fPhiW1;           // phi_W-        in X1  frame
  Double_t fQ2W1;            // q^2 of final state W-
  Double_t fCosThetaF1;      // cos(theta_f1) in W-  frame
  Double_t fPhiF1;           // phi_f1        in W-  frame
  Double_t fQ2X2;            // q^2 of final state X2
  Double_t fCosThetaW2;      // cos(theta_W+) in X2  frame
  Double_t fPhiW2;           // phi_W+        in X2  frame
  Double_t fQ2W2;            // q^2 of final state W+
  Double_t fCosThetaF2;      // cos(theta_f2) in W+  frame
  Double_t fPhiF2;           // phi_f2        in W+  frame
  Double_t fZBoost;          // p_z(cm)      in lab frame
  Double_t fEcmIP;           // Ecm after B-strahlung

  ClassDef(ETCETCSpringBuf, 1)  // ETCETCSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ETCETCSpring
// =====================
//-----------------------------------------------------------------------
class ETCETCSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ETCETCSpring(const char   *name  = "ETCETCSpring", 
                const char   *title = "ETCETC Spring test",
                ETCETCBases  *bases = 0);
   virtual ~ETCETCSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ETCETCSpring, 1)  // ETCETCSpring class
};

//_______________________________________________________________________
// =====================
//  class ETCBoson
// =====================
//-----------------------------------------------------------------------

class ETCBoson: public GENPDTEntry {
public:
   ETCBoson(Double_t m    = 180.,
            Double_t mdm  =  60.);
  virtual ~ETCBoson() {}

private:
   void     Initialize();
   Double_t GamToSV(Double_t m1, Double_t m2, Double_t a, Double_t color);

private:
   Double_t fMassDM; // daughter eta_I mass

   ClassDef(ETCBoson, 1)  // X boson class
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

   ClassDef(ETIBoson, 1)  // X boson class
};
#endif
