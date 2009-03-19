#ifndef WHWHSPRING_H
#define WHWHSPRING_H
//*****************************************************************************
//* =====================
//*  WHWHSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- -> XX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/28  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class WHBoson;
class AHBoson;

//_______________________________________________________________________
// =====================
//  class WHWHBases
// =====================
//-----------------------------------------------------------------------
class WHWHBases : public JSFBases {
friend class WHWHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WHWHBases(const char *name  = "WHWHBases", 
            const char *title = "WHWH Bases");
  virtual ~WHWHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetF         ()           const { return fF;          }
  Double_t GetKappaL    ()           const { return fKappaL;     }
  Double_t GetMass      ()           const { return fMass;       }
  Double_t GetMassDM    ()           const { return fMassDM;     }
  Double_t GetMassT     ()           const { return fMassT;      }
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
                          const HELVector  &xm,
                          const HELVector  &xp);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fF;              // SU(5)->SO(5) breaking scale
  Double_t fKappaL;         // heavy lepton mass parameter
  Double_t fMass;           // m_x    : mass  of X
  Double_t fMassDM;         // m_dm   : mass  of DM
  Double_t fMassT;          // m_t    : mass  of t-channel heavy object

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam Spread relative to nominal Ebm
  Double_t fPole;           // electron polarization
  Int_t    fWmModesLo;      // W- decay mode lo;
  Int_t    fWmModesHi;      // W- decay mode hi;
  Int_t    fWpModesLo;      // W+ decay mode lo;
  Int_t    fWpModesHi;      // W+ decay mode hi;

  // ----------------
  //  Particle Data
  // ----------------
  WHBoson      *fXBosonPtr;       //! PD table entry of "X"
  AHBoson      *fDMBosonPtr;      //! PD table entry of "DM"
  GENPDTWBoson *fW1BosonPtr;      //! PD table entry of "W1"
  GENPDTWBoson *fW2BosonPtr;      //! PD table entry of "W2"
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
  ANL4DVector    fP[6];           // [0,1,2,3] = [ah1 , fwm1, fwm2, ah2, fwp1, fwp2]
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

  ClassDef(WHWHBases, 1) // Bases for e+e- -> XX process
};

class WHWHSpring;

//_______________________________________________________________________
// =====================
//  class WHWHSpringBuf
// =====================
//-----------------------------------------------------------------------
class WHWHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  WHWHSpringBuf(const char   *name   = "WHWHSpringBuf", 
                const char   *title  = "WHWH Spring test event buffer",
	          WHWHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~WHWHSpringBuf() {}

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

  ClassDef(WHWHSpringBuf, 1)  // WHWHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class WHWHSpring
// =====================
//-----------------------------------------------------------------------
class WHWHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   WHWHSpring(const char   *name  = "WHWHSpring", 
              const char   *title = "WHWH Spring test",
                WHWHBases  *bases = 0);
   virtual ~WHWHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(WHWHSpring, 1)  // WHWHSpring class
};

//_______________________________________________________________________
// =====================
//  class WHBoson
// =====================
//-----------------------------------------------------------------------

class WHBoson: public GENPDTEntry {
public:
   WHBoson(Double_t fF = 580.);
   virtual ~WHBoson() {}

   Double_t  GetF     () const { return fF;      }
   Double_t  GetSh    () const { return fSh;     }
   Double_t  GetMassDM() const { return fMassDM; }

private:
   void     Initialize();
   Double_t GamToVV(Double_t m1, Double_t m2, Double_t a);

private:
   Double_t fF;      // SU(5)->SO(5) breaking scale
   Double_t fSh;     // sh parameter
   Double_t fMassDM; // daughter DM mass

   ClassDef(WHBoson, 1)  // X boson class
};

//_______________________________________________________________________
// =====================
//  class AHBoson
// =====================
//-----------------------------------------------------------------------

class AHBoson: public GENPDTEntry {
public:
   AHBoson(Double_t m      = 81.9);
  virtual ~AHBoson() {}

   ClassDef(AHBoson, 1)  // X boson class
};
#endif
