#ifndef __RSAXSPRING__
#define __RSAXSPRING__
//*****************************************************************************
//* =====================
//*  RSAXSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- -> AX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class RSXBoson;

//_______________________________________________________________________
// =====================
//  class RSAXBases
// =====================
//-----------------------------------------------------------------------
class RSAXBases : public JSFBases {
friend class RSAXSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  RSAXBases(const char *name  = "RSAXBases", 
            const char *title = "RSAX Bases");
  virtual ~RSAXBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetQ2AX     ()           const { return fQ2AX;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2X      ()           const { return fQ2X;       }
  Double_t GetCosThetaA()           const { return fCosThetaA; }
  Double_t GetPhiA     ()           const { return fPhiA;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass    = m;      }
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

  Double_t  DSigmaDX     (GENBranch &cmbranch);
  Double_t  AmpSquared   (GENBranch &cmbranch);
  Complex_t FullAmplitude();
  Complex_t AmpEEtoAX    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &xf,
                          const HELVector  &af);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fLambda;         // lambda
  Double_t fC0;             // C_0
  Double_t fC1;             // C_1
  Double_t fC2;             // C_2
  Double_t fC3;             // C_3

  Double_t fMass;           // m_x    : mass  of X

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fPole;           // electron polarization

  // ----------------
  //  Particle Data
  // ----------------
  RSXBoson     *fXBosonPtr;       //! PD table entry of "X"
  GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"
  GENPDTPhoton *fPhotonPtr;       //! PD table entry of "Gamma"
  GENPDTGluon  *fGluonPtr;        //! PD table entry of "Gluon"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2AX;           // q^2 of AX system
  Double_t       fQ2X;            // q^2 of final state X
  GENDecayMode  *fXModePtr;       // pointer to X decay mode
  GENPDTEntry   *f1Ptr;           // pointer to 1st X daughter
  GENPDTEntry   *f2Ptr;           // pointer to 2nd X daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [3];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fXMode;          // X decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [ax1, ax2,  a]
  Double_t       fM[3];           // [0,1,2] = [m1 , m2 , m3]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fXDecayMode;     // decay mode selector for X
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaA;      // cos(theta_a) in X   frame
  Double_t       fPhiA;           // phi_a        in X   frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[4];          //! [0,1,2,3] = (cosx, phix,
  Double_t       fXL[4];          //!              cosa, phia)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(RSAXBases, 1) // Bases for e+e- -> XX process
};

class RSAXSpring;

//_______________________________________________________________________
// =====================
//  class RSAXSpringBuf
// =====================
//-----------------------------------------------------------------------
class RSAXSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  RSAXSpringBuf(const char *name   = "RSAXSpringBuf", 
                const char *title  = "RSAX Spring test event buffer",
	        RSAXSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~RSAXSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2AX     ()           const { return fQ2AX;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2X      ()           const { return fQ2X;       }
  Double_t GetCosThetaA()           const { return fCosThetaA; }
  Double_t GetPhiA     ()           const { return fPhiA;      }
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
  Double_t fQ2AX;           // q^2 of AX system
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fQ2X;            // q^2 of final state X
  Double_t fCosThetaA;      // cos(theta_a) in X   frame
  Double_t fPhiA;           // phi_a        in X   frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(RSAXSpringBuf, 1)  // RSAXSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class RSAXSpring
// =====================
//-----------------------------------------------------------------------
class RSAXSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   RSAXSpring(const char *name  = "RSAXSpring", 
              const char *title = "RSAX Spring test",
              RSAXBases  *bases = 0);
   virtual ~RSAXSpring();

   ClassDef(RSAXSpring, 1)  // RSAXSpring class
};

//_______________________________________________________________________
// =====================
//  class RSXBoson
// =====================
//-----------------------------------------------------------------------
class RSXBoson: public GENPDTEntry {
public:
   RSXBoson(Double_t m      = 120., 
            Double_t lambda = 1000.,
	    Double_t c0     = 0.,
	    Double_t c1     = 1.,
	    Double_t c2     = 1.,
	    Double_t c3     = 1.);
   virtual ~RSXBoson() {}

private:
   void     Initialize();
   Double_t GamToVV(Double_t m1, Double_t m2, Double_t a, Double_t color);

 private:
   Double_t fLambda;
   Double_t fC0;
   Double_t fC1;
   Double_t fC2;
   Double_t fC3;

   ClassDef(RSXBoson, 1)  // X boson class
};
#endif
