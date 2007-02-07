#ifndef __RSAXSPRING__
#define __RSAXSPRING__
//*****************************************************************************
//* =====================
//*  RSAXSpring
//* =====================
//*  
//* (Description)
//*    Kaluza-Klein  ee -> hh generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "TH1.h"

#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
const Double_t   k4Pi     = 4*TMath::Pi();    // GeV to fb
const Double_t   kGeV2fb  = 0.389379292e12;   // GeV to fb
const Double_t   kAlpha   = 1./128.;          // alpha(mz)  = 1/128.
const Double_t   kAlpha0  = 1./137.0359895;   // alpha(q=0) = 1/137.
const Double_t   kM_e     = 0.510998902e-3;   // electron mass [GeV]
const Double_t   kM_z     = 91.188;           // Z mass [GeV]
const Double_t   kSin2W   = 0.23;	              // sin^2(theta_W)
const Double_t   kSinW    = TMath::Sqrt(kSin2W);      // sin(theta_W)
const Double_t   kCos2W   = (1. - kSinW)*(1. + kSinW);// cos^2(theta_W)
const Double_t   kCosW    = TMath::Sqrt(kCos2W);      // cos(theta_W)
const Double_t   kSinCosW = kSinW*kCosW;              // sin(2theta_W)/2

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

  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetCosThetaA()           const { return fCosThetaA; }
  Double_t GetPhiA     ()           const { return fPhiA;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m)       { fMass    = m;      }
  void     SetEcmInit  (Double_t ecm)     { fEcmInit = ecm;    }
  void     SetISR      (Bool_t b = kTRUE) { fISR     = b;      }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr = b;      }

  // ----------------------
  //   Base class methods
  // ----------------------
#if 0
  virtual void     Userin();   // Bases user initialization
#endif
  virtual void     Userout();  // Bases user output 

  Double_t Func();     // Bases integration function.

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

  Double_t fMass;           // m_x mass of the higgs 

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?

  // ----------------
  //  Event info
  // ----------------
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaA;      // cos(theta_a) in X   frame
  Double_t       fPhiA;           // phi_a        in X   frame
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  TLorentzVector fP1;             // 1st daughter 4-momentum
  TLorentzVector fP2;             // 2nd daughter 4-momentum
  TLorentzVector fPA;             // gamma        4-memontum

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[4];          //! [0,1,2,3] = (cosx, phix, cosa, phia)
  Double_t       fXL[4];          //! 

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
  Double_t GetCosTheta ()           const { return fCosTheta;   }
  Double_t GetPhi      ()           const { return fPhi;        }
  Double_t GetCosThetaA()           const { return fCosThetaA;  }
  Double_t GetPhiA     ()           const { return fPhiA;       }
  Double_t GetZBoost   ()           const { return fZBoost;     }
  Double_t GetEcmIP    ()           const { return fEcmIP;      }

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
//  class GENFrame
// =====================
//-----------------------------------------------------------------------
class GENFrame {
public:
   GENFrame();
   GENFrame(const TLorentzVector &q, const GENFrame &eb);
   virtual ~GENFrame() {}

   TLorentzVector Transform(const TLorentzVector &pb);

private:
   TVector3 fEV[3];	// axis vectors

   ClassDef(GENFrame, 1)  // Reference frame class
};

//_______________________________________________________________________
// =====================
//  class GENPhase2
// =====================
//-----------------------------------------------------------------------
class GENPhase2 {
public:
   GENPhase2() {}
   GENPhase2(const TLorentzVector &q,
                   Double_t       m12,
                   Double_t       m22,
             const GENFrame      &eb,
                   Double_t       costh,
                   Double_t       phi,
                   Int_t          mode = 0);
   virtual ~GENPhase2() {}

   TLorentzVector GetFourMomentum(Int_t i);
   GENFrame       GetFrame(Int_t i = 1);
   Double_t       GetBetaBar();

private:
   Double_t       Beta2(Double_t x1, Double_t x2);
   void           Update();

private:
   TLorentzVector  fQ;		// parent 4-momentum
   Double_t        fM12;	// 1st dauter mass^2
   Double_t        fM22;	// 2nd dauter mass^2
   GENFrame        fEb;		// Eb: original frame
   GENFrame        fEa;		// Ea: parent's helicity frame
   Double_t        fCosTheta;   // cos(theta_1) in Ea
   Double_t        fPhi;        // phi_1        in Ea
   TLorentzVector  fP1;		// 1st daughter 4-momentum in Eb
   TLorentzVector  fP2;		// 2nd daughter 4-momentum in Eb
   Double_t        fBetaBar;	// beta_bar
   Int_t           fMode;	// (0,1)=(no transf, transf)
   Bool_t          fDone;	// true if updated

   ClassDef(GENPhase2, 1)  // 2-body phase space class
};


#endif
