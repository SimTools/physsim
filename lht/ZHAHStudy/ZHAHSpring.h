#ifndef ZHAHSPRING_H
#define ZHAHSPRING_H
//*****************************************************************************
//* =====================
//*  ZHAHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ZH AH generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/11/29  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class ZHBoson;
class AHBoson;

//_______________________________________________________________________
// =====================
//  class ZHAHBases
// =====================
//-----------------------------------------------------------------------
class ZHAHBases : public JSFBases {
friend class ZHAHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHAHBases(const char *name  = "ZHAHBases", 
            const char *title = "ZHAH Bases");
  virtual ~ZHAHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetF         ()           const { return fF;          }
  Double_t GetKappaL    ()           const { return fKappaL;     }
  Double_t GetMass      ()           const { return fMass;       }
  Double_t GetMassDM    ()           const { return fMassDM;     }
  Double_t GetMassT     ()           const { return fMassT;      }
  Double_t GetMassHiggs ()           const { return fMassHiggs;  }
  Double_t GetEcmInit   ()           const { return fEcmInit;    }

  Double_t GetQ2XD      ()           const { return fQ2XD;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X       ()           const { return fQ2X;        }
  Double_t GetCosThetaH ()           const { return fCosThetaH;  }
  Double_t GetPhiH      ()           const { return fPhiH;       }
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
  Complex_t AmpEEtoXD    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELVector  &x,
                          const HELVector  &dm);

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
  Double_t fMassHiggs;      // m_dm   : mass  of DM

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fPole;           // electron polarization

  // ----------------
  //  Particle Data
  // ----------------
  ZHBoson      *fXBosonPtr;       //! PD table entry of "X"
  AHBoson      *fDMBosonPtr;      //! PD table entry of "DM"
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
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [3];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[3];           // [0,1,2] = [h , ah1, ah2]
  Double_t       fM[3];           // [0,1,2] = [mh, mdm, mdm]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaH;      // cos(theta_h) in X  frame
  Double_t       fPhiH;           // phi_h        in X  frame
  Double_t       fXQ2X;           // q^2 of final state X

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[4];          //! [0,1,2,3] = (cosx, phix,
  Double_t       fXL[4];          //!              cosh, phih)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ZHAHBases, 1) // Bases for e+e- -> XD process
};

class ZHAHSpring;

//_______________________________________________________________________
// =====================
//  class ZHAHSpringBuf
// =====================
//-----------------------------------------------------------------------
class ZHAHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHAHSpringBuf(const char   *name   = "ZHAHSpringBuf", 
                const char   *title  = "ZHAH Spring test event buffer",
	          ZHAHSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZHAHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2XD      ()           const { return fQ2XD;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X       ()           const { return fQ2X;        }
  Double_t GetCosThetaH ()           const { return fCosThetaH;  }
  Double_t GetPhiH      ()           const { return fPhiH;       }
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
  Double_t fCosThetaH;       // cos(theta_h) in X  frame
  Double_t fPhiH;            // phi_h        in X  frame
  Double_t fZBoost;          // p_z(cm)      in lab frame
  Double_t fEcmIP;           // Ecm after B-strahlung

  ClassDef(ZHAHSpringBuf, 1)  // ZHAHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ZHAHSpring
// =====================
//-----------------------------------------------------------------------
class ZHAHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ZHAHSpring(const char   *name  = "ZHAHSpring", 
              const char   *title = "ZHAH Spring test",
                ZHAHBases  *bases = 0);
   virtual ~ZHAHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ZHAHSpring, 1)  // ZHAHSpring class
};

//_______________________________________________________________________
// =====================
//  class ZHBoson
// =====================
//-----------------------------------------------------------------------

class ZHBoson: public GENPDTEntry {
public:
   ZHBoson(Double_t f  = 580.,
           Double_t mh = 134.);
   virtual ~ZHBoson() {}

   Double_t  GetF     () const { return fF;      }
   Double_t  GetSh    () const { return fSh;     }
   Double_t  GetCh    () const { return fCh;     }
   Double_t  GetSf    () const { return fSf;     }
   Double_t  GetCf    () const { return fCf;     }
   Double_t  GetMassDM() const { return fMassDM; }

private:
   void     Initialize();
   Double_t GamToSV(Double_t m1, Double_t m2, Double_t a);

private:
   Double_t fF;         // SU(5)->SO(5) breaking scale
   Double_t fSh;        // Sh parameter
   Double_t fCh;        // Ch parameter
   Double_t fSf;        // Sf parameter
   Double_t fCf;        // Cf parameter
   Double_t fMassDM;    // daughter DM mass
   Double_t fMassHiggs; // daughter DM mass

   ClassDef(ZHBoson, 1)  // X boson class
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
