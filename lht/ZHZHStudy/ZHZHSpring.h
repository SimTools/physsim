#ifndef ZHZHSPRING_H
#define ZHZHSPRING_H
//*****************************************************************************
//* =====================
//*  ZHZHSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ZH ZH generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2008/12/28  K.Fujii	Original version.
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
//  class ZHZHBases
// =====================
//-----------------------------------------------------------------------
class ZHZHBases : public JSFBases {
friend class ZHZHSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHZHBases(const char *name  = "ZHZHBases", 
            const char *title = "ZHZH Bases");
  virtual ~ZHZHBases();

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

  Double_t GetQ2XX      ()           const { return fQ2XX;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X1      ()           const { return fQ2X1;       }
  Double_t GetCosThetaH1()           const { return fCosThetaH1; }
  Double_t GetPhiH1     ()           const { return fPhiH1;      }
  Double_t GetQ2X2      ()           const { return fQ2X2;       }
  Double_t GetCosThetaH2()           const { return fCosThetaH2; }
  Double_t GetPhiH2     ()           const { return fPhiH2;      }
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
                          const HELVector  &x1,
                          const HELVector  &x2);

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
  Int_t    fIBType;            // Beam spread type (1: Uniform, 2: Gauss)
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam Spread relative to nominal Ebm
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
  Double_t       fQ2XX;           // q^2 of XX system
  Double_t       fQ2X1;           // q^2 of final state X1
  Double_t       fQ2X2;           // q^2 of final state X2
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [h1, ah1, h2, ah2]
  Double_t       fM[3];           // [0,1,2,3] = [mh, mdm, mh, mdm]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fCosTheta;       // cos(theta_x1) in cm frame
  Double_t       fPhi;            // phi_x1        in cm frame
  Double_t       fCosThetaH1;     // cos(theta_h1) in X1 frame
  Double_t       fPhiH1;          // phi_h1        in X1 frame
  Double_t       fXQ2X1;          // q^2 of final state X1
  Double_t       fCosThetaH2;     // cos(theta_h1) in X2 frame
  Double_t       fPhiH2;          // phi_h1        in X1 frame
  Double_t       fXQ2X2;          // q^2 of final state X2

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1,2,3,4,5] = (cosx, phix,
  Double_t       fXL[6];          //!    cosh1, phih1, cosh2, phih2)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ZHZHBases, 1) // Bases for e+e- -> XX process
};

class ZHZHSpring;

//_______________________________________________________________________
// =====================
//  class ZHZHSpringBuf
// =====================
//-----------------------------------------------------------------------
class ZHZHSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHZHSpringBuf(const char   *name   = "ZHZHSpringBuf", 
                const char   *title  = "ZHZH Spring test event buffer",
	        ZHZHSpring   *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZHZHSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2XX      ()           const { return fQ2XX;       }
  Double_t GetCosTheta  ()           const { return fCosTheta;   }
  Double_t GetPhi       ()           const { return fPhi;        }
  Double_t GetQ2X1      ()           const { return fQ2X1;       }
  Double_t GetCosThetaH1()           const { return fCosThetaH1; }
  Double_t GetPhiH1     ()           const { return fPhiH1;      }
  Double_t GetQ2X2      ()           const { return fQ2X2;       }
  Double_t GetCosThetaH2()           const { return fCosThetaH2; }
  Double_t GetPhiH2     ()           const { return fPhiH2;      }
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
  Double_t fCosThetaH1;      // cos(theta_h1) in X1 frame
  Double_t fPhiH1;           // phi_h1        in X1 frame
  Double_t fQ2X2;            // q^2 of final state X2
  Double_t fCosThetaH2;      // cos(theta_h1) in X2 frame
  Double_t fPhiH2;           // phi_h1        in X2 frame
  Double_t fZBoost;          // p_z(cm)       in lab frame
  Double_t fEcmIP;           // Ecm after B-strahlung

  ClassDef(ZHZHSpringBuf, 1)  // ZHZHSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ZHZHSpring
// =====================
//-----------------------------------------------------------------------
class ZHZHSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ZHZHSpring(const char   *name  = "ZHZHSpring", 
              const char   *title = "ZHZH Spring test",
              ZHZHBases    *bases = 0);
   virtual ~ZHZHSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ZHZHSpring, 1)  // ZHZHSpring class
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
