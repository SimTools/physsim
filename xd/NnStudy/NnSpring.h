#ifndef NNSPRING_H
#define NNSPRING_H
//*****************************************************************************
//* =====================
//*  NnSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> Nn generator
//*
//* (Update Record)
//*    2009/08/11  T.Saito	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

class RNeutrino;

//_______________________________________________________________________
// =====================
//  class NnBases
// =====================
//-----------------------------------------------------------------------
class NnBases : public JSFBases {
friend class NnSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NnBases(const char *name  = "NnBases", 
          const char *title = "NnBases");
  virtual ~NnBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetMass4    ()           const { return fMass4;     }
  Double_t GetMassnu   ()           const { return fMassnu;    }
  Int_t    GetNkk      ()           const { return fNkk;       }
  Int_t    GetGenNu    ()           const { return fGenNu;     }
  Int_t    GetGenLepton()           const { return fGenLepton; }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPole     (Double_t p      ) { fPole      = p;    }

  Double_t GetQ2XX     ()           const { return fQ2XX;      }
  Double_t GetQ2X1     ()           const { return fQ2X1;      }
  Double_t GetQ2W1     ()           const { return fQ2W1;      }

  Double_t GetCosThetaW1()     const { return fCosThetaW1;     }
  Double_t GetPhiW1     ()     const { return fPhiW1;          }

  Double_t GetCosThetaF1()     const { return fCosThetaF1;     }
  Double_t GetPhiF1     ()     const { return fPhiF1;          }


  // ----------------------
  //   Base class methods
  // ----------------------
  virtual void     Userin( );  // Bases user initialization
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
  Complex_t AmpEEtoNn    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &fnr,
			  const HELFermion &fn,
			        Double_t glznrn,
			        Double_t grznrn,
			        Double_t glrnwe,
			        Double_t grrnwe,
			        Double_t glnwe,
			        Double_t grnwe,
		  	        Double_t mrnu,
                                Double_t gamrn);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization

  Double_t fMass;           // right handed neutrino mass
  Double_t fMassnu;         // neutrino mass
  Double_t fMass4;          // 4-dim neutrino mass
  Int_t    fNkk;            // NR KK mode number
  Int_t    fGenNu;          // NR Generation
  Int_t    fGenLepton;      // Lepton Generation

  Int_t    fWpModesLo;      // W+ decay mode lo;
  Int_t    fWpModesHi;      // W+ decay mode hi;


  // ----------------
  //  Particle Data
  // ----------------
  GENPDTEntry  *fFPtr;            //! PD table entry of "f"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  RNeutrino    *fNRPtr;           //! PD table entry of "NR[i]"
  GENPDTWBoson *fW1BosonPtr;      //! PD table entry of "W1"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1] = [f, fbar]

  Double_t       fQ2XX;           // q^2 of XX system
  Double_t       fQ2X1;           // q^2 of final state X1
  Double_t       fQ2W1;           // q^2 of final state W+

  GENDecayMode  *fW1ModePtr;      // pointer to W+ decay mode
  Int_t          fW1Mode;         // W+ decay mode
  GENDecayMode  *fZ1ModePtr;      // pointer to Z decay mode
  Int_t          fZ1Mode;         // Z decay mode
  GENPDTEntry   *f1Ptr;           // point to 1st W+ daughter
  GENPDTEntry   *f2Ptr;           // point to 2nd W+ daughter

  Double_t       fM[4];           // [0,1,2,3] = [nub, lep, up, downb]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fW1DecayMode;    // decay mode selector for W+

  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1] = (cosx, phix)
  Double_t       fXL[6];          //!

  Double_t       fCosThetaW1;     // cos(theta_W+) in X1  frame
  Double_t       fPhiW1;          // phi_W+        in X1  frame
  Double_t       fXQ2X1;          // q^2 of final state X1
  Double_t       fCosThetaF1;     // cos(theta_f1) in W+  frame
  Double_t       fPhiF1;          // phi_f1        in W+  frame
  Double_t       fXQ2W1;          // q^2 of final state W+


  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(NnBases, 1) // Bases for e+e- -> XX process
};

class NnSpring;

//_______________________________________________________________________
// =====================
//  class NnSpringBuf
// =====================
//-----------------------------------------------------------------------
class NnSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  NnSpringBuf(const char *name   = "NnSpringBuf", 
              const char *title  = "NnSpring test event buffer",
	        NnSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NnSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetCosTheta  ()           const { return fCosTheta;  }
  Double_t GetPhi       ()           const { return fPhi;       }
  Double_t GetZBoost    ()           const { return fZBoost;    }
  Double_t GetEcmIP     ()           const { return fEcmIP;     }

  Double_t GetQ2XX      ()           const { return fQ2XX;      }
  Double_t GetQ2X1      ()           const { return fQ2X1;      }
  Double_t GetQ2W1      ()           const { return fQ2W1;      }

  Double_t GetCosThetaW1    ()       const { return fCosThetaW1;     }
  Double_t GetPhiW1         ()       const { return fPhiW1;          }
  Double_t GetCosThetaF1    ()       const { return fCosThetaF1;     }
  Double_t GetPhiF1         ()       const { return fPhiF1;          }

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
  Double_t fCosTheta;       // cos(theta_x2) in cm  frame
  Double_t fPhi;            // phi_x2        in cm  frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  Double_t       fCosThetaW1;     // cos(theta_W+) in X1  frame
  Double_t       fPhiW1;          // phi_W+        in X1  frame
  Double_t       fXQ2X1;          // q^2 of final state X1
  Double_t       fCosThetaF1;     // cos(theta_f1) in W+  frame
  Double_t       fPhiF1;          // phi_f1        in W+  frame
  Double_t       fXQ2W1;          // q^2 of final state W+

  Double_t       fQ2XX;           // q^2 of XX system
  Double_t       fQ2X1;           // q^2 of final state X1
  Double_t       fQ2W1;           // q^2 of final state W+

  ClassDef(NnSpringBuf, 1)  // NnSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class NnSpring
// =====================
//-----------------------------------------------------------------------
class NnSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   NnSpring(const char *name  = "NnSpring", 
            const char *title = "NnSpring test",
              NnBases  *bases = 0);
   virtual ~NnSpring();

   ClassDef(NnSpring, 1)  // NnSpring class
};

//_______________________________________________________________________
// =====================
//  class RNeutrino
// =====================
//-----------------------------------------------------------------------
class RNeutrino: public GENPDTEntry {
 public:

  RNeutrino(Double_t m      = 100.,
	    Double_t m4     = 3.3e-9,
	    Int_t    gen    = 1,
	    Int_t    kkmode = 1);
  virtual ~RNeutrino() {}

        Int_t     GetNkk()            const { return fN;            }
  const Double_t *GetGwl()            const { return &fGwl[0];      }
  const Double_t *GetGwe()            const { return &fGwe[0];      }
  const Double_t *GetGzn(Int_t j = 1) const { return &fGzn[0][j-1]; }

 private:
  void     Initialize();
  Double_t GamToFV(Double_t m1, Double_t m2, Double_t a);

 private:
  Double_t fMass4;     // 4-dim NR mass                                  
  Int_t    fGen;       // N2(basis lepton's geneation)
  Int_t    fN;         // KK mode number
  Double_t fGwl[2];    // (gL,gR)_W+NR_2l-l
  Double_t fGwe[2];    // (gL,gR)_W+NR_2l-e   
  Double_t fGzn[2][3]; // (gL,gR)_Z-NR_2l-n_j (j=e,mu,tau)

  ClassDef(RNeutrino, 1)  // NR boson class
};

#endif
