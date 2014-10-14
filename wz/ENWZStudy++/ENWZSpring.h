#ifndef ENWZSPRING_H
#define ENWZSPRING_H
//*****************************************************************************
//* =====================
//*  ENWZSpring
//* =====================
//*  
//* (Description)
//*    e+e- -> ENWZ generator
//*
//* (Update Record)
//*    2014/10/13  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include "HELLib.h"
#include "GENLib.h"

//_______________________________________________________________________
// =====================
//  class ENWZBases
// =====================
//-----------------------------------------------------------------------
class ENWZBases : public JSFBases {
friend class ENWZSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ENWZBases(const char *name  = "ENWZBases", 
            const char *title = "ENWZ Bases");
  virtual ~ENWZBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetWidth    ()           const { return fWidth;     }
  Int_t    GetCP       ()           const { return fCP;        }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2WZ     ()           const { return fQ2WZ;      }
  Double_t GetCosW     ()           const { return fCosW;      }
  Double_t GetPhiW     ()           const { return fPhiW;      }
  Double_t GetQ2W      ()           const { return fQ2W;       }
  Double_t GetCosWF    ()           const { return fCosWF;     }
  Double_t GetPhiWF    ()           const { return fPhiWF;     }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosZF    ()           const { return fCosZF;     }
  Double_t GetPhiZF    ()           const { return fPhiZF;     }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPolem    (Double_t p      ) { fPolem     = p;    }
  void     SetPolep    (Double_t p      ) { fPolep     = p;    }

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

  Double_t  DSigmaDX     (GENBranch &wwbranch);
  Double_t  AmpSquared   ();
  Complex_t FullAmplitude();
  Complex_t AmpEEtoENWZ  (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &emf,
                          const HELFermion &nbf,
                          const HELVector  &wf,
                          const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H
  Double_t fWidth;          // Gam_h  : width of H

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPolem;          // electron polarization
  Double_t fPolep;          // positron polarization
  Int_t    fWModesLo;       // W decay mode lo;
  Int_t    fWModesHi;       // W decay mode hi;
  Int_t    fZModesLo;       // Z decay mode lo;
  Int_t    fZModesHi;       // Z decay mode hi;
  Int_t    fNCALL;          // # points in a cell
  Double_t fACC1;           // accuracy for grid optimization
  Double_t fACC2;           // accuracy for integration
  Int_t    fITMX1;          // # iterations for grid optimization
  Int_t    fITMX2;          // # iterations for integration

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fXi;             // xi
  Double_t       fEta1;           // eta_1 for e- -> e-
  Double_t       fEta2;           // eta_2 for e+ -> e+
  Double_t       fPhi1;           // phi_1 for e- -> e-
  Double_t       fPhi2;           // phi_2 for e+ -> e+
  Double_t       fQ2WZ;           // q^2 of the WZ system
  Double_t       fCosW;           // cos(theta_W) in the WZ frame
  Double_t       fPhiW;           // phi_W        in the WZ frame
  Double_t       fQ2W;            // q^2 of the W system
  Double_t       fCosWF;          // cos(theta_fub) in the W frame
  Double_t       fPhiWF;          // phi_fub        in the W frame
  Double_t       fQ2Z;            // q^2 of the Z system
  Double_t       fCosZF;          // cos(theta_fu) in the Z frame
  Double_t       fPhiZF;          // phi_fu        in the Z frame
  GENDecayMode  *fWModePtr;       // pointer to W decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st W daughter (fub)
  GENPDTEntry   *f4Ptr;           // point to 2nd W daughter (fd)
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f5Ptr;           // point to 1st Z daughter (fu)
  GENPDTEntry   *f6Ptr;           // point to 2nd Z daughter (fdb)
  Int_t          fCP;             // CP flip: (-1,+1) = (nu e+ W- Z, e- nb W+ Z)
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [6];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fWMode;          // W decay mode
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[6];           // [0,1,2,3,4,5] = [e-, e+, fub, fd, fu, fdb]
  Double_t       fM[6];           // [0,1,2,3,4,5] = [m1, m2, m3, m4 , m5, m6 ]
  Double_t       fSh1;            // sh1
  Double_t       fCh1;            // ch1
  Double_t       fSh2;            // sh2
  Double_t       fCh2;            // ch2

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fXComb;          // CP & helicity combination
  Double_t       fWDecayMode;     // decay mode selector for W
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fXXi;            // xi
  Double_t       fXEta1;          // eta_1 for e- -> e-
  Double_t       fXEta2;          // eta_2 for e+ -> e+
  Double_t       fXPhi1;          // phi_1 for e- -> e-
  Double_t       fXPhi21;         // phi_2 - phi_1 for e+ -> e+
  Double_t       fXQ2WZ;          // q^2 of the WZ system
  Double_t       fXQ2W;           // q^2 of the W system
  Double_t       fXQ2Z;           // q^2 of the Z system

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1] = (coswm, phiwm, coswmf, 
  Double_t       fXL[6];          //!          phiwmf, coswpf, phiwpf)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(ENWZBases, 1) // Bases for e+e- -> XX process
};

class ENWZSpring;

//_______________________________________________________________________
// =====================
//  class ENWZSpringBuf
// =====================
//-----------------------------------------------------------------------
class ENWZSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ENWZSpringBuf(const char *name   = "ENWZSpringBuf", 
                const char *title  = "ENWZ Spring test event buffer",
	        ENWZSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ENWZSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Int_t    GetCP       ()           const { return fCP;        }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }
  Double_t GetXi       ()           const { return fXi;        }
  Double_t GetEta1     ()           const { return fEta1;      }
  Double_t GetEta2     ()           const { return fEta2;      }
  Double_t GetPhi1     ()           const { return fPhi1;      }
  Double_t GetPhi2     ()           const { return fPhi2;      }
  Double_t GetQ2WZ     ()           const { return fQ2WZ;      }
  Double_t GetCosW     ()           const { return fCosW;      }
  Double_t GetPhiW     ()           const { return fPhiW;      }
  Double_t GetQ2W      ()           const { return fQ2W;       }
  Double_t GetCosWF    ()           const { return fCosWF;     }
  Double_t GetPhiWF    ()           const { return fPhiWF;     }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosZF    ()           const { return fCosZF;     }
  Double_t GetPhiZF    ()           const { return fPhiZF;     }

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
  Int_t    fCP;             // CP flip: (-1,+1) = (nu e+ W- Z, e- nb W+ Z)
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung
  Double_t fXi;             // xi
  Double_t fEta1;           // eta_1 for e- -> e-
  Double_t fEta2;           // eta_2 for e+ -> e+
  Double_t fPhi1;           // phi_1 for e- -> e-
  Double_t fPhi2;           // phi_2 - phi_1 for e+ -> e+
  Double_t fQ2WZ;           // q^2 of the WZ system
  Double_t fCosW;           // cos(theta_W) in the WZ frame
  Double_t fPhiW;           // phi_W        in the WZ frame
  Double_t fQ2W;            // q^2 of the W system
  Double_t fCosWF;          // cos(theta_fub) in the W frame
  Double_t fPhiWF;          // phi_fub        in the W frame
  Double_t fQ2Z;            // q^2 of the Z system
  Double_t fCosZF;          // cos(theta_fu) in the Z frame
  Double_t fPhiZF;          // phi_fu        in the Z frame

  ClassDef(ENWZSpringBuf, 1)  // ENWZSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class ENWZSpring
// =====================
//-----------------------------------------------------------------------
class ENWZSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   ENWZSpring(const char *name  = "ENWZSpring", 
              const char *title = "ENWZ Spring test",
              ENWZBases  *bases = 0);
   virtual ~ENWZSpring();

   virtual Bool_t Initialize(); // overload the base class method

   ClassDef(ENWZSpring, 1)  // ENWZSpring class
};
#endif
