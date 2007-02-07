#ifndef __RSZXSPRING__
#define __RSZXSPRING__
//*****************************************************************************
//* =====================
//*  RSZXSpring
//* =====================
//*  
//* (Description)
//*    RS+SUSY e+e- -> ZX generator
//*    Based on example/FFbarSpring directory.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "TNamed.h"

#include "ANL4DVector.h"
#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"
#include "JSFBeamGeneration.h"

#include <complex>
typedef std::complex<Double_t>  Complex_t;

class GENPDTEntry;
class GENPDTZBoson;
class GENBranch;
class HELFermion;
class HELScalar;
class HELVector;

//_______________________________________________________________________
// =====================
//  class RSZXBases
// =====================
//-----------------------------------------------------------------------
class RSZXBases : public JSFBases {
friend class RSZXSpringBuf;
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  RSZXBases(const char *name  = "RSZXBases", 
            const char *title = "RSZX Bases");
  virtual ~RSZXBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetQ2ZX     ()           const { return fQ2ZX;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2X      ()           const { return fQ2X;       }
  Double_t GetCosThetaA()           const { return fCosThetaA; }
  Double_t GetPhiA     ()           const { return fPhiA;      }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
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
  Complex_t AmpEEtoZX    (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &xf,
                          const HELVector  &zf);

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
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2ZX;           // q^2 of ZX system
  Double_t       fQ2X;            // q^2 of final state X
  Double_t       fQ2Z;            // q^2 of final state Z
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [ax1, ax2, fz1, fz2]
  Double_t       fM[4];           // [0,1,2,3] = [m1 , m2 , m3 , m4 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fCosThetaA;      // cos(theta_a) in X   frame
  Double_t       fPhiA;           // phi_a        in X   frame
  Double_t       fXQ2Z;           // q^2 of final state Z
  Double_t       fCosThetaF;      // cos(theta_f) in Z   frame
  Double_t       fPhiF;           // phi_f        in Z   frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  Double_t       fXU[6];          //! [0,1,2,3,4,5] = (cosx, phix,
  Double_t       fXL[6];          //!      cosa, phia, cosf, phif)

  TFile   *fBeamFile;         //! Beamstrahlung data file
  JSFBeamGenerationCain *fBM; //! Beamstrahlung data

  ClassDef(RSZXBases, 1) // Bases for e+e- -> XX process
};

class RSZXSpring;

//_______________________________________________________________________
// =====================
//  class RSZXSpringBuf
// =====================
//-----------------------------------------------------------------------
class RSZXSpringBuf : public JSFSpringBuf {
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  RSZXSpringBuf(const char *name   = "RSZXSpringBuf", 
                const char *title  = "RSZX Spring test event buffer",
	        RSZXSpring *spring = 0)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~RSZXSpringBuf() {}

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetQ2ZX     ()           const { return fQ2ZX;      }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2X      ()           const { return fQ2X;       }
  Double_t GetCosThetaA()           const { return fCosThetaA; }
  Double_t GetPhiA     ()           const { return fPhiA;      }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
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
  Double_t fQ2ZX;           // q^2 of ZX system
  Double_t fCosTheta;       // cos(theta_x) in cm  frame
  Double_t fPhi;            // phi_x        in cm  frame
  Double_t fQ2X;            // q^2 of final state X
  Double_t fCosThetaA;      // cos(theta_a) in X   frame
  Double_t fPhiA;           // phi_a        in X   frame
  Double_t fQ2Z;            // q^2 of final state Z
  Double_t fCosThetaF;      // cos(theta_f) in Z   frame
  Double_t fPhiF;           // phi_f        in Z   frame
  Double_t fZBoost;         // p_z(cm)      in lab frame
  Double_t fEcmIP;          // Ecm after B-strahlung

  ClassDef(RSZXSpringBuf, 1)  // RSZXSpring event buffer
};

//_______________________________________________________________________
// =====================
//  class RSZXSpring
// =====================
//-----------------------------------------------------------------------
class RSZXSpring : public JSFSpring {
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
public:
   RSZXSpring(const char *name  = "RSZXSpring", 
              const char *title = "RSZX Spring test",
              RSZXBases  *bases = 0);
   virtual ~RSZXSpring();

   ClassDef(RSZXSpring, 1)  // RSZXSpring class
};


//_______________________________________________________________________
// =====================
//  class TVectorC
// =====================
//-----------------------------------------------------------------------
class TVectorC : public TObject {
public:
#if 0
  TVectorC(Int_t n = 4) : fData(n) {}
#else
  TVectorC(Int_t n = 4) { fData[0]=0.; fData[1]=0.; fData[2]=0.; fData[3]=0.; }
#endif
  TVectorC(const TVectorC &src) : fData(src.fData) {}

  Complex_t  operator[](Int_t i) const { return fData[i]; }
  Complex_t &operator[](Int_t i)       { return fData[i]; }

private:
#if 0
  std::vector<Complex_t> fData; // data
#else
  Complex_t   fData[4]; // data
#endif

  ClassDef(TVectorC, 1) 	// Complex vector class
};

//_______________________________________________________________________
// =====================
//  class HELFermion
// =====================
//-----------------------------------------------------------------------
class HELFermion : public TVectorC {
friend class HELVector;
friend class HELScalar;
public:
#if 0
   HELFermion() : TVectorC(4) {}
#else
   HELFermion() {}
#endif
   HELFermion(const ANL4DVector &p,
                    Double_t     m,
                    Int_t        hel,
                    Int_t        nsf = 1,
		    Bool_t       isincom = kFALSE);
   virtual ~HELFermion() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Double_t     GetMass        () const  { return fM;   }
   inline       Int_t        GetHelicity    () const  { return fHel; }
   inline       Int_t        GetNSF         () const  { return fNSF; }
   
private:
   ANL4DVector fP;                // 4-momentum * fNSF
   Double_t    fM;                // mass
   Int_t       fHel;              // (-1,+1) = (-1/2,+1/2)
   Int_t       fNSF;              // (-1,+1) = (antiparticle, particle)
   Bool_t      fIsIncoming;       // true if incoming

   ClassDef(HELFermion, 1)  // Incoming fermion class
};

//_______________________________________________________________________
// =====================
//  class HELVector
// =====================
//-----------------------------------------------------------------------
class HELVector: public TVectorC {
friend class HELFermion;
friend class HELScalar;
public:
#if 0
   HELVector() : TVectorC(4) {}
#else
   HELVector() {}
#endif
   HELVector(const ANL4DVector &p,
                   Double_t     m,
	           Int_t        hel,
		   Int_t        nsv = 1);
   HELVector(const HELFermion &fin, 
             const HELFermion &fout,
	           Double_t    gl,
		   Double_t    gr,
		   Double_t    m,
		   Double_t    gm);
   virtual ~HELVector() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Double_t     GetMass        () const  { return fM;   }
   inline       Int_t        GetHelicity    () const  { return fHel; }
   inline       Int_t        GetNSV         () const  { return fNSV; }
   
private:
   ANL4DVector fP;                // 4-momentum * fNSV
   Double_t    fM;                // mass
   Int_t       fHel;              // helicity
   Int_t       fNSV;              // (-1,1) = (initial, final)

   ClassDef(HELVector, 1)  // Vector boson class
};

//_______________________________________________________________________
// =====================
//  class HELScalar
// =====================
//-----------------------------------------------------------------------
class HELScalar: public Complex_t {
friend class HELFermion;
friend class HELVector;
public:
   HELScalar() {}
   HELScalar(const ANL4DVector &p,
		   Int_t        nss = 1);
   virtual ~HELScalar() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Int_t        GetNSS         () const  { return fNSS; }
   
private:
   ANL4DVector fP;                // 4-momentum * fNSS
   Int_t       fNSS;              // (-1,1) = (initial, final)

   ClassDef(HELScalar, 1)  // Scalar boson class
};

//_______________________________________________________________________
// =====================
//  class GENDecayMode
// =====================
//-----------------------------------------------------------------------
class GENDecayMode : public TObjArray {
friend class GENModePicker;
public:
   GENDecayMode(Double_t gm = 0.) : fGamma(gm), fBR(0.), fCumBR(0.) {}
   virtual ~GENDecayMode() {}

   Double_t GetGamma()             { return fGamma; }
   Double_t GetBR   ()             { return fBR;    }
   void     SetGamma(Double_t gm ) { fGamma = gm;   }
   void     SetBR   (Double_t br ) { fBR    = br;   }

   void     DebugPrint(const Option_t *opt ="");

private:
   Double_t fGamma;		// partial width
   Double_t fBR;		// branching fraction
   Double_t fCumBR;		// cumulative branching fraction

   ClassDef(GENDecayMode, 1) 	// Decay mode class
};

//_______________________________________________________________________
// =====================
//  class GENModePicker
// =====================
//-----------------------------------------------------------------------
class GENModePicker : public TObjArray {
public:
   GENModePicker() : fGamma(0.), fDone(kFALSE) {}
   virtual ~GENModePicker() {}

   using TObjArray::Add;
   virtual void     Add          (GENDecayMode *mp);
           GENDecayMode *PickMode(Double_t x, 
			          Double_t &weight,
				  Int_t &mode);

           Double_t GetWidth  () { if(!fDone) Update(); return fGamma; }

protected:
   virtual void     Update();

private:
   Double_t fGamma;		// total width [GeV]
   Bool_t   fDone;		// true if updated

   ClassDef(GENModePicker, 1) 	// Decay mode picker class
};

//_______________________________________________________________________
// =====================
//  class GENPDTEntry
// =====================
//-----------------------------------------------------------------------
class GENPDTEntry: public GENModePicker {
public:
   GENPDTEntry() {}
   GENPDTEntry(const Char_t     *name,
                     Int_t       pid,
                     Double_t    charge,
                     Double_t    spin,
                     Double_t    mass,
                     Int_t       gen   = 0,
                     Double_t    ispin = 0.,
                     Double_t    color = 1.);
   virtual ~GENPDTEntry() {}

   TString & GetName  () { return fName;    }
   Int_t     GetPID   () { return fPID;     }
   Double_t  GetCharge() { return fCharge;  }
   Double_t  GetMass  () { return fMass;    }
   Int_t     GetGenNo () { return fGen;     }
   Double_t  GetISpin () { return fIsoSpin; }
   Double_t  GetColor () { return fColor;   }

   Double_t  GetQ2BW  (Double_t    qmin,   // Q_min
                       Double_t    qmax,   // Q_max
                       Double_t       x,   // integration variable
                       Double_t &weight);  // Jacobian weight

   void      DebugPrint(const Option_t *opt = "");

protected:
   TString       fName;		//  name
   Int_t         fPID;		//  PDG ID code
   Double_t      fCharge;	//  charge
   Double_t      fSpin;		//  spin
   Double_t      fMass;		//  mass  [GeV]
   Int_t         fGen;		//  generation
   Double_t      fIsoSpin;	//  (0, 1, 2) = (0, up, down)
   Double_t      fColor;	//  color factor

   ClassDef(GENPDTEntry, 1) 	// PD table entry class
};

//_______________________________________________________________________
// =====================
//  class GENPDTZBoson
// =====================
//-----------------------------------------------------------------------
class GENPDTZBoson: public GENPDTEntry {
public:
   GENPDTZBoson();
   ~GENPDTZBoson();

private:
   void     Initialize();   
   Double_t GamToFF(Double_t t3, Double_t qf, Double_t cf, Double_t m);

   ClassDef(GENPDTZBoson, 1)  // Reference frame class
};

//_______________________________________________________________________
// =====================
//  class GENBranch
// =====================
//-----------------------------------------------------------------------
class GENBranch {
public:
   GENBranch(Double_t q2    = 0.,
             Double_t costh = 0.,
             Double_t phi   = 0.,
             Double_t m12   = 0.,
	     Double_t m22   = 0.);

   GENBranch(Double_t   q2,
             Double_t   costh,
             Double_t   phi,
             GENBranch *br1p,
             GENBranch *br2p);

   virtual ~GENBranch() {}

   inline Double_t GetQ2      ()            { return fQ2;       }
   inline Double_t GetCosTheta()            { return fCosTheta; }
   inline Double_t GetPhi     ()            { return fPhi;      }
   inline Double_t GetM12     ()            { return fM12;      }
   inline Double_t GetM22     ()            { return fM22;      }
   inline Double_t GetBetaBar ()            { return fBetaBar;  }

   inline void     SetQ2      (Double_t q2) { fQ2       = q2;   }
   inline void     SetCosTheta(Double_t cs) { fCosTheta = cs;   }
   inline void     SetPhi     (Double_t fi) { fPhi      = fi;   }

   inline GENBranch * GetBranchPtr(Int_t i) { return i ? fBR2Ptr
                                                       : fBR1Ptr; }

private:
   Double_t   fQ2;        // q^2
   Double_t   fCosTheta;  // cos(theta)
   Double_t   fPhi;       // phi
   Double_t   fM12;       // m1*m1
   Double_t   fM22;       // m2*m2
   GENBranch *fBR1Ptr;    // 1st daughter branch if any
   GENBranch *fBR2Ptr;    // 2nd daughter branch if any
   Double_t   fBetaBar;   // beta_bar

   ClassDef(GENBranch, 1)  // Branch class
};

//_______________________________________________________________________
// =====================
//  class GENFrame
// =====================
//-----------------------------------------------------------------------
class GENFrame {
public:
   GENFrame();
   GENFrame(const ANL4DVector &q, const GENFrame &eb);
   virtual ~GENFrame() {}

   ANL4DVector Transform(const ANL4DVector &pb);

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
   GENPhase2(const ANL4DVector &q,
                   Double_t     m12,
                   Double_t     m22,
             const GENFrame    &eb,
                   Double_t     costh,
                   Double_t     phi,
                   Int_t        mode = 0);
   virtual ~GENPhase2() {}

   ANL4DVector GetFourMomentum(Int_t i);
   GENFrame    GetFrame(Int_t i = 1);
   Double_t    GetBetaBar();

private:
   Double_t       Beta2(Double_t x1, Double_t x2);
   void           Update();

private:
   ANL4DVector  fQ;		// parent 4-momentum
   Double_t     fM12;		// 1st dauter mass^2
   Double_t     fM22;		// 2nd dauter mass^2
   GENFrame     fEb;		// Eb: original frame
   GENFrame     fEa;		// Ea: parent's helicity frame
   Double_t     fCosTheta;   	// cos(theta_1) in Ea
   Double_t     fPhi;        	// phi_1        in Ea
   ANL4DVector  fP1;		// 1st daughter 4-momentum in Eb
   ANL4DVector  fP2;		// 2nd daughter 4-momentum in Eb
   Double_t     fBetaBar;	// beta_bar
   Int_t        fMode;		// (0,1)=(no transf, transf)
   Bool_t       fDone;		// true if updated

   ClassDef(GENPhase2, 1)  // 2-body phase space class
};


#endif
