#ifndef __TTHSpring__
#define __TTHSpring__

//////////////////////////////////////////////////////////////////////////
//
// TTHSpring
//
// e+e- -> ttbar H
//
//////////////////////////////////////////////////////////////////////////

#include "JSFConfig.h"

#include "TNamed.h"
#include "TMath.h"
#include "TDatime.h"

#ifndef __JSFModule__
#include "JSFModule.h"
#endif
#ifndef __JSFBases__
#include "JSFBases.h"
#endif
#ifndef __JSFSpring__
#include "JSFSpring.h"
#endif


// =====================================================================
//  COMMONs for shuffling random numbers.
// =====================================================================
typedef struct {
  Int_t nzz;
  Int_t ishufl[50];
} COMMON_BSHUFL;             //  Common /BSHUFL/

extern COMMON_BSHUFL bshufl_;

// =====================================================================
//  COMMONs for SM parameters, see USMPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  alfi, alfs, amsw, amsz, amsh, amst;
} COMMON_USMPRM;             //  Common /USMPRM/

extern COMMON_USMPRM usmprm_;

// =====================================================================
//  COMMONs for TTH calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts, polebm, sgmebm;
  Int_t    isrb;
  Int_t    imd1lo;	// W- decay mode ID.
  Int_t    imd1hi;	// W- decay mode ID.
  Int_t    imd2lo;	// W+ decay mode ID.
  Int_t    imd2hi;	// W+ decay mode ID.
  Float_t  polpbm;      // positron polarization.
#ifdef WITH_DBD_STANDARD
  Int_t    ndcysl;        // W+W- decay mode selector
  Int_t    isr_lla_order; // LLA_Order for ISR calculation
  Int_t    ns_remnants;   // (0,1)=(not store, store) ISR remnants
  Int_t    ns_beams;     // (0,1)=(not store, store) e+/e- after beamstrahlung
#endif
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

typedef struct {
  Float_t roots;
  Float_t pole;
  Float_t xphase[11][3];
  Int_t   idptcl[14][3];
  Int_t   ihlcmb[14];
  Float_t pvout[14][4];
  Float_t dsgmdx;
  Float_t dpdebm;
  Float_t bmengy[5];
  Int_t   isrbm;
  Int_t   imd1mn;
  Int_t   imd1mx;
  Int_t   imd2mn;
  Int_t   imd2mx;
  Float_t polp;
#ifdef WITH_DBD_STANDARD
  Double_t z3;
  Double_t z4;
#endif

} TTZCNF_t;
extern TTZCNF_t ttzcnf_;

class TTHBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy (GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fPolPositron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fWmModesLo;	// W- decay mode ID.
  Int_t    fWmModesHi;	// W- decay mode ID.
  Int_t    fWpModesLo;	// W+ decay mode ID.
  Int_t    fWpModesHi;	// W+ decay mode ID.
#ifdef WITH_DBD_STANDARD
  Int_t    fFinalStatesMix; // If not zero, W- decay mode means one of W 
			    //  and W; decay mode means another W.  This is useful 
	                    // to generate (W-,W+)--> 
			   //	(lnu+qq) and (qq+ln) simultaneously 
  Int_t    fISR_LLA_Order; // ISR LLA order
  Int_t    fStoreRemnants; // Store remnants or not (0,1)=(no, yes)
  Int_t    fStoreBeams;    // Store e+e- beams after Beamstrahlung (0,1)=(no,yes)
  std::string  fLumiFileDirectory; // directory of lumi file.
#endif
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  TTHBases(const char *name="TTHBases", 
	     const char *title="TTHbar  Bases");
  virtual ~TTHBases(){} 

  void SetRoots(Double_t sqrts){ fRoots=sqrts; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}

  void SetPolPositron(Double_t polpbm){ fPolPositron=polpbm; }
  Double_t GetPolPositron(){ return fPolPositron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
  
  void SetAlphai(Double_t alphai){ fAlphai=alphai; }
  Double_t GetAlphai(){ return fAlphai;}
   
  void SetAlphas(Double_t alphas){ fAlphas=alphas; }
  Double_t GetAlphas(){ return fAlphas;}
   
  void SetMassW(Double_t mass){ fMassW=mass; }
  Double_t GetMassW(){ return fMassW;}
   
  void SetMassZ(Double_t mass){ fMassZ=mass; }
  Double_t GetMassZ(){ return fMassZ;}
   
  void SetMassHiggs(Double_t mass){ fMassHiggs=mass; }
  Double_t GetMassHiggs(){ return fMassHiggs;}
   
  void SetMassTop(Double_t mass){ fMassTop=mass; }
  Double_t GetMassTop(){ return fMassTop;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(); // New style not yet
  Double_t Func(Double_t x[]); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(TTHBases,1)  // Bases for e+e- -> ttbar process

};

class TTHSpring;

class TTHSpringBuf : public JSFSpringBuf {
public:
//  <<+TTHbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-TTHbar>>
public:
  TTHSpringBuf(const char *name="TTHSpringBuf", 
	     const char *title="TTHSpring test event buffer",
	     TTHSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~TTHSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(TTHSpringBuf,1)  // TTHSpring event buffer
};

class TTHSpring : public JSFSpring {
public:
   TTHSpring(const char *name="TTHSpring", 
	      const char *title="TTHSpring test",
             TTHBases *bases=NULL);
   virtual ~TTHSpring();

   virtual Bool_t Initialize(); // 
   virtual void  GetPy6frmProb(int nseq, double prob[7]); 

   ClassDef(TTHSpring,2)  // TTHSpring class
};


#endif







