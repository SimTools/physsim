#ifndef NNZHSPRING_H
#define NNZHSPRING_H

//////////////////////////////////////////////////////////////////////////
//
// NNZHSpring
//
// e+e- -> nu nu Z H
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
//  COMMONs for NNZH calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;	// sqrt(s).
  Float_t  polebm;	// electron beam polarization.
  Float_t  sgmebm;	// beam energy spread (fraction).
  Int_t    isrb;	// ISR and BM flag.
  Int_t    imdmin;	// 1st Z0 decay mode ID.
  Int_t    imdmax;	// lst Z0 decay mode ID.
#ifdef WITH_DBD_STANDARD
  Int_t    isr_lla_order; // LLA_Order for ISR calculation
  Int_t    ns_remnants;   // (0,1)=(not store, store) ISR remnants
  Int_t    ns_beams;     // (0,1)=(not store, store) e+/e- after beamstrahlung
#endif
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class NNZHBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fZModesLo;	// 1st Z0 decay mode ID.
  Int_t    fZModesHi;	// lst Z0 decay mode ID.
#ifdef WITH_DBD_STANDARD
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
  NNZHBases(const char *name="NNZHBases", 
	     const char *title="NNZHbar  Bases");
  virtual ~NNZHBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
   
  void SetZModesLo(Int_t i){ fZModesLo=i; }
  Int_t GetZModesLo(){ return fZModesLo;}
   
  void SetZModesHi(Int_t i){ fZModesHi=i; }
  Int_t GetZModesHi(){ return fZModesHi;}
   
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

  ClassDef(NNZHBases,1)  // Bases for e+e- -> NNZH process

};

class NNZHSpring;

class NNZHSpringBuf : public JSFSpringBuf {
public:
//  <<+NNZHbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-NNZHbar>>
public:
  NNZHSpringBuf(const char *name="NNZHSpringBuf", 
	     const char *title="NNZHSpring test event buffer",
	     NNZHSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~NNZHSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(NNZHSpringBuf,1)  // NNZHSpring event buffer
};

class NNZHSpring : public JSFSpring {
public:
   NNZHSpring(const char *name="NNZHSpring", 
	      const char *title="NNZHSpring test",
             NNZHBases *bases=NULL);
   virtual ~NNZHSpring();

   virtual Bool_t Initialize(); // 
   virtual void  GetPy6frmProb(int nseq, double prob[7]); 

   ClassDef(NNZHSpring,1)  // NNZHSpring class
};


#endif
