#ifndef __EEFFSpring__
#define __EEFFSpring__

//////////////////////////////////////////////////////////////////////////
//
// EEFFSpring
//
// e+e- -> EEFFbar
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
//  COMMONs for EEFF calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;	// sqrt(s).
  Float_t  polebm;	// electron beam polarization.
  Float_t  sgmebm;	// beam energy spread (fraction).
  Int_t    isrb;	// SR and BM flag.
  Int_t	   igfr;	// gneration number.
  Int_t    itfr;	// isospin, (1,2)=(up,down).
  Int_t    lqfr;	// lepton-quark switch, (1,2)=(l,q).
  Float_t  pect;	// e- energy cut (GeV).
  Float_t  cecl;	// cos(theta_e-)_min
  Float_t  cecu;	// cos(theta_e-)_max
  Float_t  ppct;	// e+ energy cut (GeV).
  Float_t  cpcl;	// cos(theta_e+)_min
  Float_t  cpcu;	// cos(theta_e+)_max
  Float_t  pfct;	// fermion energy cut (GeV).
  Float_t  cfct;	// cos(theta_f) cut.
  Float_t  wmnf;	// minimum m_ff (GeV).
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class EEFFBases  : public JSFBases {
protected:
  Double_t fRoots;		// CM Energy(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;		// beam energy spread
  Int_t    fISRBM;		// Flag for ISR & Beamstrahlung
  Double_t fAlphai;		// 1/alpha(m_Z)
  Double_t fAlphas;		// alpha_s(m_Z)
  Double_t fMassW;		// m_W
  Double_t fMassZ;		// m_Z
  Double_t fMassHiggs;		// m_H
  Double_t fMassTop;		// m_t
  Int_t    fGeneration;		// fermion generation
  Int_t    fIsospin;		// isospin, (1,2)=(up,down).
  Int_t    fLorQ;		// lepton-quark switch, (1,2)=(l,q).
  Double_t fEeCut;		// e- energy cut (GeV).
  Double_t fCoseMin;		// cos(theta_e-)_min
  Double_t fCoseMax;		// cos(theta_e-)_max
  Double_t fEpCut;		// e+ energy cut (GeV).
  Double_t fCospMin;		// cos(theta_e+)_min
  Double_t fCospMax;		// cos(theta_e+)_max
  Double_t fEfCut;		// fermion energy cut (GeV).
  Double_t fCosfCut;		// cos(theta_f) cut.
  Double_t fMassffMin;		// minimum m_ff (GeV).
  Int_t    fISHUFL[50];  	// random number shuffler

public:
  EEFFBases(const char *name="EEFFBases", 
	     const char *title="EEFFbar  Bases");
  virtual ~EEFFBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
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
  
  void SetGeneration(Int_t i){ fGeneration=i; }
  Int_t GetGeneration(){ return fGeneration;}
  
  void SetIsospin(Int_t i){ fIsospin=i; }
  Int_t GetIsospin(){ return fIsospin;}
  
  void SetLorQ(Int_t i){ fLorQ=i; }
  Int_t GetLorQ(){ return fLorQ;}
  
  void SetEeCut(Double_t e){ fEeCut = e; }
  Double_t GetEeCut(){ return fEeCut;}
   
  void SetCoseMin(Double_t e){ fCoseMin = e; }
  Double_t GetCoseMin(){ return fCoseMin;}
    
  void SetCoseMax(Double_t e){ fCoseMax = e; }
  Double_t GetCoseMax(){ return fCoseMax;}

  void SetEpCut(Double_t e){ fEpCut = e; }
  Double_t GetEpCut(){ return fEpCut;}
   
  void SetCospMin(Double_t e){ fCospMin = e; }
  Double_t GetCospMin(){ return fCospMin;}
    
  void SetCospMax(Double_t e){ fCospMax = e; }
  Double_t GetCospMax(){ return fCospMax;}

  void SetEfCut(Double_t e){ fEfCut = e; }
  Double_t GetEfCut(){ return fEfCut;}
   
  void SetCosfCut(Double_t e){ fCosfCut = e; }
  Double_t GetCosfCut(){ return fCosfCut;}

  void SetMassffMin(Double_t e){ fMassffMin = e; }
  Double_t GetMassffMin(){ return fMassffMin;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(Double_t x[]); // Bases integration function.
  Double_t Func(); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(EEFFBases,1)  // Bases for e+e- -> EEFFbar process

};

class EEFFSpring;

class EEFFSpringBuf : public JSFSpringBuf {
public:
//  <<+EEFFbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-EEFFbar>>
public:
  EEFFSpringBuf(const char *name="EEFFSpringBuf", 
	     const char *title="EEFFSpring test event buffer",
	     EEFFSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEFFSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(EEFFSpringBuf,1)  // EEFFSpring event buffer
};

class EEFFSpring : public JSFSpring {
public:
   EEFFSpring(const char *name="EEFFSpring", 
	      const char *title="EEFFSpring test",
             EEFFBases *bases=NULL);
   virtual ~EEFFSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(EEFFSpring,1)  // EEFFSpring class
};


#endif










