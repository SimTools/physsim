//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  TTZSpring
//  
//  e+e- -> ttbar Z
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : ISR
//     ( 4) : e- helicity
//     ( 5) : helicity combination for final states.
//     ( 6) : m(t-bar)**2
//     ( 7) : m(t)**2
//     ( 8) : m(J)**2
//     ( 9) : m(t-t_bar)**2
//     (10) : m(W-)**2
//     (11) : m(W+)**2
//     (12) : cos(theta_H)
//     (13) : phi_H
//     (14) : cos(theta_t)     in t-t_bar rest frame
//     (15) : phi_t            in t-t_bar rest frame
//     (16) : cos(theta_q_bar) in t_bar rest frame
//     (17) : phi_q_bar        in t_bar rest frame
//     (18) : cos(theta_f)     in W- rest frame
//     (19) : phi_f            in W- rest frame
//     (20) : cos(theta_q)     in t rest frame
//     (21) : phi_q            in t rest frame
//     (22) : cos(theta_f_bar) in W+ rest frame
//     (23) : phi_f_bar        in W+ rest frame
//     (24) : cos(theta_f_bar) in Z rest frame
//     (25) : phi_f_bar        in Z rest frame
//     (26) : final state combination.
//     (27) : e- beam gaussian spread
//     (28) : e+ beam gaussian spread
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "TTZSpring.h"

ClassImp(TTZSpring)
ClassImp(TTZSpringBuf)
ClassImp(TTZBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
TTZSpring::TTZSpring(const char *name, const char *title,
			 TTZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TTZSpringBuf("TTZSpringBuf", 
	 "TTZSpring event buffer", this);
  if( !bases ) { 
    TTZBases *bs=new TTZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TTZSpring::~TTZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t TTZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void TTZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TTZBases::TTZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TTZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TTZBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TTZBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TTZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TTZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TTZBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 23;
  	fNWILD = 7;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 24;
  	fNWILD = 8;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 28;
  	fNWILD = 11;
  	fIOFF  = 0;
  } else {
  	printf("TTZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }
#ifdef ZEROWIDTH
        fNWILD -= 6;
#endif

//
// Get ttbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"TTZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TTZBases.Roots","700."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("TTZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TTZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("TTZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TTZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassHiggs","120."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("TTZBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("TTZBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void TTZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TTZ generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t TTZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TTZBases::Userin()
{
//
//   Initialize User parameters for Bases
//
  JSFBases::Userin();  // Call JSFBases::Userin() for standard setup.

  // Copy class data member into common /usmprm/
  usmprm_.alfi = fAlphai;
  usmprm_.alfs = fAlphas;
  usmprm_.amsw = fMassW;
  usmprm_.amsz = fMassZ;
  usmprm_.amsh = fMassHiggs;
  usmprm_.amst = fMassTop;

  // Copy class data member into common /usrprm/
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit( 1, -1.0, 1.0, 50,"cos(theta_Z) ");
  Xhinit( 2,  0.0,360., 50,"phi_Z        ");
  Xhinit( 3,  0.0, 1.0, 50,"m(t-tb)/roots");
  Xhinit( 4, -1.0, 1.0, 50,"cos(theta_t) ");
  Xhinit( 5,  0.0,360., 50,"phi_t        ");
  Xhinit( 6, -1.0, 1.0, 50,"cos(theta_t)_lab    ");
  Xhinit( 7,  1.0,  5.,  4,"helicity combination");
  Dhinit( 9,0.,1.,50,-1.,1.,50,"E_Z/E_bm-cos(th_Z)");
  Dhinit(10,0.,1.,50, 0.,1.,50,"E_t/E_bm-E_tb/E_bm");
  return ;
}

//_____________________________________________________________________________
void TTZBases::Userout()
{
  usrout_();
}























