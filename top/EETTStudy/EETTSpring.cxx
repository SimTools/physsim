//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  EETTSpring
//  
//  e+e- -> e+ettbar
//  
//  Integration variables.
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : helicity
//     ( 5) : xi   
//     ( 6) : eta_1
//     ( 7) : eta_2
//     ( 8) : m(ttbar)**2
//     ( 9) : m(t)**2
//     (10) : m(tbar)**2
//     (11) : m(W+)**2
//     (12) : m(W-)**2
//     (13) : phi_1
//     (14) : phi_2 - phi_1
//     (15) : cos_t        in tt frame
//     (16) : phi_t        in tt frame
//     (17) : cos_b        in t  frame
//     (18) : phi_b   	   in t  frame
//     (19) : cos_bb       in tb frame
//     (20) : phi_bb       in tb frame
//     (21) : cos(theta_f_bar) in W+ rest frame
//     (22) : phi_f_bar        in W+ rest frame
//     (23) : cos(theta_f)     in W- rest frame
//     (24) : phi_f            in W- rest frame
//     (25) : final state combination.
//     (26) : E_beam spread of E-
//     (27) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "EETTSpring.h"

ClassImp(EETTSpring)
ClassImp(EETTSpringBuf)
ClassImp(EETTBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
EETTSpring::EETTSpring(const char *name, const char *title,
			 EETTBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new EETTSpringBuf("EETTSpringBuf", 
	 "EETTSpring event buffer", this);
  if( !bases ) { 
    EETTBases *bs=new EETTBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
EETTSpring::~EETTSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t EETTSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void EETTSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
EETTBases::EETTBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("EETTBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("EETTBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("EETTBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("EETTBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("EETTBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("EETTBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 22;
  	fNWILD = 12;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 23;
  	fNWILD = 13;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 27;
  	fNWILD = 15;
  	fIOFF  = 0;
  } else {
  	printf("EETTBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get EETT specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"EETTBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("EETTBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("EETTBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("EETTBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("EETTBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("EETTBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("EETTBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("EETTBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("EETTBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("EETTBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("EETTBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("EETTBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("EETTBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("EETTBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("EETTBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("EETTBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void EETTBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->EETT generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  W- Decey Mode Lo =%d\n",fWmModesLo);
  printf("                Hi =%d\n",fWmModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t EETTBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void EETTBases::Userin()
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
  usrprm_.sqrts = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imd1lo = fWmModesLo;
  usrprm_.imd1hi = fWmModesHi;
  usrprm_.imd2lo = fWpModesLo;
  usrprm_.imd2hi = fWpModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms
  Xhinit( 1,   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit( 2, -1.0, 1.0, 50,"cos_1         ");
  Xhinit( 3,  0.0,360., 50,"phi_1         ");
  Xhinit( 4,   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit( 5, -1.0, 1.0, 50,"cos_2         ");
  Xhinit( 6,  0.0,360., 50,"phi_2         ");
  Xhinit( 7,   0., 1.0, 50,"TT/ROOTS      ");
  Xhinit( 8, -1.0, 1.0, 50,"cos_t         ");
  Xhinit( 9,  0.0,360., 50,"phi_t         ");
  Xhinit(10, 100.,200., 50,"m_t           ");
  Xhinit(11, -1.0, 1.0, 50,"cos(theta_b ) ");
  Xhinit(12,  0.0,360., 50,"phi_b         ");
  Xhinit(13, 100.,200., 50,"m_tb          ");
  Xhinit(14, -1.0, 1.0, 50,"cos(theta_bb) ");
  Xhinit(15,  0.0,360., 50,"phi_bb        ");
  Xhinit(16,  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit(17,  1.0,25.0, 24,"W decay mode  ");
  Xhinit(18,  1.0, 9.0,  8,"helicity      ");
  Xhinit(19, -20.,20.0, 50,"eta_1         ");
  Xhinit(20, -20.,20.0, 50,"eta_2         ");
  Xhinit(21,   0., 1.0, 50,"PT_t/E_bm     ");
  Xhinit(22,   0., 1.0, 50,"PT_tb/E_bm    ");
  Xhinit(23,   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit(24,   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit(25,   0., 1.0, 50,"E_t/E_bm      ");
  Xhinit(26,   0., 1.0, 50,"E_tb/E_bm     ");

  return ;
}

//_____________________________________________________________________________
void EETTBases::Userout()
{
  usrout_();
}




















