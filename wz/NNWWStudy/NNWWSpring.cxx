//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  NNWWSpring
//  
//  e+e- -> nu nubar W+ W-
//  
//  Integration variables.
//    Z( 1) : e- E_spread
//     ( 2) : e- beamstrahlung
//     ( 3) : e+ E_spread
//     ( 4) : e+ beamstrahlung
//     ( 5) : bremsstrahlng
//     ( 6) : helicity for initial states
//     ( 7) : m(W-)**2
//     ( 8) : m(W+)**2
//     ( 9) : cos(theta_W-)
//     (10) : phi_W-
//     (11) : cos(theta_fd)     in W-   rest frame
//     (12) : phi_fd            in W-   rest frame
//     (13) : cos(theta_fd_bar) in W+   rest frame
//     (14) : phi_fd_bar        in W+   rest frame
//     (15) : final state combination.
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "NNWWSpring.h"

ClassImp(NNWWSpring)
ClassImp(NNWWSpringBuf)
ClassImp(NNWWBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
NNWWSpring::NNWWSpring(const char *name, const char *title,
			 NNWWBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new NNWWSpringBuf("NNWWSpringBuf", 
	 "NNWWSpring event buffer", this);
  if( !bases ) { 
    NNWWBases *bs=new NNWWBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
NNWWSpring::~NNWWSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t NNWWSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void NNWWSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
NNWWBases::NNWWBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("NNWWBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("NNWWBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("NNWWBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("NNWWBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("NNWWBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("NNWWBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 16;
  	fNWILD = 11;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 17;
  	fNWILD = 12;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 21;
  	fNWILD = 14;
  	fIOFF  = 0;
  } else {
  	printf("NNWWBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get NNWW specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"NNWWBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("NNWWBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("NNWWBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("NNWWBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("NNWWBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("NNWWBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("NNWWBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("NNWWBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("NNWWBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("NNWWBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("NNWWBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("NNWWBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("NNWWBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("NNWWBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("NNWWBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("NNWWBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void NNWWBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->NNWW generator\n");
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
Double_t NNWWBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void NNWWBases::Userin()
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
  Xhinit( 7,   0., 1.0, 50,"WW/ROOTS      ");
  Xhinit( 8, -1.0, 1.0, 50,"cos_W-        ");
  Xhinit( 9,  0.0,360., 50,"phi_W-        ");
  Xhinit(10,  60.,110., 50,"m(W-)         ");
  Xhinit(11, -1.0, 1.0, 50,"cos(theta_fd) ");
  Xhinit(12,  0.0,360., 50,"phi_fd        ");
  Xhinit(13,  60.,110., 50,"m(W+)         ");
  Xhinit(14, -1.0, 1.0, 50,"cos(theta_fdb)");
  Xhinit(15,  0.0,360., 50,"phi_fdb       ");
  Xhinit(16,  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit(17,  1.0,25.0, 24,"W decay mode  ");
  Xhinit(18,  1.0,17.0, 16,"helicity      ");
  Xhinit(19, -20.,20.0, 50,"eta_1         ");
  Xhinit(20, -20.,20.0, 50,"eta_2         ");
  Xhinit(21,   0., 1.0, 50,"PT_W-/E_bm    ");
  Xhinit(22,   0., 1.0, 50,"PT_W+/E_bm    ");
  Xhinit(23,   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit(24,   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit(25,   0., 1.0, 50,"E_W-/E_bm     ");
  Xhinit(26,   0., 1.0, 50,"E_W+/E_bm     ");

  return ;
}

//_____________________________________________________________________________
void NNWWBases::Userout()
{
  usrout_();
}



















