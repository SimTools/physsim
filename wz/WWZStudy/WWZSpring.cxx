//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  WWZSpring
//  
//  e+e- -> W+W-Z
//  
//  Integration variables.
//    Z( 1) : e- E_spread
//     ( 2) : e- beamstrahlung
//     ( 3) : e+ E_spread
//     ( 4) : e+ beamstrahlung
//     ( 5) : bremsstrahlng
//     ( 6) : helicity for initial states
//     ( 7) : helicity combination for final states
//     ( 8) : m(W-)**2
//     ( 9) : m(W+)**2
//     (10) : m(Z)**2
//     (11) : m(W-W+)**2
//     (12) : cos(theta_Z)
//     (13) : phi_Z
//     (14) : cos(theta_W-)     in W-W+ rest frame
//     (15) : phi_W-            in W-W+ rest frame
//     (16) : cos(theta_fd)     in W-   rest frame
//     (17) : phi_fd            in W-   rest frame
//     (18) : cos(theta_fd_bar) in W+   rest frame
//     (19) : phi_fd_bar        in W+   rest frame
//     (20) : cos(theta_f_bar)  in Z    rest frame
//     (21) : phi_f_bar         in Z    rest frame
//     (22) : final state combination.
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "WWZSpring.h"

ClassImp(WWZSpring)
ClassImp(WWZSpringBuf)
ClassImp(WWZBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
WWZSpring::WWZSpring(const char *name, const char *title,
			 WWZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new WWZSpringBuf("WWZSpringBuf", 
	 "WWZSpring event buffer", this);
  if( !bases ) { 
    WWZBases *bs=new WWZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
WWZSpring::~WWZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t WWZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void WWZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
WWZBases::WWZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("WWZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("WWZBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("WWZBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("WWZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("WWZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("WWZBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 17;
  	fNWILD = 6;
  	fIOFF  = 5;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 18;
  	fNWILD = 7;
  	fIOFF  = 4;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 22;
  	fNWILD = 11;
  	fIOFF  = 0;
  } else {
  	printf("WWZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }
#ifdef ZEROWIDTH
	fNWILD -= 3;
#endif

//
// Get WWZ specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"WWZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("WWZBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("WWZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("WWZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("WWZBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("WWZBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("WWZBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("WWZBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("WWZBases.ZModesLo","1"),"%d",&fZModesLo);
  sscanf(gJSF->Env()->GetValue("WWZBases.ZModesHi","12"),"%d",&fZModesHi);
  sscanf(gJSF->Env()->GetValue("WWZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("WWZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("WWZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("WWZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("WWZBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("WWZBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("WWZBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("WWZBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void WWZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->WWZ generator\n");
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
  printf("  Z0 Decey Mode Lo =%d\n",fZModesLo);
  printf("                Hi =%d\n",fZModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t WWZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void WWZBases::Userin()
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
  usrprm_.imd1lo = fWmModesLo;
  usrprm_.imd1hi = fWmModesHi;
  usrprm_.imd2lo = fWpModesLo;
  usrprm_.imd2hi = fWpModesHi;
  usrprm_.imd3lo = fZModesLo;
  usrprm_.imd3hi = fZModesHi;

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
      Xhinit( 3,  0.0, 1.0, 50,"m(W-W+)/roots");
      Xhinit( 4, -1.0, 1.0, 50,"cos(theta_W-)");
      Xhinit( 5,  0.0,360., 50,"phi_W-       ");
      Xhinit( 7,  1.0,  5.,  4,"helicity combination  ");
      Dhinit( 9,0.,1.,50,-1.,1.,50,"E_Z/E_bm-cos(theta_Z)");
      Dhinit(10,0.,1.,50, 0.,1.,50,"E_W-/E_bm-E_W+/E_bm ");
  return ;
}

//_____________________________________________________________________________
void WWZBases::Userout()
{
  usrout_();
}
















