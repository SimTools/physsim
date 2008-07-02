//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  EEWWSpring
//  
//  e+e- -> e+e-W+W-
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
#include "EEWWSpring.h"

ClassImp(EEWWSpring)
ClassImp(EEWWSpringBuf)
ClassImp(EEWWBases)

// Variable conversion.

#define fNDIM  ndim 
#define fNWILD nwild
#define fIG    ig
#define fXL    xl
#define fXU    xu
#define fNCALL ncall
#define fACC1  acc1
#define fACC2  acc2
#define fITMX1 itmx1
#define fITMX2 itmx2

#define Xhinit(id,xlo,xhi,n,title) H1Init(id,title,n,xlo,xhi)
#define Dhinit(id,xlo,xhi,nx,ylo,yhi,ny,title) H2Init(id,title,nx,xlo,xhi,ny,ylo,yhi)

extern "C" {
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
JSFBases *bso;			// need to map xhfill to h1fill
void xhfill_(char *t, double *x, double *w, int len)
{
	char tmp[1024];
	int i;
	for (i=0; i<len; i++) tmp[i] = t[i];
	tmp[len] = '\0';
	bso->H1Fill(tmp,*x,*w);
}
void dhfill_(char *t, double *x, double *y, double *w, int len)
{
	char tmp[1024];
	int i;
	for (i=0; i<len; i++) tmp[i] = t[i];
	tmp[len] = '\0';
	bso->H2Fill(tmp,*x,*y,*w);
}
};

//_____________________________________________________________________________
EEWWSpring::EEWWSpring(const char *name, const char *title,
			 EEWWBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new EEWWSpringBuf("EEWWSpringBuf", 
	 "EEWWSpring event buffer", this);
  if( !bases ) { 
    EEWWBases *bs=new EEWWBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
EEWWSpring::~EEWWSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t EEWWSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    EEWWBases *bs = (EEWWBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> EEWWBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void EEWWSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
EEWWBases::EEWWBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("EEWWBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("EEWWBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("EEWWBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("EEWWBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("EEWWBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("EEWWBases.NCALL","80000"),"%d",&fNCALL);

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
  	printf("EEWWBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get EEWW specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"EEWWBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("EEWWBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("EEWWBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("EEWWBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("EEWWBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("EEWWBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("EEWWBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("EEWWBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("EEWWBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("EEWWBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("EEWWBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("EEWWBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("EEWWBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("EEWWBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void EEWWBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->EEWW generator\n");
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


}

//_____________________________________________________________________________
Double_t EEWWBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  EEWWBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t EEWWBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void EEWWBases::Userin()
{
//
//   Initialize User parameters for Bases
//

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
  Xhinit("h01",   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit("h02", -1.0, 1.0, 50,"cos_1         ");
  Xhinit("h03",  0.0,360., 50,"phi_1         ");
  Xhinit("h04",   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit("h05", -1.0, 1.0, 50,"cos_2         ");
  Xhinit("h06",  0.0,360., 50,"phi_2         ");
  Xhinit("h07",   0., 1.0, 50,"WW/ROOTS      ");
  Xhinit("h08", -1.0, 1.0, 50,"cos_W-        ");
  Xhinit("h09",  0.0,360., 50,"phi_W-        ");
  Xhinit("h10",  60.,110., 50,"m(W-)         ");
  Xhinit("h11", -1.0, 1.0, 50,"cos(theta_fd) ");
  Xhinit("h12",  0.0,360., 50,"phi_fd        ");
  Xhinit("h13",  60.,110., 50,"m(W+)         ");
  Xhinit("h14", -1.0, 1.0, 50,"cos(theta_fdb)");
  Xhinit("h15",  0.0,360., 50,"phi_fdb       ");
  Xhinit("h16",  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit("h17",  1.0,25.0, 24,"W decay mode  ");
  Xhinit("h18",  1.0,17.0, 16,"helicity      ");
  Xhinit("h19", -20.,20.0, 50,"eta_1         ");
  Xhinit("h20", -20.,20.0, 50,"eta_2         ");
  Xhinit("h21",   0., 1.0, 50,"PT_W-/E_bm    ");
  Xhinit("h22",   0., 1.0, 50,"PT_W+/E_bm    ");
  Xhinit("h23",   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit("h24",   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit("h25",   0., 1.0, 50,"E_W-/E_bm     ");
  Xhinit("h26",   0., 1.0, 50,"E_W+/E_bm     ");

}

//_____________________________________________________________________________
void EEWWBases::Userout()
{
  printf("End of EEWWBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}



















