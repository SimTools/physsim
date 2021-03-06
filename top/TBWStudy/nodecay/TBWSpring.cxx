//////////////////////////////////////////////////////////////////
//
//  TBWSpring
//  
//  e+e- -> tbwbar
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : e- helicity
//     ( 5) : final state helicity combination
//     ( 6) : m(t_bar)**2
//     ( 7) : cos(theta_t)
//     ( 8) : phi_t
//     ( 9) : cos(theta_q_bar) in t_bar rest frame
//     (10) : phi_q_bar        in t_bar rest frame
//     (11) : e- beam gaussian spread
//     (12) : e+ beam gaussian spread
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "TBWSpring.h"

ClassImp(TBWSpring)
ClassImp(TBWSpringBuf)
ClassImp(TBWBases)

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
TBWSpring::TBWSpring(const char *name, const char *title,
			 TBWBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TBWSpringBuf("TBWSpringBuf", 
	 "TBWSpring event buffer", this);
  if( !bases ) { 
    TBWBases *bs=new TBWBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TBWSpring::~TBWSpring()
{
  if( !fEventBuf ) { delete fEventBuf; fEventBuf = 0; }
}


//_____________________________________________________________________________
Bool_t TBWSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    TBWBases *bs = (TBWBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> TBWBases writbwen to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void TBWSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TBWBases::TBWBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TBWBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TBWBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TBWBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TBWBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TBWBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TBWBases.NCALL","40000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 7;
  	fNWILD = 4;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 8;
  	fNWILD = 5;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 12;
  	fNWILD = 7;
  	fIOFF  = 0;
  } else {
  	printf("TBWBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get tbwbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"TBWBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TBWBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("TBWBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TBWBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TBWBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TBWBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TBWBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TBWBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TBWBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void TBWBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TBW generator\n");
  printf("  Roots       = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t TBWBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  TBWBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t TBWBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TBWBases::Userin()
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
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.isrb   = fISRBM;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit("h01", -1.0,  1.0, 50,"cos(theta_t)  ");
  Xhinit("h02",  0.0, 360., 50,"phi_t         ");
  Xhinit("h03", -1.0,  1.0, 50,"cos(theta_b_bar) ");
  Xhinit("h04",  0.0, 360., 50,"phi_b_bar        ");
  Xhinit("h07",  0.0,   1., 50,"RS/ROOTS      ");
  Xhinit("h08",  0.0, 500., 50,"M(tbbar)      ");
  Xhinit("h09",  0.0, 500., 50,"M(bW+)        ");
  Xhinit("h10",  0.0, 500., 50,"M(bbarW-)     ");
  Xhinit("h11",  0.0, 250., 50,"|p|_cm        ");
  Xhinit("h12",  0.0, 250., 50,"|p|_lab       ");
  Xhinit("h13",  0.0, 250., 50,"|p|_lab       ");
  Xhinit("h14",  1.0,  13., 12,"helicity      ");
  Xhinit("h15", -1.0,  1.0, 50,"cos(theta_W-) ");
}

//_____________________________________________________________________________
void TBWBases::Userout()
{
  printf("End of TBWBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}
