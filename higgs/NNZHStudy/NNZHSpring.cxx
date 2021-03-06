//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  NNZHSpring
//  
//  e+e- -> nu nu Z
//  
//  Integration variables.
//     Z( 1) : e- beam
//      ( 2) : e+ beam
//      ( 3) : bremsstrahlng
//      ( 4) : helicity
//      ( 5) : xi   (E_nu)
//      ( 6) : m(ZH)**2
//      ( 7) : m(Z)**2
//      ( 8) : cos_Z        in (Z H)  frame
//      ( 9) : phi_Z        in (Z H)  frame
//      (10) : eta  (cos_nu)
//      (11) : zeta (cos_nb) in (nb Z H) frame
//      (12) : phi_nu
//      (13) : phi_nb  	   in (nb Z H) frame
//      (14) : cos_fb       in Z      frame
//      (15) : phi_f        in Z      frame
//      (16) : final state combination.
//      (17) : E_beam spread of E-
//      (18) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "NNZHSpring.h"

ClassImp(NNZHSpring)
ClassImp(NNZHSpringBuf)
ClassImp(NNZHBases)

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
NNZHSpring::NNZHSpring(const char *name, const char *title,
			 NNZHBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new NNZHSpringBuf("NNZHSpringBuf", 
	 "NNZHSpring event buffer", this);
  if( !bases ) { 
    NNZHBases *bs=new NNZHBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
NNZHSpring::~NNZHSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t NNZHSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    NNZHBases *bs = (NNZHBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> NNZHBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void NNZHSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }


//_____________________________________________________________________________
NNZHBases::NNZHBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("NNZHBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("NNZHBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 13;
  	fNWILD = 5;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 14;
  	fNWILD = 6;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 18;
  	fNWILD = 8;
  	fIOFF  = 0;
  } else {
  	printf("NNZHBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }
#ifdef ZEROWIDTH
	fNWILD -= 3;
#endif

//
// Get NNZH specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"NNZHBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("NNZHBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("NNZHBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("NNZHBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ZModesLo","1"),"%d",&fZModesLo);
  sscanf(gJSF->Env()->GetValue("NNZHBases.ZModesHi","12"),"%d",&fZModesHi);
  sscanf(gJSF->Env()->GetValue("NNZHBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("NNZHBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("NNZHBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("NNZHBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("NNZHBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("NNZHBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void NNZHBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->NNZH generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Z0 Decey Mode Lo =%d\n",fZModesLo);
  printf("                Hi =%d\n",fZModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t NNZHBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  NNZHBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t NNZHBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  //  cerr << "func =  " << val << endl;
  return val;

}

//_____________________________________________________________________________
void NNZHBases::Userin()
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
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imdmin = fZModesLo;
  usrprm_.imdmax = fZModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms
  Xhinit("h01", -1.0, 1.0, 50,"cos(n)        ");
  Xhinit("h02",  0.0,360., 50,"phi_n         ");
  Xhinit("h03",   0., 1.0, 50,"E_n/E_bm      ");
  Xhinit("h04",  60.,110., 50,"m(Z)          ");
  Xhinit("h05", -1.0, 1.0, 50,"cos(nb)       ");
  Xhinit("h06",  0.0,360., 50,"phi_nb        ");
  Xhinit("h07", -1.0, 1.0, 50,"cos(theta_fb) ");
  Xhinit("h08",  0.0,360., 50,"phi_fb        ");
  Xhinit("h09",  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit("h10",  1.0,13.0, 12,"Z decay mode  ");
  Xhinit("h11",  1.0,17.0, 16,"helicity      ");
  Xhinit("h12", -10.,20.0, 40,"eta_n         ");
  Xhinit("h13",   0., 1.0, 50,"PT_Z/E_bm     ");
  Xhinit("h14",   0., 1.0, 50,"E_n/E_bm      ");
  Xhinit("h15",   0., 1.0, 50,"E_nb/E_bm     ");
  Xhinit("h16",   0., 1.0, 50,"E_Z/E_bm      ");
  Xhinit("h17",  200.,500.,50,"m(ZH)          ");
}

//_____________________________________________________________________________
void NNZHBases::Userout()
{
  printf("End of NNZHBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}






















