//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  ZHHSpring
//  
//  e+e- -> ZHH
//  
//  Integration variables.
//    Z( 1) : e- beam
//     ( 2) : e- beam gaussian spread
//     ( 3) : e+ beam
//     ( 4) : e+ beam gaussian spread
//     ( 5) : bremsstrahlng
//     ( 6) : e- helicity
//     ( 7) : helicity combination for final states.
//     ( 8) : m(HH)**2
//     ( 9) : m(Z)**2
//     (10) : cos(theta_hh)
//     (11) : phi_hh
//     (12) : cos(theta_h)     in hh rest frame
//     (13) : phi_h            in hh rest frame
//     (14) : cos(theta_fb)    in Z  rest frame
//     (15) : phi_fb           in Z  rest frame
//     (16) : final state combination
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "ZHHSpring.h"

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

ClassImp(ZHHSpring)
ClassImp(ZHHSpringBuf)
ClassImp(ZHHBases)

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
};

//_____________________________________________________________________________
ZHHSpring::ZHHSpring(const char *name, const char *title,
			 ZHHBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new ZHHSpringBuf("ZHHSpringBuf", 
	 "ZHHSpring event buffer", this);
  if( !bases ) { 
    ZHHBases *bs=new ZHHBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
ZHHSpring::~ZHHSpring()
{
  //if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t ZHHSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    ZHHBases *bs = (ZHHBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> ZHHBases written to file" << endl;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void ZHHSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
ZHHBases::ZHHBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;				// set this for xhfill_

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("ZHHBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("ZHHBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("ZHHBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("ZHHBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("ZHHBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("ZHHBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 11;
  	fNWILD = 4;
  	fIOFF  = 5;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 12;
  	fNWILD = 5;
  	fIOFF  = 4;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 16;
  	fNWILD = 9;
  	fIOFF  = 0;
  } else {
  	printf("ZHHBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get ZHH specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"ZHHBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("ZHHBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("ZHHBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("ZHHBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("ZHHBases.Z0ModesLo","1"),"%d",&fZ0ModesLo);
  sscanf(gJSF->Env()->GetValue("ZHHBases.Z0ModesHi","12"),"%d",&fZ0ModesHi);
  sscanf(gJSF->Env()->GetValue("ZHHBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("ZHHBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("ZHHBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("ZHHBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("ZHHBases.MassHiggs","120."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("ZHHBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void ZHHBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->ZHH generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Z0 Decey Mode Lo =%d\n",fZ0ModesLo);
  printf("                Hi =%d\n",fZ0ModesHi);
  printf("  H0 Decey Mode = currently bbar only\n");
//  printf("  H0 Decey Mode Lo =%d\n",fH0ModesLo);
//  printf("                Hi =%d\n",fH0ModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);
}

//_____________________________________________________________________________
Double_t ZHHBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
Double_t ZHHBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  ZHHBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
void ZHHBases::Userin()
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
  usrprm_.imd1lo = fZ0ModesLo;
  usrprm_.imd1hi = fZ0ModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nZHH = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit("h01",  0.0,  1., 50,"rs/roots      ");
  Xhinit("h02",  60.,110., 50,"m(Z0)         ");
  Xhinit("h03", 200.,500., 50,"m(HH)         ");
  Xhinit("h04",  0.0,  1., 50,"miss/roots    ");
  Xhinit("h05", -1.0, 1.0, 50,"cos(theta_HH) ");
  Xhinit("h06",  0.0,360., 50,"phi_HH        ");
  Xhinit("h07", -1.0, 1.0, 50,"cos(theta_H)  ");
  Xhinit("h08",  0.0,360., 50,"phi_H         ");
  Xhinit("h09", -1.0, 1.0, 50,"cos(theta_fb) ");
  Xhinit("h10",  0.0,360., 50,"phi_fb        ");
  Xhinit("h11",  1.0,  9.,  8,"helicity      ");
  Xhinit("h12",  1.0, 13., 12,"Z decay mode  ");
}

//_____________________________________________________________________________
void ZHHBases::Userout()
{
  printf("End of Bases of ee --> sf sf process\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iteration  = %d\n",GetNoOfIterate());  
}
