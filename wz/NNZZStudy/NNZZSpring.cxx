//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  NNZZSpring
//  
//  e+e- -> nu nubar Z Z
//  
//  Integration variables.
//     Z( 1) : e- beam
//      ( 2) : e+ beam
//      ( 3) : bremsstrahlng
//      ( 4) : helicity
//      ( 5) : xi   
//      ( 6) : eta_1
//      ( 7) : eta_2
//      ( 8) : m(ZZ)**2
//      ( 9) : m(Z1)**2
//      (10) : m(Z2)**2
//      (11) : phi_1
//      (12) : phi_2 - phi_1
//      (13) : cos_Z1       in ZZ frame
//      (14) : phi_Z1       in ZZ frame
//      (15) : cos_f        in Z1 frame
//      (16) : phi_f        in Z1 frame
//      (17) : cos_f        in Z2 frame
//      (18) : phi_f        in Z2 frame
//      (19) : final state combination.
//      (20) : E_beam spread of E-
//      (21) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "NNZZSpring.h"

ClassImp(NNZZSpring)
ClassImp(NNZZSpringBuf)
ClassImp(NNZZBases)

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
NNZZSpring::NNZZSpring(const char *name, const char *title,
			 NNZZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new NNZZSpringBuf("NNZZSpringBuf", 
	 "NNZZSpring event buffer", this);
  if( !bases ) { 
    NNZZBases *bs=new NNZZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
NNZZSpring::~NNZZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t NNZZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    NNZZBases *bs = (NNZZBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> NNZZBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void NNZZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
NNZZBases::NNZZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("NNZZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("NNZZBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("NNZZBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("NNZZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("NNZZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("NNZZBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 16;
  	fNWILD = 7;
  	fIOFF  = 3;
        fISHUFL[ 0] =  4;
        fISHUFL[ 1] =  8;
        fISHUFL[ 2] =  5;
        fISHUFL[ 3] =  6;
        fISHUFL[ 4] =  7;
        fISHUFL[ 5] = 13;
        fISHUFL[ 6] = 14;
        fISHUFL[ 7] = 11;
        fISHUFL[ 8] = 12;
        fISHUFL[ 9] =  9;
        fISHUFL[10] = 10;
        for (Int_t i=11; i<fNDIM; i++) {
           fISHUFL[i] = i + fIOFF + 1;
        }
        for (Int_t i=0; i<fNDIM; i++) {
          fIG[i] = 1;
        }
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 17;
  	fNWILD =  8;
  	fIOFF  = 2;
        fISHUFL[ 0] =  4;
        fISHUFL[ 1] =  8;
        fISHUFL[ 2] =  5;
        fISHUFL[ 3] =  6;
        fISHUFL[ 4] =  7;
        fISHUFL[ 5] =  3;
        fISHUFL[ 6] = 13;
        fISHUFL[ 7] = 14;
        fISHUFL[ 8] = 11;
        fISHUFL[ 9] = 12;
        fISHUFL[10] =  9;
        fISHUFL[11] = 10;
        for (Int_t i=12; i<fNDIM; i++) {
           fISHUFL[i] = i + fIOFF + 1;
        }
        for (Int_t i=0; i<fNDIM; i++) {
          fIG[i] = 1;
        }
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 21;
  	fNWILD = 10;
  	fIOFF  = 0;
        fISHUFL[ 0] =  4; // helicity
        fISHUFL[ 1] =  8; // m_V1V2^2
        fISHUFL[ 2] =  5; // xi
        fISHUFL[ 3] =  6; // eta_1
        fISHUFL[ 4] =  7; // eta_2
        fISHUFL[ 5] =  3; // ISR
        fISHUFL[ 6] =  1; // e- beamstrahlung
        fISHUFL[ 7] =  2; // e+ beamstrahlung
        fISHUFL[ 8] = 13; // cos_V1
        fISHUFL[ 9] = 14; // phi_V1
        fISHUFL[10] = 11; // phi_1
        fISHUFL[11] = 12; // phi_2 - phi_1
        fISHUFL[12] =  9; // m_V1^2
        fISHUFL[13] = 10; // m_V2-2
        for (Int_t i=14; i<fNDIM; i++) {
           fISHUFL[i] = i + fIOFF + 1;
        }
        for (Int_t i=0; i<fNDIM; i++) {
          fIG[i] = 1;
        }
  } else {
  	printf("NNZZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get NNZZ specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"NNZZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
  }

  sscanf(gJSF->Env()->GetValue("NNZZBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("NNZZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("NNZZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Z1ModesLo","1"),"%d",&fZ1ModesLo);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Z1ModesHi","12"),"%d",&fZ1ModesHi);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Z2ModesLo","1"),"%d",&fZ2ModesLo);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Z2ModesHi","12"),"%d",&fZ2ModesHi);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("NNZZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("NNZZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("NNZZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("NNZZBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("NNZZBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void NNZZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->NNZZ generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Z Decey Mode Lo =%d\n",fZ1ModesLo);
  printf("                Hi =%d\n",fZ1ModesHi);
  printf("  Z Decey Mode Lo =%d\n",fZ2ModesLo);
  printf("                Hi =%d\n",fZ2ModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t NNZZBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  NNZZBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t NNZZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void NNZZBases::Userin()
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
  usrprm_.imd1lo = fZ1ModesLo;
  usrprm_.imd1hi = fZ1ModesHi;
  usrprm_.imd2lo = fZ2ModesLo;
  usrprm_.imd2hi = fZ2ModesHi;

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
  Xhinit("h07",   0., 1.0, 50,"ZZ/ROOTS      ");
  Xhinit("h08", -1.0, 1.0, 50,"cos_Z1        ");
  Xhinit("h09",  0.0,360., 50,"phi_Z1        ");
  Xhinit("h10",  60.,110., 50,"m(Z1)         ");
  Xhinit("h11", -1.0, 1.0, 50,"cos(theta_fd) ");
  Xhinit("h12",  0.0,360., 50,"phi_fd        ");
  Xhinit("h13",  60.,110., 50,"m(Z2)         ");
  Xhinit("h14", -1.0, 1.0, 50,"cos(theta_fdb)");
  Xhinit("h15",  0.0,360., 50,"phi_fdb       ");
  Xhinit("h16",  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit("h17",  1.0,25.0, 24,"Z decay mode  ");
  Xhinit("h18",  1.0,17.0, 16,"helicity      ");
  Xhinit("h19", -20.,20.0, 50,"eta_1         ");
  Xhinit("h20", -20.,20.0, 50,"eta_2         ");
  Xhinit("h21",   0., 1.0, 50,"PT_Z1/E_bm    ");
  Xhinit("h22",   0., 1.0, 50,"PT_Z2/E_bm    ");
  Xhinit("h23",   0., 1.0, 50,"E_1/E_bm      ");
  Xhinit("h24",   0., 1.0, 50,"E_2/E_bm      ");
  Xhinit("h25",   0., 1.0, 50,"E_Z1/E_bm     ");
  Xhinit("h26",   0., 1.0, 50,"E_Z2/E_bm     ");

}

//_____________________________________________________________________________
void NNZZBases::Userout()
{
  printf("End of NNZZBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}



















