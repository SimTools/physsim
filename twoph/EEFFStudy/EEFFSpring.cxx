//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  EEFFSpring
//  
//  e+e- -> EEFFbar
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "EEFFSpring.h"

ClassImp(EEFFSpring)
ClassImp(EEFFSpringBuf)
ClassImp(EEFFBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
extern int snprintf ( char *, size_t, const char *, ... );
};

//_____________________________________________________________________________
EEFFSpring::EEFFSpring(const char *name, const char *title,
			 EEFFBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new EEFFSpringBuf("EEFFSpringBuf", 
	 "EEFFSpring event buffer", this);
  if( !bases ) { 
    EEFFBases *bs=new EEFFBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
EEFFSpring::~EEFFSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t EEFFSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void EEFFSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________

EEFFBases::EEFFBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("EEFFBases.ISRBM","1"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("EEFFBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("EEFFBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("EEFFBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("EEFFBases.ITMX2","20"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("EEFFBases.NCALL","100000"),"%d",&fNCALL);

  Int_t fIOFF = 0;

  fNDIM  = 7;
  fNWILD = 7;

//   if ( fISRBM == 1 ) {
//   	fNDIM  = 9;
//   	fNWILD = 3;
//   	fIOFF  = 3;
//   } else if ( fISRBM == 2 ) {
//    	fNDIM  = 10;
//   	fNWILD = 4;
//   	fIOFF  = 2;
//   } else if ( fISRBM == 3 ) {
//   	fNDIM  = 14;
//   	fNWILD = 6;
//   	fIOFF  = 0;
//   } else {
//   	printf("EEFFBases: Invalid ISRBM = %d\n",fISRBM);
//   	printf("    Will STOP immediately\n");
//   	exit(1);
//   }

//
// Get EEFFbar specific parameters.
//
  Char_t pname[40], pvalue[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"EEFFBases.X%i2.2Range",i+1);
    switch (i) {
      case 1: {
        sscanf(gJSF->Env()->GetValue(pname,"-1.0 1.0"),
		 "%lg%lg",&fXL[i],&fXU[i]);
        break;
        }
      case 2: {
        Double_t x3l = 0.0;
        Double_t x3u = 2*TMath::Pi();
        sprintf(pvalue,"%g %g",x3l,x3u);
        sscanf(gJSF->Env()->GetValue(pname,pvalue),
               "%lg%lg",&fXL[i],&fXU[i]);
        break;
        }
      default: {
        sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),
               "%lg%lg",&fXL[i],&fXU[i]);
        break;
        }
    }
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("EEFFBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("EEFFBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("EEFFBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("EEFFBases.Alphai","137."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("EEFFBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("EEFFBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("EEFFBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("EEFFBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("EEFFBases.MassTop","170."),"%lg",&fMassTop);  
  sscanf(gJSF->Env()->GetValue("EEFFBases.Generation","2"),"%d",&fGeneration);
  sscanf(gJSF->Env()->GetValue("EEFFBases.Isospin","2"),"%d",&fIsospin);
  sscanf(gJSF->Env()->GetValue("EEFFBases.LorQ","1"),"%d",&fLorQ);
  sscanf(gJSF->Env()->GetValue("EEFFBases.EeCut","0."),"%lg",&fEeCut);
  sscanf(gJSF->Env()->GetValue("EEFFBases.CoseMin","-1."),"%lg",&fCoseMin);
  sscanf(gJSF->Env()->GetValue("EEFFBases.CoseMax","1."),"%lg",&fCoseMax);
  sscanf(gJSF->Env()->GetValue("EEFFBases.EpCut","0."),"%lg",&fEpCut);
  sscanf(gJSF->Env()->GetValue("EEFFBases.CospMin","-1."),"%lg",&fCospMin);
  sscanf(gJSF->Env()->GetValue("EEFFBases.CospMax","1."),"%lg",&fCospMax);
  sscanf(gJSF->Env()->GetValue("EEFFBases.EfCut","0."),"%lg",&fEfCut);
  sscanf(gJSF->Env()->GetValue("EEFFBases.CosfCut","0.8"),"%lg",&fCosfCut);
  sscanf(gJSF->Env()->GetValue("EEFFBases.MassffMin","5."),"%lg",&fMassffMin);

  fPrintInfo = gJSF->Env()->GetValue("EEFFBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("EEFFBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void EEFFBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->EEFF generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Generation                = %d\n",fGeneration);
  printf("  Isospin((1,2)=(up,down))  = %d\n",fIsospin);
  printf("  Lepton/Quark((1,2)=(l,q)) = %d\n",fLorQ);
  printf("  Ee- cut      = %g\n",fEeCut);
  printf("  cos(e-)_min  = %g\n",fCoseMin);
  printf("  cos(e-)_max  = %g\n",fCoseMax);
  printf("  Ee+ cut      = %g\n",fEpCut);
  printf("  cos(e+)_min  = %g\n",fCospMin);
  printf("  cos(e+)_max  = %g\n",fCospMax);
  printf("  Ef cut       = %g\n",fEfCut);
  printf("  |cos(f)|_cut = %g\n",fCosfCut);
  printf("  m(ff)_min    = %g\n",fMassffMin);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t EEFFBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void EEFFBases::Userin()
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
  usrprm_.igfr   = fGeneration;
  usrprm_.itfr   = fIsospin;
  usrprm_.lqfr   = fLorQ;
  usrprm_.pect   = fEeCut;
  usrprm_.cecl   = fCoseMin;
  usrprm_.cecu   = fCoseMax;
  usrprm_.ppct   = fEpCut;
  usrprm_.cpcl   = fCospMin;
  usrprm_.cpcu   = fCospMax;
  usrprm_.pfct   = fEfCut;
  usrprm_.cfct   = fCosfCut;
  usrprm_.wmnf   = fMassffMin;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

      Double_t qmx  = fRoots/4.;
      Xhinit( 1, -1.0, 1.0, 50,"cos(theta_e-)      ");
      Xhinit( 2, -1.0, 1.0, 50,"cos(theta_e+)      ");
      Xhinit( 3,  0.0, 1.0, 50,"E_e-/E_bm          ");
      Xhinit( 4,  0.0, 1.0, 50,"E_e+/E_bm          ");
      Xhinit( 5,  0.0, 4.0, 50,"(M_ff/E_bm)**2     ");
      Xhinit( 6, -1.0, 1.0, 50,"cos(theta_f)       ");
      Xhinit( 7,  0.0, 1.0, 50,"E_f/E_bm           ");
      Xhinit( 8,  0.0, 2.0, 50,"M_ff/E_bm          ");
      Xhinit( 9,  0.0, qmx, 50,"M_ff (GeV)         ");
  return ;
}

//_____________________________________________________________________________
void EEFFBases::Userout()
{
  usrout_();
}



















