#include <iostream>

#include "TNtuple.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TH1D.h"

extern "C" {
  void isr_function_( double *factor, double *x, double *eps, int *LLA_order);
  void isr_remnant_(double *x, double *x0, double *xpt, double *xphi, double *sqrts, double pmom[4]);
};

int main()
{
  double ene=1000.0;
  double ame=0.00051;
  double alphai=137.0;
  double alpha=1.0/alphai;
  double pi=TMath::ACos(-1.0);
  int LLA_order=3;
  double eps=alpha/pi*2*TMath::Log(ene/ame);

  TFile *f=new TFile("funcplot.root","RECREATE");
  TNtuple *nt=new TNtuple("nt","strfun","x:fi:e:pt");
  TH1D *hisr=new TH1D("hisr","ISR energy",200,0.0,1.1);

  int nevt=10000;
  double x0=1.0e0;
  double pmom[4];
  for (int i=0;i<nevt;i++) {
    double x=gRandom->Rndm();
    double xpt=gRandom->Rndm();
    double xphi=gRandom->Rndm();
    double factor=1.0e0;

    isr_function_(&factor, &x, &eps, &LLA_order);
    isr_remnant_(&x, &x0, &xpt, &xphi, &ene, pmom);
//    std::cout << " factor=" << factor << " x=" << x << std::endl;
    double eisr=pmom[0];
    double ptisr=TMath::Sqrt(pmom[1]*pmom[1]+pmom[2]*pmom[2]);

    nt->Fill((float)x, (float)factor, eisr, ptisr);
    hisr->Fill((float)x, (float)factor);

  }
  f->Write();
  f->Close();
}

