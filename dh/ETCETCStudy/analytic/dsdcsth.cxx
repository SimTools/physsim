#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

double dsgdcsth(double csth, double rs, double M)
{
  const double kAlpha = 1./128.;
  const double kSin2W = 0.230;
  const double kMz    = 91.18;
  const double kGmz   = 2.5;
  const double kPi    = acos(-1.);
  const double kCos2W = 1. - kSin2W;
  const double kSinW  = sqrt(kSin2W);
  const double kCosW  = sqrt(kCos2W);
  const double kGeV2fb  = 0.389379292e12;   // GeV to fb

  double e   = sqrt(4.*kPi*kAlpha);
  double gw  = e/kSinW;
  double gz  = gw/kCosW;
  double gle = gz * (-0.5 + kSin2W);
  double gre = gz * (       kSin2W);

  double s = rs*rs;
  const complex<double> kI(0.,1.);
  double beta2 = 1. - 4*M*M/s;
  double beta  = sqrt(beta2);
  double lr  = pow(abs(e*e + gle*gle * s/(s - kMz*kMz + kI*kMz*kGmz)),2)*beta2*(1-csth)*(1+csth);
  double rl  = pow(abs(e*e + gle*gre * s/(s - kMz*kMz + kI*kMz*kGmz)),2)*beta2*(1-csth)*(1+csth);
  double fac = (1./2./s)*(1./4.)*(beta/8./kPi)*(1./4./kPi)*2*kPi;
  double dsg = fac*(lr + rl)*kGeV2fb;
  return dsg;
}

int main()
{
  double M  = 180.;
  double rs = 500.;

  double sgtot = 0.;
  int np = 100;
  double dcsth = 2./np;
  for (int icsth=0; icsth<=np; icsth++) {
    double csth = -1 + dcsth*icsth;
    double dsg  = dsgdcsth(csth, rs, M);
    cerr << csth << " " << dsg << endl;
    sgtot += dsg*dcsth;
  }
  cerr << " sgtot = " << sgtot << endl;
}
