#include <cmath>
#include <iostream>

using namespace std;
//_____________________________________________________________________________
// --------------------------
//  GamToSV
// --------------------------
double GamToSV(double M,  // parent mass
               double m1, // 1st daughter mass
               double m2, // 2nd daughter mass
               double a,  // coupling
               double cf) // color factor
{
   static const double kPi = acos(-1.);
   double x1   = pow(m1/M,2);
   double x2   = pow(m2/M,2);
   double beta = 1. - 2.*(x1+x2) + pow((x1-x2),2);

   if (beta <= 0.) return 0.;
   beta = sqrt(beta);

   double p1p2 = (M*M - m1*m1 - m2*m2)/2.;
   double tta = a*a*(-(4*m1*m1+m2*m2+4*p1p2)+pow(m2*m2+2*p1p2,2)/m2/m2);
   double fac = 1./(16.*kPi)/M;
   double gam = fac*tta*beta*cf;

   return gam;
}

int main()
{
   static const double kPi = acos(-1.);

   double s2w   = 0.230;
   double alpha = 1./128.;
   double e     = sqrt(4*kPi*alpha);
   double gw    = e/sqrt(s2w);
   double gz    = e/sqrt(s2w*(1.-s2w));

   double a  = gz/2;
   double M  = 180.;
   double m1 =  60.;
   double m2 =  91.18;
   double cf =   1.;

   double gam = GamToSV(M, m1, m2, a, cf);
   cerr << " gam = " << gam << endl;
}
