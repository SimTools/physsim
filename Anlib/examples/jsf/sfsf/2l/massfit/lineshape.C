void CompEndPoints(Double_t ep, Double_t mp, Double_t md, 
                   Double_t &emin, Double_t &emax)
{
   Double_t gamma = ep/mp;
   Double_t beta  = TMath::Sqrt((gamma-1)*(gamma+1))/gamma;
   Double_t x     = md/mp;
   Double_t pmax  = (mp/2)*(1-x)*(1+x);
   Double_t emx   = TMath::Sqrt(pmax*pmax + md*md);

   emin           = ep - gamma*(emx + beta*pmax);
   emax           = ep - gamma*(emx - beta*pmax);
}

#if 1
Double_t dNdEl(Double_t *elp, Double_t *x)
{
   const Double_t a     = 3.e-4;
   const Double_t sqrt2 = TMath::Sqrt(2.);
   const Double_t ep    = 175.;

   Double_t el    = *elp;
   Double_t mp    = x[0];
   Double_t md    = x[1];

   Double_t emin;
   Double_t emax;
   CompEndPoints(ep, mp, md, emin, emax);

   Double_t sgmin = sqrt2*a*emin*emin;
   Double_t sgmax = sqrt2*a*emax*emax;
#if 1
   // tanb = +3.
   //
#if 0
   // Rebin(4)
#if 1
   Double_t par[9] = { 1.42607e+04,
                      -9.23741e-02,  4.28024e-03,  1.77691e-05,
                      -7.63172e-06,  2.31496e-07, -2.24428e-09,
                      -2.02253e-01, -9.47806e-01 };
#else
   Double_t par[9] = { 5.07738e+03,
                       2.96392e-02, -1.48248e-03, -1.68521e-05,
                       1.09916e-06,  3.15101e-08, -9.85088e-10,
                      -8.01102e-01, -6.86451e-01 };
#endif
#else
#if 0
   // Rebin(5)
   Double_t par[9] = { 1.91468e+04, 
                      -7.39607e-02,  1.00549e-03,  2.13598e-04,
                      -1.27711e-05,  2.87112e-07, -2.38770e-09,
                      -2.31772e-01, -1.35356e+00 };
#else
#if 0
   // Rebin(2)
   Double_t par[9] = { 3.26613e+03,
                      -9.54246e-03, -8.14102e-04,  1.54395e-04,
                      -1.00100e-05,  2.85620e-07, -2.98910e-09,
                      -3.46236e-01, -8.13385e-01 };
#else
   // Rebin(1)
   Double_t par[9] = { 1.63705e+04, 
                      -1.08021e-01,  5.46464e-03, -3.61827e-06, 
                      -8.48439e-06,  2.67534e-07, -2.59704e-09,
                      -2.04807e-01, -7.48093e-01 };
#endif
#endif
#endif
#else
   // tanb = +10.
   Double_t par[9] = { 1., 
                       4.50391e-02, -6.99964e-03, 2.34414e-04,
                       1.44801e-06, -1.61092e-07, 1.60443e-09,
                       -4.02155e-01, -7.03154e-01 };
#endif
   par[0] = x[2];
#else
Double_t dNdEl(Double_t *elp, Double_t *par)
{
   const Double_t a     = 3.e-4;
   const Double_t sqrt2 = TMath::Sqrt(2.);

   Double_t el    = *elp;
#if 1
   // (m0,mu,m2,tanb) = (70,400,250,+3)
   // mslr = 143.169296
   // msz1 = 119.233253
   //
   Double_t emin  = 11.3934;
   Double_t emax  = 42.2305;
#else
   // (m0,mu,m2,tanb) = (70,400,250,+10)
   // mslr = 144.368652
   // msz1 = 121.750183
   // 
   Double_t emin  = 10.9877;
   Double_t emax  = 39.5519;
#endif
   Double_t sgmin = sqrt2*a*emin*emin;
   Double_t sgmax = sqrt2*a*emax*emax;
#endif

   Double_t val = 0.25*(1+TMath::Erf((el-emin)/sgmin))
                      *(1-TMath::Erf((el-emax)/sgmax))
                      *par[0]*(1 + el*(par[1] 
                                 + el*(par[2]
                                 + el*(par[3]
                                 + el*(par[4] 
                                 + el*(par[5] + el*par[6])))))
             + par[7]*TMath::Exp(-TMath::Power((el-emin-par[8]/2)/par[8],2)/2));
   return val;
}

