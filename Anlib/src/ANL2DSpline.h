#ifndef __ANL2DSPLINE__
#define __ANL2DSPLINE__
//*************************************************************************
//* ===================
//*  ANL2DSpline Class
//* ===================
//*
//* (Description)
//*    2-dimensional splining class
//* (Requires)
//*     class TSpline
//* (Provides)
//*     class ANL2DSpline
//* (Update Recored)
//*    2000/12/27  K.Ikematsu   Original version.
//*
//*************************************************************************
//

#include "TSpline.h"
//_____________________________________________________________________
//  -----------------
//  ANL2DSpline Class
//  -----------------
//

class ANL2DSpline : public TNamed {
public:
  ANL2DSpline(Int_t nxbins = 0, Double_t xmin = 0., Double_t xmax = 0.,
	      Int_t nybins = 0, Double_t ymin = 0., Double_t ymax = 0.,
	      Double_t **data = 0,
	      const char *name = "T2DSpline",
	      const char *title = "2-dimensional spline");
  virtual ~ANL2DSpline();

  void     Initialize(Int_t nxbins, Double_t xmin, Double_t xmax,
		      Int_t nybins, Double_t ymin, Double_t ymax,
		      Double_t **data);
  Double_t GetValueAt(Double_t x, Double_t y);

private:
  Int_t    fNxbins;     //! No. of x-bins
  Double_t fXmin;       //! x-bin min
  Double_t fXmax;       //! x-bin max
  Int_t    fNybins;     //! No. of y-bins
  Double_t fYmin;       //! y-bin min
  Double_t fYmax;       //! y-bin max

  Double_t **fData;     //! pointer to 2D data
  TSpline  **fYSpline;  //! pointer to y-splines

  ClassDef(ANL2DSpline,1) // ANL2DSpline class
};

#endif
