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
#include "ANL2DSpline.h"
//_____________________________________________________________________
//  -----------------
//  ANL2DSpline Class
//  -----------------
//

ClassImp(ANL2DSpline)

//*--
//*  Constructors
//*--
ANL2DSpline::ANL2DSpline(Int_t nxbins, Double_t xmin, Double_t xmax,
			 Int_t nybins, Double_t ymin, Double_t ymax,
			 Double_t **data,
			 const char *name, const char *title)
  : TNamed(name,title),
    fNxbins(nxbins), fXmin(xmin), fXmax(xmax),
    fNybins(nybins), fYmin(ymin), fYmax(ymax),
    fData(data)
{
  Initialize(nxbins, xmin, xmax,
	     nybins, ymin, ymax, data);
}

//*--
//*  Destructor
//*--
ANL2DSpline::~ANL2DSpline() {

}

//*--
//*  Setters
//*--
void ANL2DSpline::Initialize(Int_t nxbins, Double_t xmin, Double_t xmax,
			     Int_t nybins, Double_t ymin, Double_t ymax,
			     Double_t **data) {
  if (!fData) return;
  fYSpline = new TSpline *[fNybins];
  for (Int_t i = 0; i < fNybins; i++) {
    Char_t spname[256];
    sprintf(spname,"yspline%5.5d",i);
    fYSpline[i] = new TSpline3(spname, fXmin, fXmax, fData[i], fNxbins);
  }
}

//*--
//*  Gettters
//*--
Double_t ANL2DSpline::GetValueAt(Double_t x, Double_t y) {
  Double_t data[fNybins];
  for (Int_t i = 0; i < fNybins; i++) {
    data[i] = fYSpline[i]->Eval(x);
  }
  TSpline3 xspline("Temp", fYmin, fYmax, data, fNybins);

  return xspline.Eval(y);
}
