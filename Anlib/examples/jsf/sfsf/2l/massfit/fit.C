{
  gROOT->LoadMacro("massfit/lineshape.C");

  TF1 *func = new TF1("dNdEl",dNdEl,0.,50.,3);      
#if 1
// tanb = +3.0
#if 1
  func->SetParameters(143.169296,119.233253,200.);     

  hElFinal->Rebin(4);

  hElFinal->Fit("dNdEl","ri");

  TH2F *frame = new TH2F("Frame","Cont", 2, 142.5, 144.5, 2, 118.5, 120.5);
#else
  func->SetParameters(143.169296,119.233253,752.);     

  hElFinal->Rebin(5);

  hElFinal->Fit("dNdEl","ri");

  TH2F *frame = new TH2F("Frame","Cont", 2, 142.5, 144.2, 2, 118.5, 120.5);
#endif
  frame->Draw();
#else
// tanb = +10.0
#if 0
  func->SetParameters(144.3687,121.7502,95.94);     

  hElFinal->Rebin(4);
#else
  func->SetParameters(144.3687,121.7502,250.);     

  hElFinal->Rebin(4);
#endif

  hElFinal->Fit("dNdEl","ri");

  TH2F *frame = new TH2F("Frame","Cont", 2, 143.0, 146.5, 2, 120.0, 123.0);
  frame->Draw();
#endif

  gMinuit->SetErrorDef(4.61);
  TGraph *gr3 = (TGraph*)gMinuit->Contour(40,0,1); 
  gr3->SetFillColor(42);
  gr3->Draw("lf");

  gMinuit->SetErrorDef(2.28);
  TGraph *gr2 = (TGraph*)gMinuit->Contour(40,0,1); 
  gr2->SetFillColor(38);
  gr2->Draw("lf");

  gMinuit->SetErrorDef(1.);
  TGraph *gr1 = (TGraph*)gMinuit->Contour(40,0,1); 
  gr1->SetFillColor(3);
  gr1->Draw("lf");



}
