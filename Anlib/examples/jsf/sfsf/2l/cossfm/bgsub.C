{
   gROOT->Reset();
#if 0
//
// BG for tanb = +10
//
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  p0           2.84561e+02   5.33443e+00   6.58365e-03   3.22008e-07
//   2  p1           3.81648e+00   9.17746e+00   1.13266e-02   6.25138e-07
//
   Double_t p0 =  2.84561e+02;
   Double_t p1 =  3.81648e+00;
#else
//
// BG for tanb = +10
//
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//  1  p0           3.03380e+02   5.50799e+00   1.17854e-02   7.86874e-08
//  2  p1          -2.58374e+00   9.40403e+00   2.01218e-02  -4.01674e-07
//
   Double_t p0 =  3.03380e+02;
   Double_t p1 =  -2.58374e+00;
#endif

   TCanvas c1("c1","cos(theta_sf)",10,10,710,710);
   c1.Range(0,0,1,1);
   c1.SetBorderSize(2);
   c1.SetFrameFillColor(0);

   c1.cd();
   TPad pad1("pad1","pad1",0.01, 0.50, 0.99,0.99);
   pad1.Range(-1.2, 0.,1.2,1400.);
   pad1.SetBorderMode(0);
   pad1.SetBorderSize(0);
   pad1.SetBottomMargin(0.);
   pad1.SetTopMargin(0.25);
   pad1.Draw();
   pad1.cd();

   TH1F h(*hCosSFmFinal);
   h.Rebin(5);
   TH1F hg(h);

   h.SetTitle("");
   h.SetMarkerColor(50);
   h.SetMarkerStyle(8);
   h.SetMarkerSize(1.5);
   h.SetStats(0);
   h.SetMaximum(890.);

   TH1F ho(h);
   gStyle->SetErrorX(0.);
   ho.SetStats(0);
   ho.SetTitle("");
   ho.SetMaximum(890.);
   ho.GetXaxis()->SetLabelSize(0.);
   ho.GetYaxis()->SetLabelSize(0.06);
   ho.Draw("e1");

   hCosSFm1G->ProjectionY();
   hCosSFm1G_py->Rebin(5);
   hCosSFm1G_py->Draw("same");

   TLatex *   tex = new TLatex(-0.936255,791.984,"a)");
   tex->SetTextSize(0.0735429);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.394329,770.166,"#sqrt{s} = 350GeV");
   tex->SetTextSize(0.0625115);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();


   c1.cd();
   TPad pad2("pad2","pad2",0.01, 0.01, 0.99, 0.50);
   pad2.Range(-1.2,-400.,1.2,1000.);
   pad2.SetBorderMode(0);
   pad2.SetBorderSize(0);
   pad2.SetBottomMargin(0.25);
   pad2.SetTopMargin(0.);
   pad2.Draw();
   pad2.cd();

   Double_t x, y, dy;
   Double_t sum = 0;
   for (Int_t i=1; i<h.GetNbinsX()+1; i++) {
      x  = h.GetBinCenter(i);
      y  = h.GetBinContent(i) - p0 - p1*x;
      dy = h.GetBinError(i);
      sum += y;
      cerr << i << " " << y << " " << dy << endl;
      h.SetBinContent(i,y);
      h.SetBinError(i,dy);
   }
   h.GetYaxis()->SetLabelSize(0.06);
   h.GetXaxis()->SetLabelSize(0.06);
   h.GetXaxis()->SetTitle("cos#theta");
   h.GetXaxis()->CenterTitle(kTRUE);
   h.GetXaxis()->SetTitleSize(0.08);
   h.Draw("e1");

   for (Int_t i=1; i<hg.GetNbinsX()+1; i++) {
      Double_t xl = hg.GetBinLowEdge(i);
      Double_t xh = xl + hg.GetBinWidth(i);
      y = (3*sum/4)*(xh-xl - (1/3.)*(xh*xh*xh-xl*xl*xl));
      hg.SetBinContent(i,y);
   }
   hg.Draw("h same");
   

         tex = new TLatex(-0.923038,770.002,"b)");
   tex->SetTextSize(0.0735429);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.515419,771.626,"100fb^{-1}");
   tex->SetTextSize(0.0625115);
   tex->SetLineWidth(2);
:   tex->Draw();
      tex = new TLatex(0.484582,644.483,"Pol.e^{-} = +0.90");
   tex->SetTextSize(0.0625115);
   tex->SetLineWidth(2);
   tex->Draw();
   pad2->Modified();

}

