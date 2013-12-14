TGraph *Plot(Char_t *infile = "xsection.dat")
{
   TCanvas *c1 = new TCanvas("c1", "c1",9,31,699,550);
   gStyle->SetOptStat(0);
   c1->ToggleEventStatus();
   c1->Range(347.5,-0.8764369,1074.664,5.646552);
   c1->SetBorderSize(2);
   c1->SetBottomMargin(0.1343612);

   TH2D *framep = new TH2D("h","",10,340.,352.,10,0.,5000.);
   TH2D &frame = *framep;
   frame.SetStats(0);
   frame.Draw();
   frame.GetXaxis()->SetTitle("#sqrt{s} [GeV]");
   frame.GetXaxis()->SetLabelSize(0.05);
   frame.GetXaxis()->SetTitleSize(0.06);
   frame.GetYaxis()->SetTitle("#sigma [fb]");
   frame.GetYaxis()->SetLabelSize(0.05);
   frame.GetYaxis()->SetTitleSize(0.06);
   frame.GetYaxis()->SetTitleOffset(0.77);
   frame.Draw("");

   double x[1000], y[1000], dy[1000];
   int npt = 0;
   ifstream in(infile);
   while ((in >> x[npt] >> y[npt] >> dy[npt])) npt++;
   TGraph *grp = new TGraph(npt, x, y);
   grp->Draw("same");
   return grp;
}
