{
   TCanvas *c1 = new TCanvas("c1", "c1",9,53,591,561);
   gStyle->SetOptStat(0);
   c1->ToggleEventStatus();
   c1->Range(-16.32353,-0.05370844,1672.941,0.3383632);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.1517368);
   c1->SetBottomMargin(0.1369863);

   TH2D frame("h","",10,240.,1500.,10,0.,1.0);
   frame.SetStats(0);
   frame.Draw();
   frame.GetXaxis()->SetTitle("#sqrt{s} [GeV]");
   frame.GetXaxis()->SetLabelSize(0.05);
   frame.GetXaxis()->SetTitleSize(0.06);
   frame.GetYaxis()->SetTitle("#sigma [fb]");
   frame.GetYaxis()->SetLabelSize(0.05);
   frame.GetYaxis()->SetTitleSize(0.06);
   frame.GetYaxis()->SetTitleOffset(1.15);
   frame.Draw("");

   double x[1000], y[1000], dy[1000];
   int npt = 0;
   ifstream in("xsection.zzz.dat");
   while ((in >> x[npt] >> y[npt] >> dy[npt])) npt++;
   TGraph gr(npt, x, y);
   gr.Draw("pc same");
#if 0
   npt = 0;
   in.close();
   in.open("xsection.zhh.dat");
   while ((in >> x[npt] >> y[npt] >> dy[npt])) npt++;
   TGraph gr2(npt, x, y);
   gr2.Draw("same");
#endif
}
