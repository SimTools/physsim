Int_t formula(){
const Double_t   kGev2fb=0.389379292e12;
  Double_t cosT[201];
  Double_t xs[3][201];
  Double_t Ms=2;
  Double_t Mh=0.120;
  Double_t s=1;
  Double_t beta=TMath::Sqrt(1-((4*Mh*Mh)/s));
  Int_t iloop=0;
  Int_t mloop=0;

  for (mloop=0;mloop<=2;mloop++) {
    for (iloop=0;iloop<=200;iloop++) {
      cosT[iloop]=((iloop/(float)100)-1);
      if (mloop==0) {
	Ms=2;
      }
      if (mloop==1) {
	Ms=2.5;
      }
      if (mloop==2) {
	Ms=3;
      }
      xs[mloop][iloop]=(3.14159/((Double_t)128*(Ms*Ms*Ms*Ms*Ms*Ms*Ms*Ms)))*(s*s*s)*(beta*beta*beta*beta*beta)*((1-(cosT[iloop]*cosT[iloop]))*(cosT[iloop]*cosT[iloop]))*(kGev2fb/1000000);
      // cout << "iloop: " << iloop << " / " << mloop <<"  " << cosT[iloop] << " " << xs[mloop][iloop] << endl;
    }
  }

    //select postscript output type
  //  Int_t type = 111;   //portrait  ps
      Int_t type = 112;   //landscape ps
    // Int_t type = 113;   //eps

  //create a postscript file and set the paper size
  TPostScript ps("formula.ps",type);
  ps.Range(20,26);  //set x,y of printed page
  TCanvas c1("c1","canvas",800,600);

  gROOT->SetStyle("Plain");
  gPad->SetFillStyle(4000);
  gPad->SetFillColor(kWhite);
  /*
  gStyle->SetOptStat(0000000);
  gPad->SetFillStyle(4000);
  */
  TGraph  *distrib0= new TGraph(201,cosT,xs[0]);
  TGraph  *distrib1= new TGraph(201,cosT,xs[1]);
  TGraph  *distrib2= new TGraph(201,cosT,xs[2]);
  TMultiGraph *mg= new TMultiGraph();
  

  distrib0->SetLineColor(kRed);
  distrib0->SetLineStyle(0);
  distrib0->SetMarkerColor(kWhite);
  distrib1->SetLineColor(kGreen);
  distrib1->SetLineStyle(1);
  distrib2->SetLineColor(kBlue);
  distrib2->SetLineStyle(2);
  mg->SetMinimum(0);
  mg->Add(distrib0);
  
  mg->Add(distrib1);
  mg->Add(distrib2);
  
  mg->Draw("ALP");
  
  mg->GetXaxis()->SetTitle("cos #theta_{jj}");
  mg->GetYaxis()->SetTitle("#frac{d#sigma}{dcos #theta_{jj}} [fb]");
  mg->GetXaxis()->SetLimits(-1,1);

  distrib0->SetFillStyle(000);
  distrib1->SetFillStyle(000);
  distrib2->SetFillStyle(000);

  TText *t = new TText();
  t->SetTextFont(32);
  t->SetTextSize(0.06);
  t->SetTextAlign(33);
  t->SetTextColor(kRed);
  t->DrawText(0.45,6,"Ms=2 TeV");
  t->SetTextColor(kGreen);
  t->DrawText(0.9,2.1,"Ms=2.5 TeV");
  t->SetTextColor(kBlue);
  t->DrawText(0.85,.9,"Ms=3 TeV");


  /*
  TLegend *theLegend=new TLegend(.85,.75,.98,.98);
  theLegend->AddEntry(distrib0,"Ms=2");
  theLegend->AddEntry(distrib1,"Ms=2.5");
  theLegend->AddEntry(distrib2,"Ms=3");
  theLegend->Draw();
  */

  c1->Update();
  
  c1->Print("formula.gif");
  
  ps.Close();

  return 0;
} // formula
