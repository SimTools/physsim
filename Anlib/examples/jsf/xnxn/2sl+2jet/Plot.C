//*----------------------------------------------------------------------------
//* ================
//*  Macro Plot.C
//* ================
//*
//* (Description)
//*	A sample macro to display ntuple contents.
//* (Usage)
//*	$ root -l
//*	......
//*	root [0] .L Plot.C
//*	root [1] Plot("jsf.root"); 
//*
//*	This will display distributions of some useful observables 
//*
//* (Update Record)
//*	2010/10/12	K.Fujii		Original version.
//*
//*----------------------------------------------------------------------------

const Double_t kPi  = TMath::Pi();
const Double_t kMx  = 200.0;
const Double_t kMst = 150.0;

void Plot(Char_t *filen = "jsf.root")
{
  //gROOT->Reset();
  Int_t nZoneX = 5;
  Int_t nZoneY = 3;

  TCanvas *c1 = new TCanvas("c1","",0,0, 300*nZoneX, 300*nZoneY);
  c1->SetHighLightColor(5);
  c1->SetFillColor(19);
  c1->Divide(nZoneX,nZoneY);

  gStyle->SetOptFit();

  cerr << filen << endl;
  TFile   *filep = new TFile(filen);
  TNtupleD *tup   = (TNtupleD *)gROOT->FindObject("hEvt");

  //--
  // Loop over preselected events
  //--
  Double_t kToDeg = 180./TMath::Pi();

  Int_t id = 0;
  id++; c1->cd(id); tup->Draw("pl:pt","");
  id++; c1->cd(id); tup->Draw("pt","");
  id++; c1->cd(id); tup->Draw("pl","");
  id++; c1->cd(id); tup->Draw("evis","");
  id++; c1->cd(id); tup->Draw("mx1:mx2","170.<mx1&&mx1<230.&&170.<mx2&&mx2<230.");
  id++; c1->cd(id); tup->Draw("mx1","170.<mx1&&mx1<230.&&170.<mx2&&mx2<230.");
  id++; c1->cd(id); tup->Draw("mx2","170.<mx1&&mx1<230.&&170.<mx2&&mx2<230.");
  id++; c1->cd(id); tup->Draw("csx1:csx2","-1<csx1&&csx1<1&&-1<csx2&&csx2<1");
  id++; c1->cd(id); tup->Draw("csx1","-1<csx1&&csx1<1&&-1<csx2&&csx2<1");
  id++; c1->cd(id); tup->Draw("csx2","-1<csx1&&csx1<1&&-1<csx2&&csx2<1");
  id++; c1->cd(id); tup->Draw("cstau1h:cstau2h","-1<cstau1h&&cstau1h<1&&-1<cstau2h&&cstau2h<1");
  id++; c1->cd(id); tup->Draw("cstau1h","-1<cstau1h&&cstau1h<1&&-1<cstau2h&&cstau2h<1");
  id++; c1->cd(id); tup->Draw("cstau2h","-1<cstau1h&&cstau1h<1&&-1<cstau2h&&cstau2h<1");
  id++; c1->cd(id); tup->Draw("fitau1h","");
  id++; c1->cd(id); tup->Draw("fitau2h","");
}
