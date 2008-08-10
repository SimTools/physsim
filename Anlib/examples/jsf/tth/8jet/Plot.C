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
//*	2008/07/13	K.Fujii		Original version.
//*
//*----------------------------------------------------------------------------

TH2D *hMhMw      = new TH2D("hMhMw"    , "", 120,   60.,  150.,
                                             100,   20.,  120.);
TH2D *hMwMw      = new TH2D("hMwMw"    , "", 100,   20.,  120.,
                                             100,   20.,  120.);
TH2D *hMtMt      = new TH2D("hMtMt"    , "", 120,  100.,  220.,
                                             120,  100.,  220.);
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);

TH2D *hPtPl      = new TH2D("hPtPl"    , "", 100,    0.,  200.,
                                             200, -200., +200.);
TH2D *hMttMh     = new TH2D("hMttMh"   , "", 120,  200.,  440.,
                                             120,   60.,  150.);

void Plot(Char_t *filen = "jsf.root")
{
  //gROOT->Reset();
  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");

  Int_t nZoneX = 4;
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
  // Preselection
  //--
#if 0
  TString sel("evis>0.");
#else
  TString sel("pt<20.&&abs(pl)<30.&&ycut>0.0015");
#endif
  tup->Draw(">>elist",sel.Data(),"goff");

  TEventList *elist = static_cast<TEventList*>(gROOT->FindObject("elist"));
  Int_t nlist = elist->GetN();
  cerr << "Nevent = " << nlist << endl;

  //--
  // Loop over preselected events
  //--
  Double_t kToDeg = 180./TMath::Pi();
  for (Int_t i=0; i<nlist; i++) {
#if 0
        Int_t event = elist->GetEntry(i);
        cerr << "i = " << i << " event = " << event << endl;
#else
        Int_t event = i;
#endif
        Double_t ntracks;
        Double_t evis;
        Double_t pt;
        Double_t pl;
        Double_t ycut;
        Double_t chi2;
        Double_t mh;
        Double_t mw1;
        Double_t mw2;
        Double_t mt1;
        Double_t mt2;
        Double_t mtt;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("mh",&mh);
        tup->SetBranchAddress("mw1",&mw1);
        tup->SetBranchAddress("mw2",&mw2);
        tup->SetBranchAddress("mt1",&mt1);
        tup->SetBranchAddress("mt2",&mt2);
        tup->SetBranchAddress("mtt",&mtt);

        tup->GetEntry(event);

        hMhMw->Fill(mh,  mw1, 1.);
        hMwMw->Fill(mw1, mw2, 1.);
        hMtMt->Fill(mt1, mt2, 1.);
	hEvis->Fill(evis, 1.);
	hPtPl->Fill(pt, pl, 1.);
        hMttMh->Fill(mtt,  mh, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hMhMw->Draw();
  id++; c1->cd(id); hMhMw->ProjectionX()->Draw();
  id++; c1->cd(id); hMhMw->ProjectionY()->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hMwMw->Draw();
  id++; c1->cd(id); hMwMw->ProjectionX()->Draw();
  id++; c1->cd(id); hMwMw->ProjectionY()->Draw();
  id++; c1->cd(id); hMtMt->Draw();
  id++; c1->cd(id); hMtMt->ProjectionX()->Draw();
  id++; c1->cd(id); hMtMt->ProjectionY()->Draw();
  id++; c1->cd(id); hMttMh->ProjectionX()->Draw();
}
