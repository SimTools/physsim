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
//*	2008/11/18	K.Fujii		Original version.
//*
//*----------------------------------------------------------------------------

const Double_t kPi = TMath::Pi();

TH2D *hMwMw      = new TH2D("hMwMw"    , "", 130,   10.,  140.,
                                             130,   10.,  140.);
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);

TH2D *hPtPl      = new TH2D("hPtPl"    , "", 100,    0.,  200.,
                                             200, -200., +200.);
TH1D *hElmx      = new TH1D("hElmx"    , "", 100,    0.,  100.);
TH1D *hCjmx      = new TH1D("hCjmx"    , "", 100,   -1.,   +1.);
TH2D *hCwCw      = new TH2D("hCwCw"    , "", 100,   -1.,   +1.,
                                             100,   -1.,   +1.);
TH1D *hCosW      = new TH1D("hCosW"    , "", 100,   -1.,   +1.);
TH1D *hCsjh      = new TH1D("hCsjh"    , "", 100,   -1.,   +1.);
TH1D *hFijh      = new TH1D("hFijh"    , "", 100,  -kPi,  +kPi);
TH1D *hEw        = new TH1D("hEw"      , "", 125,    0.,  250.);

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
  TString sel("evis>0.");
  tup->Draw(">>elist",sel.Data(),"goff");

  TEventList *elist = static_cast<TEventList*>(gROOT->FindObject("elist"));
  Int_t nlist = elist->GetN();
  cerr << "Nevent = " << nlist << endl;

  //--
  // Loop over preselected events
  //--
  Double_t kToDeg = 180./TMath::Pi();
  for (Int_t i=0; i<nlist; i++) {
#if 1
        Int_t event = elist->GetEntry(i);
        cerr << "i = " << i << " event = " << event << endl;
#else
        Int_t event = i;
#endif
        Double_t ntracks;
        Double_t evis;
        Double_t pt;
        Double_t pl;
        Double_t elmax;
        Double_t ycut;
        Double_t chi2;
        Double_t nsols;
        Double_t mw1;
        Double_t mw2;
        Double_t csjmax;
        Double_t csw1;
        Double_t csw2;
        Double_t ew1;
        Double_t ew2;
        Double_t mm;
        Double_t acop;
	Double_t csj11h;
	Double_t fij11h;
	Double_t csj12h;
	Double_t fij12h;
	Double_t csj21h;
	Double_t fij21h;
	Double_t csj22h;
	Double_t fij22h;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("elmax",&elmax);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("nsols",&nsols);
        tup->SetBranchAddress("mw1",&mw1);
        tup->SetBranchAddress("mw2",&mw2);
        tup->SetBranchAddress("ew1",&ew1);
        tup->SetBranchAddress("ew2",&ew2);
        tup->SetBranchAddress("csjmax",&csjmax);
        tup->SetBranchAddress("csw1",&csw1);
        tup->SetBranchAddress("csw2",&csw2);
        tup->SetBranchAddress("mm",&mm);
        tup->SetBranchAddress("acop",&acop);
        tup->SetBranchAddress("csj11h",&csj11h);
        tup->SetBranchAddress("fij11h",&fij11h);
        tup->SetBranchAddress("csj12h",&csj12h);
        tup->SetBranchAddress("fij12h",&fij12h);
        tup->SetBranchAddress("csj21h",&csj21h);
        tup->SetBranchAddress("fij21h",&fij21h);
        tup->SetBranchAddress("csj22h",&csj22h);
        tup->SetBranchAddress("fij22h",&fij22h);

        tup->GetEntry(event);

        hMwMw->Fill(mw1, mw2, 1.);
	hEvis->Fill(evis, 1.);
	hElmx->Fill(elmax, 1.);
	hCjmx->Fill(csjmax, 1.);
	hCwCw->Fill(csw1, csw2, 1.);
	hCosW->Fill(csw1, 1.);
	hCosW->Fill(csw2, 1.);
	hEw  ->Fill(ew1 , 1.);
	hEw  ->Fill(ew2 , 1.);
	hPtPl->Fill(pt, pl, 1.);
	hCsjh->Fill(csj11h, 1.);
	hCsjh->Fill(csj12h, 1.);
	hCsjh->Fill(csj21h, 1.);
	hCsjh->Fill(csj22h, 1.);
	hFijh->Fill(fij11h, 1.);
	hFijh->Fill(fij12h, 1.);
	hFijh->Fill(fij21h, 1.);
	hFijh->Fill(fij22h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hPtPl->Draw();
  id++; c1->cd(id); hElmx->Draw();
  id++; c1->cd(id); hMwMw->Draw();
  id++; c1->cd(id); hMwMw->ProjectionX()->Draw();
  id++; c1->cd(id); hMwMw->ProjectionY()->Draw();
  id++; c1->cd(id); hEw ->Draw();
  id++; c1->cd(id); hCwCw->Draw();
  hCosW->SetMinimum(0.);
  id++; c1->cd(id); hCosW->Draw();
  hCsjh->SetMinimum(0.);
  id++; c1->cd(id); hCsjh->Draw();
  hFijh->SetMinimum(0.);
  id++; c1->cd(id); hFijh->Draw();
}
