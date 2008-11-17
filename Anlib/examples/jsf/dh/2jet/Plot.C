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

TH1D *hMz        = new TH1D("hMz"      , "",  80,   70.,  110.);
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);
TH2D *hPtPl      = new TH2D("hPtPl"    , "", 100,    0.,  200.,
                                             200, -200., +200.);
TH1D *hElmx      = new TH1D("hElmx"    , "", 100,    0.,  100.);
TH1D *hCjmx      = new TH1D("hCjmx"    , "", 100,   -1.,   +1.);
TH1D *hCosZ      = new TH1D("hCosZ"    , "", 100,   -1.,   +1.);
TH1D *hCsjh      = new TH1D("hCsjh"    , "", 100,   -1.,   +1.);
TH1D *hFijh      = new TH1D("hFijh"    , "", 100,  -kPi,  +kPi);
TH1D *hEz        = new TH1D("hEz"      , "", 150,    0.,  300.);
TH1D *hMM        = new TH1D("hMM"      , "", 250,    0.,  500.);

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
        Double_t mz;
        Double_t csjmax;
        Double_t csz;
        Double_t ez;
        Double_t mm;
	Double_t csj1h;
	Double_t fij1h;
	Double_t csj2h;
	Double_t fij2h;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("elmax",&elmax);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("mz",&mz);
        tup->SetBranchAddress("ez",&ez);
        tup->SetBranchAddress("csjmax",&csjmax);
        tup->SetBranchAddress("csz",&csz);
        tup->SetBranchAddress("mm",&mm);
        tup->SetBranchAddress("csj1h",&csj1h);
        tup->SetBranchAddress("fij1h",&fij1h);
        tup->SetBranchAddress("csj2h",&csj2h);
        tup->SetBranchAddress("fij2h",&fij2h);

        tup->GetEntry(event);

	hEvis->Fill(evis, 1.);
	hElmx->Fill(elmax, 1.);
	hCjmx->Fill(csjmax, 1.);
	hCosZ->Fill(csz, 1.);
	hEz  ->Fill(ez , 1.);
	hMz  ->Fill(mz , 1.);
	hMM  ->Fill(mm , 1.);
	hPtPl->Fill(pt, pl, 1.);
	hCsjh->Fill(csj1h, 1.);
	hCsjh->Fill(csj2h, 1.);
	hFijh->Fill(fij1h, 1.);
	hFijh->Fill(fij2h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hPtPl->Draw();
  id++; c1->cd(id); hElmx->Draw();
  id++; c1->cd(id); hMz  ->Draw();
  id++; c1->cd(id); hEz  ->Draw();
  id++; c1->cd(id); hMM  ->Draw();
  hCosZ->SetMinimum(0.);
  id++; c1->cd(id); hCosZ->Draw();
  hCsjh->SetMinimum(0.);
  id++; c1->cd(id); hCsjh->Draw();
  hFijh->SetMinimum(0.);
  id++; c1->cd(id); hFijh->Draw();
}
