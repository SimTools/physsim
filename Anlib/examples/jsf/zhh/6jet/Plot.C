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
const Double_t kMh = 120.0;

TH2D *hMhMh      = new TH2D("hMhMh"    , "", 130,   30.,  160.,
                                             130,   30.,  160.);
TH1D *hMz        = new TH1D("hMz"      , "", 160,   40.,  120.);
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);

TH2D *hPtPl      = new TH2D("hPtPl"    , "", 100,    0.,  200.,
                                             200, -200., +200.);
TH1D *hElmx      = new TH1D("hElmx"    , "", 100,    0.,  100.);
TH1D *hCjmx      = new TH1D("hCjmx"    , "", 100,   -1.,   +1.);
TH2D *hChCh      = new TH2D("hChCh"    , "", 100,   -1.,   +1.,
                                             100,   -1.,   +1.);
TH1D *hCosH      = new TH1D("hCosH"    , "", 100,   -1.,   +1.);
TH1D *hCsjh      = new TH1D("hCsjh"    , "", 100,   -1.,   +1.);
TH1D *hFijh      = new TH1D("hFijh"    , "", 100,  -kPi,  +kPi);
TH1D *hEh        = new TH1D("hEh"      , "", 125,    0.,  250.);
TH1D *hCsjz      = new TH1D("hCsjz"    , "", 100,   -1.,   +1.);
TH1D *hFijz      = new TH1D("hFijz"    , "", 100,  -kPi,  +kPi);
TH1D *hPh        = new TH1D("hPh"      , "", 125,    0.,  250.);
TH1D *hMM        = new TH1D("hMM"      , "", 250,    0.,  500.);

void Plot(Char_t *filen = "jsf.root")
{
  //gROOT->Reset();
  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");

  Int_t nZoneX = 6;
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
        Double_t mh1;
        Double_t mh2;
        Double_t mz;
        Double_t csjmax;
        Double_t csh1;
        Double_t csh2;
        Double_t eh1;
        Double_t eh2;
        Double_t mm;
	Double_t csj11h;
	Double_t fij11h;
	Double_t csj12h;
	Double_t fij12h;
	Double_t csj21h;
	Double_t fij21h;
	Double_t csj22h;
	Double_t fij22h;
	Double_t csj31h;
	Double_t fij31h;
	Double_t csj32h;
	Double_t fij32h;
	Double_t pj1e;
	Double_t pj1x;
	Double_t pj1y;
	Double_t pj1z;
	Double_t pj2e;
	Double_t pj2x;
	Double_t pj2y;
	Double_t pj2z;
	Double_t pj3e;
	Double_t pj3x;
	Double_t pj3y;
	Double_t pj3z;
	Double_t pj4e;
	Double_t pj4x;
	Double_t pj4y;
	Double_t pj4z;
	Double_t pj5e;
	Double_t pj5x;
	Double_t pj5y;
	Double_t pj5z;
	Double_t pj6e;
	Double_t pj6x;
	Double_t pj6y;
	Double_t pj6z;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("elmax",&elmax);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("nsols",&nsols);
        tup->SetBranchAddress("mh1",&mh1);
        tup->SetBranchAddress("mh2",&mh2);
        tup->SetBranchAddress("mz",&mz);
        tup->SetBranchAddress("eh1",&eh1);
        tup->SetBranchAddress("eh2",&eh2);
        tup->SetBranchAddress("csjmax",&csjmax);
        tup->SetBranchAddress("csh1",&csh1);
        tup->SetBranchAddress("csh2",&csh2);
        tup->SetBranchAddress("mm",&mm);
        tup->SetBranchAddress("csj11h",&csj11h);
        tup->SetBranchAddress("fij11h",&fij11h);
        tup->SetBranchAddress("csj12h",&csj12h);
        tup->SetBranchAddress("fij12h",&fij12h);
        tup->SetBranchAddress("csj21h",&csj21h);
        tup->SetBranchAddress("fij21h",&fij21h);
        tup->SetBranchAddress("csj22h",&csj22h);
        tup->SetBranchAddress("fij22h",&fij22h);
        tup->SetBranchAddress("csj31h",&csj31h);
        tup->SetBranchAddress("fij31h",&fij31h);
        tup->SetBranchAddress("csj32h",&csj32h);
        tup->SetBranchAddress("fij32h",&fij32h);
        tup->SetBranchAddress("pj1e",&pj1e);
        tup->SetBranchAddress("pj1x",&pj1x);
        tup->SetBranchAddress("pj1y",&pj1y);
        tup->SetBranchAddress("pj1z",&pj1z);
        tup->SetBranchAddress("pj2e",&pj2e);
        tup->SetBranchAddress("pj2x",&pj2x);
        tup->SetBranchAddress("pj2y",&pj2y);
        tup->SetBranchAddress("pj2z",&pj2z);
        tup->SetBranchAddress("pj3e",&pj3e);
        tup->SetBranchAddress("pj3x",&pj3x);
        tup->SetBranchAddress("pj3y",&pj3y);
        tup->SetBranchAddress("pj3z",&pj3z);
        tup->SetBranchAddress("pj4e",&pj4e);
        tup->SetBranchAddress("pj4x",&pj4x);
        tup->SetBranchAddress("pj4y",&pj4y);
        tup->SetBranchAddress("pj4z",&pj4z);
        tup->SetBranchAddress("pj5e",&pj5e);
        tup->SetBranchAddress("pj5x",&pj5x);
        tup->SetBranchAddress("pj5y",&pj5y);
        tup->SetBranchAddress("pj5z",&pj5z);
        tup->SetBranchAddress("pj6e",&pj6e);
        tup->SetBranchAddress("pj6x",&pj6x);
        tup->SetBranchAddress("pj6y",&pj6y);
        tup->SetBranchAddress("pj6z",&pj6z);

        tup->GetEntry(event);

	ANL4DVector qj11(pj1e,pj1x,pj1y,pj1z);
	ANL4DVector qj12(pj2e,pj2x,pj2y,pj2z);
	ANL4DVector qj21(pj3e,pj3x,pj3y,pj3z);
	ANL4DVector qj22(pj4e,pj4x,pj4y,pj4z);
	ANL4DVector qh1 = qj11 + qj12;
	ANL4DVector qh2 = qj21 + qj22;
	Double_t ph1  = qh1.Vect().Mag();
	Double_t ph2  = qh2.Vect().Mag();
	Double_t eh1p = TMath::Sqrt(ph1*ph1+kMh*kMh);
	Double_t eh2p = TMath::Sqrt(ph2*ph2+kMh*kMh);

        hMhMh->Fill(mh1, mh2, 1.);
        hMz  ->Fill(mz, 1.);
	hEvis->Fill(evis, 1.);
	hElmx->Fill(elmax, 1.);
	hCjmx->Fill(csjmax, 1.);
	hChCh->Fill(csh1, csh2, 1.);
	hCosH->Fill(csh1, 1.);
	hCosH->Fill(csh2, 1.);
#if 1
	hEh  ->Fill(eh1 , 1.);
	hEh  ->Fill(eh2 , 1.);
#else
	hEh  ->Fill(eh1p, 1.);
	hEh  ->Fill(eh2p, 1.);
#endif
	hPh  ->Fill(ph1 , 1.);
	hPh  ->Fill(ph2 , 1.);
	hMM  ->Fill(mm  , 1.);
	hPtPl->Fill(pt, pl, 1.);
	hCsjh->Fill(csj11h, 1.);
	hCsjh->Fill(csj12h, 1.);
	hCsjh->Fill(csj21h, 1.);
	hCsjh->Fill(csj22h, 1.);
	hCsjz->Fill(csj31h, 1.);
	hCsjz->Fill(csj32h, 1.);
	hFijh->Fill(fij11h, 1.);
	hFijh->Fill(fij12h, 1.);
	hFijh->Fill(fij21h, 1.);
	hFijh->Fill(fij22h, 1.);
	hFijz->Fill(fij31h, 1.);
	hFijz->Fill(fij32h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hPtPl->Draw();
  id++; c1->cd(id); hElmx->Draw();
  id++; c1->cd(id); hMhMh->Draw();
  id++; c1->cd(id); hMhMh->ProjectionX()->Draw();
  id++; c1->cd(id); hMhMh->ProjectionY()->Draw();
  id++; c1->cd(id); hMz ->Draw();
  id++; c1->cd(id); hEh ->Draw();
  id++; c1->cd(id); hPh ->Draw();
  id++; c1->cd(id); hMM ->Draw();
  //id++; c1->cd(id); hChCh->Draw();
  hCosH->SetMinimum(0.);
  id++; c1->cd(id); hCosH->Draw();
  hCsjh->SetMinimum(0.);
  id++; c1->cd(id); hCsjh->Draw();
  hFijh->SetMinimum(0.);
  id++; c1->cd(id); hFijh->Draw();
  hCsjz->SetMinimum(0.);
  id++; c1->cd(id); hCsjz->Draw();
  hFijz->SetMinimum(0.);
  id++; c1->cd(id); hFijz->Draw();
}
