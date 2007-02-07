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
//*	2007/02/04	K.Fujii		Original version.
//*
//*----------------------------------------------------------------------------

TH2D *hMxEa      = new TH2D("hMxEa"    , "", 150,   75.,  165.,
                                             150,  160.,  250.); 
TH2D *hCosXCosA  = new TH2D("hCosXCosA", "", 100,   -1.,   +1.,
                                             100,   -1.,   +1.); 
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);
TH1D *hPt        = new TH1D("hPt"      , "",  50,    0.,   50.);
TH1D *hPl        = new TH1D("hPl"      , "", 100,  -50.,   50.);
TH1D *hCosAh     = new TH1D("hCosAh"   , "", 100,   -1.,   +1.);
TH1D *hPhiAh     = new TH1D("hPhiAh"   , "", 180, -180., +180.);

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
#if 0
        Int_t event = elist->GetEntry(i);
        cerr << "i = " << i << " event = " << event << endl;
#else
        Int_t event = i;
#endif
        Double_t chi2;
        Double_t evis;
        Double_t pt;
        Double_t pl;
        Double_t cosx;
        Double_t phix;
        Double_t cosa;
        Double_t phia;
        Double_t mx;
        Double_t ex;
        Double_t ea1;
        Double_t pa1x;
        Double_t pa1y;
        Double_t pa1z;
        Double_t ea2;
        Double_t pa2x;
        Double_t pa2y;
        Double_t pa2z;
        Double_t ea;
        Double_t pax;
        Double_t pay;
        Double_t paz;
        Double_t cosa1h;
        Double_t phia1h;
        Double_t cosa2h;
        Double_t phia2h;

        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("cosx",&cosx);
        tup->SetBranchAddress("phix",&phix);
        tup->SetBranchAddress("cosa",&cosa);
        tup->SetBranchAddress("phia",&phia);
        tup->SetBranchAddress("mx",&mx);
        tup->SetBranchAddress("ex",&ex);
        tup->SetBranchAddress("ea1",&ea);
        tup->SetBranchAddress("pa1x",&pax);
        tup->SetBranchAddress("pa1y",&pay);
        tup->SetBranchAddress("pa1z",&paz);
        tup->SetBranchAddress("ea2",&ea);
        tup->SetBranchAddress("pa2x",&pax);
        tup->SetBranchAddress("pa2y",&pay);
        tup->SetBranchAddress("pa2z",&paz);
        tup->SetBranchAddress("ea",&ea);
        tup->SetBranchAddress("pax",&pax);
        tup->SetBranchAddress("pay",&pay);
        tup->SetBranchAddress("paz",&paz);
        tup->SetBranchAddress("cosa1h",&cosa1h);
        tup->SetBranchAddress("phia1h",&phia1h);
        tup->SetBranchAddress("cosa2h",&cosa2h);
        tup->SetBranchAddress("phia2h",&phia2h);

        tup->GetEntry(event);

        phia1h *= kToDeg;
        phia2h *= kToDeg;
        hMxEa    ->Fill(mx, ea, 1.);
        hCosXCosA->Fill(cosx, cosa, 1.);
        hEvis    ->Fill(evis  , 1.);
        hPt      ->Fill(pt    , 1.);
        hPl      ->Fill(pl    , 1.);
        hCosAh   ->Fill(cosa1h, 1.);
        hCosAh   ->Fill(cosa2h, 1.);
        hPhiAh   ->Fill(phia1h, 1.);
        hPhiAh   ->Fill(phia2h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hMxEa->Draw();
  id++; c1->cd(id); hMxEa->ProjectionX();
                    hMxEa_px->Draw();
  id++; c1->cd(id); hMxEa->ProjectionY();
                    hMxEa_py->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hPt->Draw();
  id++; c1->cd(id); hPl->Draw();
  id++; c1->cd(id); hCosXCosA->Draw();
  id++; c1->cd(id); hCosXCosA->ProjectionX();
                    hCosXCosA_px->SetMinimum(0.);
                    hCosXCosA_px->Draw();
  id++; c1->cd(id); hCosXCosA->ProjectionY();
                    hCosXCosA_py->SetMinimum(0.);
                    hCosXCosA_py->Draw();
  id++; c1->cd(id); hCosAh->SetMinimum(0.);
                    hCosAh->Draw();
  id++; c1->cd(id); hPhiAh->SetMinimum(0.);
                    hPhiAh->Draw();
}
