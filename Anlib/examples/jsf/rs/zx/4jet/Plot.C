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

TH2D *hMxMz      = new TH2D("hMxMz"    , "", 120,   90.,  150.,
                                             120,   60.,  120.);
TH2D *hCosXCosZ  = new TH2D("hCosXCosZ", "", 100,   -1.,   +1.,
                                             100,   -1.,   +1.); 
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);
TH1D *hCosJh     = new TH1D("hCosJh"   , "", 100,   -1.,   +1.);
TH1D *hPhiJh     = new TH1D("hPhiJh"   , "", 180, -180., +180.);
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
        Double_t ntracks;
        Double_t evis;
        Double_t pt;
        Double_t pl;
        Double_t ycut;
        Double_t chi2;
        Double_t npg1;
        Double_t eg1;
        Double_t pg1x;
        Double_t pg1y;
        Double_t pg1z;
        Double_t npg2;
        Double_t eg2;
        Double_t pg2x;
        Double_t pg2y;
        Double_t pg2z;
        Double_t npj1;
        Double_t ej1;
        Double_t pj1x;
        Double_t pj1y;
        Double_t pj1z;
        Double_t npj2;
        Double_t ej2;
        Double_t pj2x;
        Double_t pj2y;
        Double_t pj2z;
        Double_t cosg1h;
        Double_t phig1h;
        Double_t cosg2h;
        Double_t phig2h;
        Double_t cosj1h;
        Double_t phij1h;
        Double_t cosj2h;
        Double_t phij2h;
        Double_t xmass;
        Double_t zmass;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("npg1",&npg1);
        tup->SetBranchAddress("eg1",&eg1);
        tup->SetBranchAddress("pg1x",&pg1x);
        tup->SetBranchAddress("pg1y",&pg1y);
        tup->SetBranchAddress("pg1z",&pg1z);
        tup->SetBranchAddress("npg2",&npg2);
        tup->SetBranchAddress("eg2",&eg2);
        tup->SetBranchAddress("pg2x",&pg2x);
        tup->SetBranchAddress("pg2y",&pg2y);
        tup->SetBranchAddress("pg2z",&pg2z);
        tup->SetBranchAddress("npj1",&npj1);
        tup->SetBranchAddress("ej1",&ej1);
        tup->SetBranchAddress("pj1x",&pj1x);
        tup->SetBranchAddress("pj1y",&pj1y);
        tup->SetBranchAddress("pj1z",&pj1z);
        tup->SetBranchAddress("npj2",&npj2);
        tup->SetBranchAddress("ej2",&ej2);
        tup->SetBranchAddress("pj2x",&pj2x);
        tup->SetBranchAddress("pj2y",&pj2y);
        tup->SetBranchAddress("pj2z",&pj2z);
        tup->SetBranchAddress("cosg1h",&cosg1h);
        tup->SetBranchAddress("phig1h",&phig1h);
        tup->SetBranchAddress("cosg2h",&cosg2h);
        tup->SetBranchAddress("phig2h",&phig2h);
        tup->SetBranchAddress("cosj1h",&cosj1h);
        tup->SetBranchAddress("phij1h",&phij1h);
        tup->SetBranchAddress("cosj2h",&cosj2h);
        tup->SetBranchAddress("phij2h",&phij2h);
        tup->SetBranchAddress("xmass",&xmass);
        tup->SetBranchAddress("zmass",&zmass);

        tup->GetEntry(event);
#if 0
  cerr << " pg1 = (" << eg1  << ", "
                     << pg1x << ", "
                     << pg1y << ", "
                     << pg1z << ") " << endl;
#endif
        ANL4DVector pg1(eg1, pg1x, pg1y, pg1z);
        ANL4DVector pg2(eg2, pg2x, pg2y, pg2z);
        ANL4DVector pj1(ej1, pj1x, pj1y, pj1z);
        ANL4DVector pj2(ej2, pj2x, pj2y, pj2z);
        ANL4DVector px = pg1 + pg2;
        ANL4DVector pz = pj1 + pj2;
#if 0 
  cerr << " xmass = " << px.GetMass() << " zmass = " << pz.GetMass() << endl;
#endif

        phij1h *= kToDeg;
        phij2h *= kToDeg;
        phig1h *= kToDeg;
        phig2h *= kToDeg;
        hMxMz    ->Fill(xmass, zmass, 1.);
        hCosXCosZ->Fill(px.CosTheta(), pz.CosTheta(), 1.);
        hEvis    ->Fill(evis  , 1.);
        hCosJh   ->Fill(cosj1h, 1.);
        hCosJh   ->Fill(cosj2h, 1.);
        hPhiJh   ->Fill(phij1h, 1.);
        hPhiJh   ->Fill(phij2h, 1.);
        hCosAh   ->Fill(cosg1h, 1.);
        hCosAh   ->Fill(cosg2h, 1.);
        hPhiAh   ->Fill(phig1h, 1.);
        hPhiAh   ->Fill(phig2h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hMxMz->Draw();
  id++; c1->cd(id); hMxMz->ProjectionX()->Draw();
  id++; c1->cd(id); hMxMz->ProjectionY()->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hCosXCosZ->Draw();
  id++; c1->cd(id); hCosXCosZ->ProjectionX();
                    hCosXCosZ_px->SetMinimum(0.);
                    hCosXCosZ_px->Draw();
  id++; c1->cd(id); hCosXCosZ->ProjectionY();
                    hCosXCosZ_py->SetMinimum(0.);
                    hCosXCosZ_py->Draw();
  id++; c1->cd(id); hCosJh->SetMinimum(0.); hCosJh->Draw();
  id++; c1->cd(id); hPhiJh->SetMinimum(0.); hPhiJh->Draw();
  id++; c1->cd(id); hCosAh->SetMinimum(0.); hCosAh->Draw();
  id++; c1->cd(id); hPhiAh->SetMinimum(0.); hPhiAh->Draw();
}
