{
// Macro example for bases calculation.

  gROOT->Reset();

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  TFile file(jsf->Env()->GetValue("JSFGUI.Spring.BasesFile","bases.root"),"RECREATE");

  bases = new XN2XN2Bases();

//  bases->SetNoOfSample(5000);
//  bases->SetIteration1( 0.2, 1);
//  bases->SetIteration2( 0.1, 1);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();

  bases->Write();
  file->Write();

}



