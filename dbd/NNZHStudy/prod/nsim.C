#include <iomanip>
void nsim(Int_t irun)
{ 
  gROOT->Reset();

  stringstream sfname;
  sfname << "./data/nnzhsim" << setw(3) << setfill('0') << irun << ".root" << ends;
  //  sfname << "./data/nnzhsim_02_" << setw(3) << setfill('0') << irun << ".root" << ends;
  TFile file(sfname.str().data(),"RECREATE");  // Output file
  cerr << "-------------------------------" << endl;
  cerr << " Opened " << sfname.str().data() << endl;

  JSFSteer jsf;
  JSFLCFULL full;
  NNZHSpring spring("NNZHSpring");

  stringstream slname;
  if (irun == 1) {
    slname << "./data/nnzhsim.root" << ends;
    //    slname << "./data/nnzhsim100.root" << ends;
    //    slname << "./data/nnzhsim200.root" << ends;
  } else {
    slname << "./data/nnzhsim" << setw(3) << setfill('0') << irun-1 << ".root" << ends;
    //    slname << "./data/nnzhsim_02_" << setw(3) << setfill('0') << irun-1 << ".root" << ends;
  }
  cerr << " Read "   << slname.str().data() << endl;
  cerr << "-------------------------------" << endl;
  TFile flast(slname.str().data(),"READ");  // Last run file
  jsf.GetLastRunInfo(&flast);
  flast.Close();

  JSFHadronizer hdr;
  //  JSFQuickSim sim;
  stdhep = new JSFWriteStdHep();

  spring.ReadBases("bases.root");

  Int_t maxevt=10100;      // Number of events.
  //  Int_t maxevt=100;        // Number of events.
  jsf.Initialize();

  printf(" Roots is %g\n",((NNZHBases*)spring.GetBases())->GetRoots());
  jsf.BeginRun(30+irun);      // Set run number to 30.

  Int_t ev = 1;
  while (1) {
    printf(" start event %d ",ev);
    if (jsf.Process(ev)) {
      printf("Processed event %d ",ev);

      jsf.FillTree();
      jsf.Clear();

      printf(" End event %d \n",ev);
      if (ev >= maxevt) break;
      ev++;
    }
  }
  jsf.Terminate();
  //file.Write();
}
