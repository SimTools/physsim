//****************************************************
//*
//*  Dummy UserAnalysis Script to be used with gui.C 
//*  Functions in this scripts are called by GUIMainMacro.C  
//*  
//*  In this script, global functions, UserInitialize() nad UserAnalysis()
//*  Must be defined.  UserInitialize() is called at job initialization 
//*  to initialize histograms, etc.  It is also called when reset historgam
//*  menu is selected from JSF Control Panel.  UserAnalysis() is called 
//*  during the event loop after executing Process functions of JSF modules
//*  and display event daia but before display histogram.  See GetEvent()
//*  function in GUIMainMacro.C, where UserAnalysis() is called.
//*
//*  UserSetOptions() and UserTerminate() may be defined in this file.
//*  UserSetOptions() are called before declaration of JSF modules.
//*  It can be used to set parameters optional for user analsis.
//*  UserTerminate() is called during the JSF termination process.
//*
//*  When runmode=4, UserModuleDefine() must be defined in this file.
//*  It is used to define JSF modules specific to user analysis.  
//*  
//*$Id$
//*
//****************************************************


//_________________________________________________________
void UserModuleDefine()
{
    Char_t *inputFileName=jsf->Env()->GetValue("JSFGUI.InputFileName","");
    Char_t *outputFileName=jsf->Env()->GetValue("JSFGUI.OutputFileName",
					      "jsfj4.root");

    ofile= new TFile(outputFileName,"RECREATE");
    file = new TFile(inputFileName);  // Input file
    jsf->SetIOFiles();
    jsf->SetOutput(*ofile);

    
    TString liblist(gSystem->GetLibraries("","D"));
    if( liblist.Contains("libJSFJupiter.so") ) {
      jsfj4 = new JSFJupiter();
    }
    gJSFJ4 = jsfj4;

}

//_________________________________________________________
void UserInitialize()
{
  //  This function is called at the begining of the job or when
  //  "reset hist" action is selected in the gui menu.
  //  This is used to define/reset histograms.

}


//_________________________________________________________
void UserAnalysis()
{
  std::cerr << "Read event number=" << jsf->GetEventNumber() << std::endl;

}

//_________________________________________________________
void DrawHist()
{
  //  This function is called to draw histograms during the interactive 
  //  session.  Thus you can see the accumulation of the histogram
  //  interactively.  

}
