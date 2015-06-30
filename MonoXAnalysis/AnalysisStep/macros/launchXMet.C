{

  // Load the XMetAnalysis class (TSelector)
  gROOT->ProcessLine(".L XMetAnalysis_C.so");

  // Instantiate an XMetAnalysis with name of sub-directory in directory "plots/"
  XMetAnalysis x("test"); // will put results in "./plots/test/"

  // Launch the QCDKiller study
  // produce distribution of 5 qcd killers in all jet multiplicity bins
  x.StudyQCDKiller();

  // Launch the analysis
  // produce PFMETNoMu distribution from all processes
  //x.Analysis();     

  gSystem->Exit(0);

}
