{

  gROOT->ProcessLine(".L XMetAnalysis_C.so");
  gROOT->ProcessLine(".L launchXMet_C.so");

  launchXMet("v6_AN15_Spring15MC25ns_FWD_Met50_FixJMCut_NoWgtMC_NoJetCounters_Nadir22Oct","QCD");

  gSystem->Exit(0);

}
