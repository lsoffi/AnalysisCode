{

  gROOT->ProcessLine(".L XMetAnalysis_C.so");
  gROOT->ProcessLine(".L launchXMet_C.so");

  //launchXMet("v6_AN15_Spring15MC25ns_FWD_Met50_FixJMCut_NoWgtMC_NoJetCounters_Nadir22Oct","QCD");

  launchXMet("v15_AN15_Met50_Data05Oct_Lumi210_MaxRun257599_QCDSF_JetHT_FixOfflineSel","QCD");

  gSystem->Exit(0);

}
