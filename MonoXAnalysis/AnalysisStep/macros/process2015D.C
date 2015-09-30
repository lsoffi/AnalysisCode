{

  gROOT->ProcessLine(".L MyTrigger_C.so");
  
  MyTrigger m("v11_SingleMuPrompt2015D_25ns_ETM50_Prompt_HLTMuTightAna",
	      "HLTMu_TightMuon_Ana", 
	      "2015D", "25ns", "ETM50", "Prompt", 
	      "38T", "noskim", "NoHBHE", "tune");
  
  m.ProdHistos();
  
  gSystem->Exit(0);

}
