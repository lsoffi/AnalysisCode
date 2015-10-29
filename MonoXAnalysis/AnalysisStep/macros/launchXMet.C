#include "XMetAnalysis.h"

Int_t launchXMet(TString tag, TString dir)
{

  XMetAnalysis x(tag,dir);

  //x.CheckForwardJets(true);
  //x.StudyQCDKiller();
  //x.Analysis();     

  x.QCDScaleFactor(true);

  return 0;
}
