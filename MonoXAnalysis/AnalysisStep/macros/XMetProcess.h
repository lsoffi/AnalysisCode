#ifndef XMETPROCESS
#define XMETPROCESS

#include "myIncludes.h"

using namespace std;

class XMetProcess {

 public:

  // Constructors, destructor
  XMetProcess();
  XMetProcess(TString nameProcess, Int_t col);
  XMetProcess(TString nameProcess, Int_t col, TString nameFile);
  ~XMetProcess();

  // Setters
  Int_t AddDir(TString dir);
  Int_t SetNameTree(TString nameTree);
  Int_t SetNameFile(TString nameFile);
  Int_t SetPath(TString path);
  Int_t AddTrees();

  // Getters
  Int_t GetColor();

  // Analysis tools
  Int_t Skim(TString select, TCut cut);
  Int_t Draw(TH1F* h, TString var, TCut cut, TCut weight);
    
 private:
  TChain _chain;
  vector<TString> _dir;
  Int_t _col;

  TString _path;
  TString _nameFile;
  TString _nameProcess;

  map<TString, TEntryList*> _mapSkim;
};

XMetProcess::~XMetProcess()
{
}

XMetProcess::XMetProcess()
{
  _nameProcess = "";
  _path = "";
  _nameFile = "";
  _col = kBlack;
}

XMetProcess::XMetProcess(TString nameProcess, Int_t col)
{
  _nameProcess = nameProcess;
  _col = col;
}

XMetProcess::XMetProcess(TString nameProcess, Int_t col, TString nameFile)
{
  _nameProcess = nameProcess;
  _nameFile = nameFile;
  _col = col;
}


Int_t XMetProcess::AddDir(TString dir)
{
  _dir.push_back(dir);
  return 0;
}

Int_t XMetProcess::GetColor()
{
  return _col;
}

Int_t XMetProcess::SetNameTree(TString nameTree)
{
  _chain.SetName(nameTree);
  return 0;
}

Int_t XMetProcess::SetPath(TString path)
{
  _path = path;
  return 0;
}

Int_t XMetProcess::SetNameFile(TString nameFile)
{
  _nameFile = nameFile;
  return 0;
}

Int_t XMetProcess::AddTrees()
{
  for(UInt_t iD=0 ; iD<_dir.size() ; iD++) {
    _chain.Add(_path+"/"+_dir[iD]+"/"+_nameFile);
  }
  return 0;
}

Int_t XMetProcess::Skim(TString select, TCut cut)
{
  _chain.SetEntryList(0);

  TString tskim="skim_"+_nameProcess+"_"+select;
  _chain.Draw(">>+"+tskim, cut, "entrylist");
  TEntryList* skim = (TEntryList*)gDirectory->Get(tskim);

  _mapSkim[select] = skim;
  _chain.SetEntryList(skim);

  cout << "--- produced skim : " << tskim 
       << " : " << skim->GetN() << " entries" 
       << endl;

  return 0;
}

Int_t XMetProcess::Draw(TH1F* h, TString var, TCut cut, TCut weight)
{
  _chain.Draw(var+">>"+TString(h->GetName()), cut*weight);
  return 0;
}

#endif
