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
  XMetProcess(TString nameProcess, Int_t col, TString nameFile, 
	      Double_t xsec, Int_t ngen);
  ~XMetProcess();

  // Adders
  Int_t AddDir(TString dir);
  Int_t AddTrees();

  // Setters
  Int_t SetNameTree(TString nameTree);
  Int_t SetNameFile(TString nameFile);
  Int_t SetPath(TString path);
  Int_t SetXSec(Double_t xsec);
  Int_t SetWeight(Double_t w);
  Int_t SetNGen( Int_t ngen);
  Int_t SetColor(Int_t color);
  Int_t SetStyle(Int_t style);
  Int_t SetSize( Int_t size);

  // Getters
  Int_t    GetColor();
  Int_t    GetStyle();
  Int_t    GetSize();
  TString  GetNameFile();
  Double_t GetXSec();
  Double_t GetWeight();
  Int_t    GetNGen();

  // Analysis tools
  Int_t Skim(TString select, TCut cut, TString reset);
  Int_t Draw(TH1* h, TString var, TCut cut, TCut weight);
    
 private:
  TChain* _chain;
  vector<TString> _dir;

  TString  _path, _nameFile, _nameProcess;
  Double_t _xsec, _weight;
  Int_t    _ngen, _col, _style, _size;

  map<TString, TEntryList*> _mapSkim;
  map<TString, TEntryList*>::iterator _itmapSkim;
};


// Destructor //
XMetProcess::~XMetProcess()
{
  /*
  delete _chain;
  for(_itmapSkim=_mapSkim.begin();_itmapSkim!=_mapSkim.end();_itmapSkim++) {
    delete (_itmapSkim->second);
  }
  */
}

// Constructors //
XMetProcess::XMetProcess()
{
  _nameProcess = "";
  _path = "";
  _nameFile = "";
  _col   = kBlack;
  _style = kFullCircle;
  _size  = 1.0;
  _chain = new TChain();
  _xsec=1;
  _weight=1;
  _ngen=1;
}

XMetProcess::XMetProcess(TString nameProcess, Int_t col)
{
  _nameProcess = nameProcess;
  _path = "";
  _nameFile = "";
  _col = col;
  _chain = new TChain();
  _xsec=1;
  _weight=1;
  _ngen=1;
}

XMetProcess::XMetProcess(TString nameProcess, Int_t col, TString nameFile)
{
  _nameProcess = nameProcess;
  _path = "";
  _nameFile = nameFile;
  _col = col;
  _chain = new TChain();
  _xsec=1;
  _weight=1;
  _ngen=1;
}

XMetProcess::XMetProcess(TString nameProcess, Int_t col, TString nameFile, 
			 Double_t xsec, Int_t ngen)
{
  _nameProcess = nameProcess;
  _path = "";
  _nameFile = nameFile;
  _col = col;
  _chain = new TChain();
  _xsec=xsec;
  _ngen=ngen;

  _weight = ngen!=0 ? (1000.0*xsec)/ngen : 1.0;
}

// Adders //
Int_t XMetProcess::AddDir(TString dir)
{
  _dir.push_back(dir);
  return 0;
}

Int_t XMetProcess::AddTrees()
{
  for(UInt_t iD=0 ; iD<_dir.size() ; iD++) {
    _chain->Add(_path+"/"+_dir[iD]+"/"+_nameFile);
  }

  if(_chain->IsZombie()) {
    cout << "ERROR: IsZombie() "
	 << _nameProcess << endl;
    return -1;
  }

  return 0;
}

// Getters //
Double_t XMetProcess::GetXSec()
{
  return _xsec;
}

Double_t XMetProcess::GetWeight()
{
  return _weight;
}

Int_t XMetProcess::GetNGen()
{
  return _ngen;
}

Int_t XMetProcess::GetColor()
{
  return _col;
}

Int_t XMetProcess::GetStyle()
{
  return _style;
}

Int_t XMetProcess::GetSize()
{
  return _size;
}

TString XMetProcess::GetNameFile()
{
  return _nameFile;
}

// Setters //
Int_t XMetProcess::SetNameTree(TString nameTree)
{
  _chain->SetName(nameTree);
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

Int_t XMetProcess::SetXSec(Double_t xsec)
{
  _xsec = xsec;
  return 0;
}

Int_t XMetProcess::SetWeight(Double_t w)
{
  _weight = w;
  return 0;
}

Int_t XMetProcess::SetNGen(Int_t ngen)
{
  _ngen=ngen;
  return 0;
}

Int_t XMetProcess::SetColor(Int_t color)
{
  _col = color;
  return 0;
}

Int_t XMetProcess::SetStyle(Int_t style)
{
  _style = style;
  return 0;
}

Int_t XMetProcess::SetSize(Int_t size)
{
  _size = size;
  return 0;
}

// Analysis tools //
Int_t XMetProcess::Skim(TString select, TCut cut, TString option)
{

  TString tskim="skim_"+_nameProcess+"_"+select;
  TEntryList* skim; 

  if(option.Contains("Reset")) {
    _chain->SetEntryList(0);
  }

  if(option.Contains("Produce")) {
    _chain->Draw(">>+"+tskim, cut, "entrylist");
    //_chain->Draw(">>+"+tskim, cut, "entrylist",1000); // FIXME
  }

  skim = (TEntryList*)gDirectory->Get(tskim);
  _mapSkim[select] = skim;
  _chain->SetEntryList(skim);

  cout << "--- produced skim : " << tskim 
       << " : " << skim->GetN() << " entries" 
       << endl;

  return 0;
}

Int_t XMetProcess::Draw(TH1* h, TString var, TCut cut, TCut weight)
{
  _chain->Draw(var+">>"+TString(h->GetName()), cut*weight);
  //_chain->Draw(var+">>"+TString(h->GetName()), cut*weight, "", 1000); // FIXME
  return 0;
}

#endif
