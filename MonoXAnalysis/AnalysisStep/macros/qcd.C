#include "myIncludes.h"

bool FAST=true;

using namespace std;
typedef map<TString, map<TString, map<TString, TH1F*> > > M_SEL_PROC_VAR_H1D;
typedef map<TString, map<TString, map<TString, TH2F*> > > M_SEL_PROC_VAR_H2D;

typedef map<TString,TChain*> M_PROC_CH;

typedef map<TString, TH1F*>                M_VAR_H1D;
typedef map<TString, M_VAR_H1D >           M_CUT_VAR_H1D;
typedef map<TString, M_CUT_VAR_H1D >       M_PROC_CUT_VAR_H1D;
typedef map<TString, M_PROC_CUT_VAR_H1D >  M_SEL_PROC_CUT_VAR_H1D;

typedef map<TString, TH2F*>                M_VAR_H2D;
typedef map<TString, M_VAR_H2D >           M_CUT_VAR_H2D;
typedef map<TString, M_CUT_VAR_H2D >       M_PROC_CUT_VAR_H2D;
typedef map<TString, M_PROC_CUT_VAR_H2D >  M_SEL_PROC_CUT_VAR_H2D;

// Declare functions //
Int_t defineChains(M_PROC_CH &chains);
TCut  defineWeight(TString process);
TCut  defineCut(TString sample, TString region);
Int_t defineHistos1D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H1D &histos1D,
		     const UInt_t nS,    const UInt_t nC,
		     TString* selection, TString* extra);
Int_t defineHistos2D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H2D &histos2D,
		     const UInt_t nS,    const UInt_t nC,
		     TString* selection, TString* extra);
Int_t analyze(Float_t lumi, TString tag, M_PROC_CH chains, M_SEL_PROC_CUT_VAR_H1D &histos1D, M_SEL_PROC_CUT_VAR_H2D &histos2D,
	      const UInt_t nS,    const UInt_t nC, TString* selection, TString* extra);
vector<TString> GetVarName(M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D);
Int_t GetNGen(TString process);
Int_t writeHistos(TString tag, M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D);

// MAIN //
Int_t qcd(TString tag)
{

  Float_t lumi=0.210;
  
  TString dirOut="/user/ndaci/Results/Monojet/QCD/"+tag+"/";
  ofstream outlog(dirOut+"yields_"+tag+".txt",ios::out);

  // Define Chains
  M_PROC_CH chains;
  defineChains(chains);

  // Define selections to be explored
  // const UInt_t nS=5;
  // TString selection[nS] = {"inclusive","1jet","2jet","3jet","4jet"};
  const UInt_t nS=2;
  TString selection[nS] = {"inclusive","1jet"};

  const UInt_t nC=2;
  TString extra[nC]={"None","FwdVeto"};

  // Define histograms
  M_SEL_PROC_CUT_VAR_H1D histos1D;
  if(tag.Contains("1D")) defineHistos1D(chains, histos1D, nS, nC, selection, extra);

  M_SEL_PROC_CUT_VAR_H2D histos2D;
  if(tag.Contains("2D")) defineHistos2D(chains, histos2D, nS, nC, selection, extra);

  // 2D plots: region="HT"
  // 1D stack: region="Met200_1jet"
  analyze(lumi, tag, chains, histos1D, histos2D, nS, nC, selection, extra);

  writeHistos(tag, histos1D, histos2D);

  //doStackPlots(histos1D, nS, nC, selection, extra);
  //doTransferPlots(histos2D, nS, nC, selection, extra);

  return 0;
}

Int_t writeHistos(TString tag, M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D)
{

  TString dirOut="/user/ndaci/Results/Monojet/QCD/";
  TFile* outfile = new TFile(dirOut+"/"+tag+"/plots_"+tag+".root","recreate");
  outfile->cd();

  M_SEL_PROC_CUT_VAR_H1D::iterator itSelProcCutVarH1D;
  M_PROC_CUT_VAR_H1D::iterator     itProcCutVarH1D;
  M_CUT_VAR_H1D::iterator          itCutVarH1D;
  M_VAR_H1D::iterator              itVarH1D;

  M_SEL_PROC_CUT_VAR_H2D::iterator itSelProcCutVarH2D;
  M_PROC_CUT_VAR_H2D::iterator     itProcCutVarH2D;
  M_CUT_VAR_H2D::iterator          itCutVarH2D;
  M_VAR_H2D::iterator              itVarH2D;

  for((itSelProcCutVarH1D=histos1D.begin()) ; 
      (itSelProcCutVarH1D!=histos1D.end()) ; 
      itSelProcCutVarH1D++) {
    for((itProcCutVarH1D=(itSelProcCutVarH1D->second).begin()) ; 
	(itProcCutVarH1D!=(itSelProcCutVarH1D->second).end()) ; 
	itProcCutVarH1D++) {
      for((itCutVarH1D=(itProcCutVarH1D->second).begin()) ; 
	  (itCutVarH1D!=(itProcCutVarH1D->second).end()) ; 
	  itCutVarH1D++) {
	for((itVarH1D=(itCutVarH1D->second).begin()) ; 
	    (itVarH1D!=(itCutVarH1D->second).end()) ; 
	    itVarH1D++) {
	  (itVarH1D->second)->Write();
	}
      }
    }
  }

  for((itSelProcCutVarH2D=histos2D.begin()) ; 
      (itSelProcCutVarH2D!=histos2D.end()) ; 
      itSelProcCutVarH2D++) {
    for((itProcCutVarH2D=(itSelProcCutVarH2D->second).begin()) ; 
	(itProcCutVarH2D!=(itSelProcCutVarH2D->second).end()) ; 
	itProcCutVarH2D++) {
      for((itCutVarH2D=(itProcCutVarH2D->second).begin()) ; 
	  (itCutVarH2D!=(itProcCutVarH2D->second).end()) ; 
	  itCutVarH2D++) {
	for((itVarH2D=(itCutVarH2D->second).begin()) ; 
	    (itVarH2D!=(itCutVarH2D->second).end()) ; 
	    itVarH2D++) {
	  (itVarH2D->second)->Write();
	}
      }
    }
  }

  outfile->Write();
  return 0;
}

Int_t analyze(Float_t lumi, TString tag, M_PROC_CH chains, M_SEL_PROC_CUT_VAR_H1D &histos1D, M_SEL_PROC_CUT_VAR_H2D &histos2D,
	      const UInt_t nS,    const UInt_t nC, TString* selection, TString* extra)
{

  vector<TString> var = GetVarName(histos1D, histos2D);
  const UInt_t nV=var.size();

  // Loop over chains
  TChain* ch;
  TEntryList* skim;
  TH1F* hTemp1;
  TH2F* hTemp2;
  TString process, tskim, select, hname, locVar;
  M_PROC_CH::iterator itCh;
  TCut theCut, theWeight;
  Float_t theScale;
  Int_t   nGen;
  //
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {

    process = itCh->first;
    ch      = itCh->second;

    cout << "- " << process << endl;

    theScale = 1;
    nGen     = GetNGen(process);
    if(!process.Contains("data")) {
      theScale = nGen!=0 ? lumi/nGen : 1;
    }
    if(process.Contains("zll") || process.Contains("znn")) {
      theScale *= 1.23;
    }
    if(process.Contains("wln")) {
      theScale *= 1.21;
    }

    for(UInt_t iS=0 ; iS<nS ; iS++) {

      cout << "-- " << selection[iS] << endl;

      for(UInt_t iC=0 ; iC<nC ; iC++) {

	cout << "--- " << extra[iC] << endl;
	
	// Skim the chain
	select = selection[iS]+"_"+extra[iC];
	theCut = defineCut(process, select);
	theWeight = defineWeight(process);
	tskim  = "skim_"+process+"_"+select;
	//
	ch->SetEntryList(0);
	if(FAST) ch->Draw(">>+"+tskim, theCut, "entrylist", 1);
	else     ch->Draw(">>+"+tskim, theCut, "entrylist");
	skim = (TEntryList*)gDirectory->Get(tskim);
	if(skim) ch->SetEntryList(skim);

	// Draw the variable
	for(UInt_t iV=0 ; iV<nV ; iV++) {

	  cout << "----- " << process << " " << selection[iS] << " " << extra[iC] << " " << var[iV] << endl;
	  
	  locVar=var[iV];
	  if(     var[iV]=="jetmetdphimin_vs_t1mumet")    locVar="jetmetdphimin:t1mumet";
	  else if(var[iV]=="incjetmetdphimin_vs_t1mumet") locVar="incjetmetdphimin:t1mumet";

	  if(var[iV].Contains("_vs_")) {
	    hTemp2 = histos2D[selection[iS]][process][extra[iC]][var[iV]];
	    hname  = hTemp2->GetName();
	    theCut = defineCut(process, select+"_HT");
	    ch->Draw(locVar+">>"+hname, theCut*theWeight);
	    hTemp2->Scale(theScale);
	  }
	  else {
	    hTemp1 = histos1D[selection[iS]][process][extra[iC]][var[iV]];
	    hname  = hTemp1->GetName();
	    theCut = defineCut(process, select+"_Met200");
	    ch->Draw(locVar+">>"+hname, theCut*theWeight);
	    hTemp1->Scale(theScale);
	  }

	} // end loop: var
      } // end loop: cuts
    } // end loop: selection
  } // end loop: processes

  return 0;
}

TCut defineWeight(TString process)
{

  TCut weight="1";

  if(!process.Contains("data")) {
    weight = "xsec*puweight*(wgt/wgtsum)";
  }

  return weight;
}

Int_t defineHistos1D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H1D &histos1D,
		     const UInt_t nS,    const UInt_t nC, 
		     TString* selection, TString* extra)
{

  /// Stack 1D plots
  const UInt_t nV=3;
  TString var[nV]     ={"t1mumet", "jetmetdphimin", "incjetmetdphimin"};
  TString nameAxisX[nV]={"Type1 PFMETNoMu [GeV]", "Min #Delta#phi(M,J_{i}^{C})", "Min #Delta#phi(M,J_{i})"};
  TString nameAxisY[nV]={"Events","Events","Events"};

  UInt_t   nBins[nV]={ 100,   64,   64};
  Float_t xFirst[nV]={ 200,    0,    0};
  Float_t xLast[ nV]={1000,  3.2,  3.2};

  vector<UInt_t>   v_nBins;
  vector<Float_t*> v_bins;

  Float_t bins_met[8] = {200, 250, 300, 350, 400, 500, 600, 1000};
  v_nBins.push_back(8);
  v_bins .push_back(bins_met);
  
  Float_t bins_phi[64];
  for(UInt_t iB=0 ; iB<64 ; iB++) {
    bins_phi[iB] = xFirst[1] + (iB*(xLast[1]-xFirst[1])/nBins[1]);
  }
  for(UInt_t i=1 ; i<nV ; i++) {
    v_nBins.push_back(64);
    v_bins .push_back(bins_phi);
  }

  TString name,title;

  // Get process names from chains map
  vector<TString> process;
  M_PROC_CH::iterator itCh;
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {
    process.push_back( itCh->first );
  }
  const UInt_t nP = process.size();

  // Create histos and fill the map
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      for(UInt_t iC=0 ; iC<nC ; iC++) {

	for(UInt_t iV=0 ; iV<nV ; iV++) { 
	  //for(UInt_t iV=1 ; iV<nV ; iV++) {

	  name ="h1D_"+selection[iS]+"_"+process[iP]+"_"+extra[iC]+"_"+var[iV];
	  title=selection[iS]+" "+process[iP]+" "+extra[iC]+" "+var[iV];

	  histos1D[selection[iS]][process[iP]][extra[iC]][var[iV]] = 
	    new TH1F(name,title,v_nBins[iV]-1,v_bins[iV]);
	}
      }
    }
  }

  return 0;
}

Int_t defineHistos2D(M_PROC_CH chains,   M_SEL_PROC_CUT_VAR_H2D &histos2D,
		     const UInt_t nS,    const UInt_t nC, 
		     TString* selection, TString* extra)
{

  /// Stack 2D plots
  const UInt_t nV=2;
  TString var[nV]      ={"jetmetdphimin_vs_t1mumet", "incjetmetdphimin_vs_t1mumet"};
  TString nameAxisX[nV]={"Type1 PFMETNoMu [GeV]", "Type1 PFMETNoMu [GeV]"};
  TString nameAxisY[nV]={"Min #Delta#phi(M,J_{i}^{C})", "Min #Delta#phi(M,J_{i})"};

  UInt_t   nBinsY[nV]={64,   64};
  Float_t yFirst[ nV]={ 0,    0};
  Float_t yLast[  nV]={ 3.2,  3.2};

  vector<UInt_t>   v_nBinsX, v_nBinsY;
  vector<Float_t*> v_xbins,  v_ybins;

  const UInt_t nBinsMet=8;
  Float_t bins_met[nBinsMet] = {200, 250, 300, 350, 400, 500, 600, 1000};

  const UInt_t nBinsPhi=64;
  Float_t bins_phi[nBinsPhi];
  for(UInt_t iB=0 ; iB<nBinsPhi ; iB++) {
    bins_phi[iB] = yFirst[1] + (iB*(yLast[1]-yFirst[1])/nBinsY[1]);
  }

  for(UInt_t i=0 ; i<nV ; i++) {
    v_nBinsX.push_back(nBinsMet);
    v_xbins .push_back(bins_met);
    v_nBinsY.push_back(nBinsPhi);
    v_ybins .push_back(bins_phi);
  }

  TString name,title;

  // Get process names from chains map
  vector<TString> process;
  M_PROC_CH::iterator itCh;
  for(itCh = chains.begin() ; itCh!=chains.end() ; itCh++) {
    process.push_back( itCh->first );
  }
  const UInt_t nP = process.size();

  // Create histos and fill the map
  for(UInt_t iP=0 ; iP<nP ; iP++) {
    for(UInt_t iS=0 ; iS<nS ; iS++) {
      for(UInt_t iC=0 ; iC<nC ; iC++) {

	for(UInt_t iV=0 ; iV<nV ; iV++) {

	  name ="h2D_"+selection[iS]+"_"+process[iP]+"_"+extra[iC]+"_"+var[iV];
	  title=selection[iS]+" "+process[iP]+" "+extra[iC]+" "+var[iV];

	  histos2D[selection[iS]][process[iP]][extra[iC]][var[iV]] = 
	    new TH2F(name,title,v_nBinsX[iV]-1,v_xbins[iV],v_nBinsY[iV]-1,v_ybins[iV]);
	}
      }
    }
  }

  return 0;
}

Int_t defineChains(M_PROC_CH &chains)
{

  TString pathMC="/user/ndaci/Data/XMET/AdishTrees_29Oct2015/Spring15MC_25ns/";
  TString nameFileMC="skim.root";
  TString nameTreeMC= "tree/tree";

  TString pathData= "/user/ndaci/Data/XMET/AdishTrees_29Oct2015/Run2015D/";
  TString pathDataHT= "/user/ndaci/Data/XMET/NadirTrees_26Oct2015/Run2015D/";
  TString nameFileData="tree.root";
  TString nameTreeData= "tree/tree";

  const UInt_t nMC=28;
  TString dirMC[nMC]={"znn100to200","znn200to400","znn400to600","znn600toinf",
		      "wln100to200","wln200to400","wln400to600","wln600toinf",
		      "zll100to200","zll200to400","zll400to600","zll600toinf",
		      "ttbar","ww","wz","zz",
		      "singletbart","singletbarw","singlett","singletw",
		      "qcdht1000to1500","qcdht100to200","qcdht1500to2000","qcdht2000toinf",
		      "qcdht200to300","qcdht300to500","qcdht500to700","qcdht700to1000"};

  const UInt_t nData=2;
  TString nameData[nData]={"data_jetht","data_met"};
  TString dirData[nData] ={"jetht_05Oct2015","met"};

  // Define chains
  for(UInt_t i=0 ; i<nMC ; i++) {
    chains[dirMC[i]] = new TChain("tree/tree");
    chains[dirMC[i]]->Add(pathMC+"/"+dirMC[i]+"/"+nameFileMC);
  }
  //
  for(UInt_t i=0 ; i<nData ; i++) {
    chains[nameData[i]] = new TChain("tree/tree");
    if(nameData[i].Contains("ht"))
      chains[nameData[i]]->Add(pathDataHT+"/"+dirData[i]+"/skimJSON_t1mumet50.root");
    else
      chains[nameData[i]]->Add(pathData+"/"+dirData[i]+"/"+nameFileData);
  }

  return 0;
}

TCut defineCut(TString sample, TString region)
{

  TCut trig ="hltmet90>0 || hltmet120>0";
  TCut noise="flaghbheloose>0 && flagcsctight>0 && flageebadsc>0";
  TCut maxrun="run<257599";
  TCut metCut="";
  TCut leptons="nelectrons == 0 && ntaus == 0 && nmuons == 0";
  TCut photons="nphotons==0";
  TCut bveto="nbjets==0";
  TCut jetID="signaljetpt>100 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && signaljetNHfrac < 0.8";
  TCut jetBin="njets>=1";
  TCut noqcd="";
  TCut fwdveto="";

  if(region.Contains("FwdVeto")) {
    fwdveto="abs(leadingjeteta)<2.5";
    if(sample.Contains("jetht")) fwdveto="(signaljetpt>Max$(jet_pt * (abs(jet_eta)>2.5)))";
  }

  if(region.Contains("HT")) {
    if(sample.Contains("jetht")) trig = "hltpfht200 && ht>300";
    else                         trig = "ht>300";
  }  
  
  if(region.Contains("Met200")) metCut="t1mumet>200";

  // Jet ID //
  TCut jetID1 = "(signaljetpt>100 && abs(signaljeteta)<2.5 && signaljetCHfrac > 0.1 && signaljetNHfrac < 0.8)";
  TCut jetID2 = "(secondjetpt>30  && abs(secondjeteta)<2.5 && secondjetCHfrac > 0.1 && secondjetNHfrac < 0.8)";
  TCut jetID3 = "( thirdjetpt>30  && abs(thirdjeteta)<2.5  &&  thirdjetCHfrac > 0.1 &&  thirdjetNHfrac < 0.8)";

  if(     region.Contains("1jet"))   {jetBin="njets==1"; jetID = jetID1;}
  else if(region.Contains("2jet"))   {jetBin="njets==2"; jetID = jetID1*jetID2;}
  else if(region.Contains("3jet"))   {jetBin="njets==3"; jetID = jetID1*jetID2*jetID3;}
  else if(region.Contains("4jet"))   {jetBin="njets>=4"; jetID = jetID1*jetID2*jetID3;}

  return trig*noise*maxrun*metCut*leptons*photons*bveto*jetID*jetBin*noqcd*fwdveto;

}

vector<TString> GetVarName(M_SEL_PROC_CUT_VAR_H1D histos1D, M_SEL_PROC_CUT_VAR_H2D histos2D)
{

  vector<TString> var;

  M_SEL_PROC_CUT_VAR_H1D::iterator itSelProcCutVarH1D;
  M_PROC_CUT_VAR_H1D::iterator     itProcCutVarH1D;
  M_CUT_VAR_H1D::iterator          itCutVarH1D;
  M_VAR_H1D::iterator              itVarH1D;

  M_SEL_PROC_CUT_VAR_H2D::iterator itSelProcCutVarH2D;
  M_PROC_CUT_VAR_H2D::iterator     itProcCutVarH2D;
  M_CUT_VAR_H2D::iterator          itCutVarH2D;
  M_VAR_H2D::iterator              itVarH2D;

  for((itSelProcCutVarH1D=histos1D.begin()) ; 
      (itSelProcCutVarH1D!=histos1D.end()) ; 
      itSelProcCutVarH1D++) {
    for((itProcCutVarH1D=(itSelProcCutVarH1D->second).begin()) ; 
	(itProcCutVarH1D!=(itSelProcCutVarH1D->second).end()) ; 
	itProcCutVarH1D++) {
      for((itCutVarH1D=(itProcCutVarH1D->second).begin()) ; 
	  (itCutVarH1D!=(itProcCutVarH1D->second).end()) ; 
	  itCutVarH1D++) {
	for((itVarH1D=(itCutVarH1D->second).begin()) ; 
	    (itVarH1D!=(itCutVarH1D->second).end()) ; 
	    itVarH1D++) {
	  var.push_back(itVarH1D->first);
	}
      }
    }
  }

  for((itSelProcCutVarH2D=histos2D.begin()) ; 
      (itSelProcCutVarH2D!=histos2D.end()) ; 
      itSelProcCutVarH2D++) {
    for((itProcCutVarH2D=(itSelProcCutVarH2D->second).begin()) ; 
	(itProcCutVarH2D!=(itSelProcCutVarH2D->second).end()) ; 
	itProcCutVarH2D++) {
      for((itCutVarH2D=(itProcCutVarH2D->second).begin()) ; 
	  (itCutVarH2D!=(itProcCutVarH2D->second).end()) ; 
	  itCutVarH2D++) {
	for((itVarH2D=(itCutVarH2D->second).begin()) ; 
	    (itVarH2D!=(itCutVarH2D->second).end()) ; 
	    itVarH2D++) {
	  var.push_back(itVarH2D->first);
	}
      }
    }
  }

  return var;

}

Int_t GetNGen(TString sample)
{
  
  if(     sample=="qcdht1000to1500") return 4963895;
  else if(sample=="qcdht100to200")   return 81637494;
  else if(sample=="qcdht1500to2000") return 3868886;
  else if(sample=="qcdht2000toinf")  return 1912529;
  else if(sample=="qcdht200to300")   return 18718905;
  else if(sample=="qcdht300to500")   return 19826197;
  else if(sample=="qcdht500to700")   return 19664159;
  else if(sample=="qcdht700to1000")  return 15356448;
  else if(sample=="singletbarw")     return 988500;
  else if(sample=="singletw")        return 995600;
  else if(sample=="singletbart")     return 1680200;
  else if(sample=="singlett")        return 3299800;
  else if(sample=="ttbar")           return 11344206;
  else if(sample=="wln100to200")     return 10152718;
  else if(sample=="wln200to400")     return 5221599;
  else if(sample=="wln400to600")     return 1745914;
  else if(sample=="wln600toinf")     return 1039152;
  else if(sample=="ww")              return 993640;
  else if(sample=="wz")              return 978512;
  else if(sample=="zll100to200")     return 2725655;
  else if(sample=="zll200to400")     return 973937;
  else if(sample=="zll400to600")     return 1067758;
  else if(sample=="zll600toinf")     return 998912;
  else if(sample=="znn100to200")     return 5154824;
  else if(sample=="znn200to400")     return 4998316;
  else if(sample=="znn400to600")     return 1018882;
  else if(sample=="znn600toinf")     return 1008333;
  else if(sample=="zz")              return 996944;

  else return 1;
}
